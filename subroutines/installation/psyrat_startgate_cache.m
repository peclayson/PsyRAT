function varargout = psyrat_startgate_cache(action, varargin)
%Session store for the launcher's dependency and release verdicts (G91).
%
% [hit, cached, key, why] = psyrat_startgate_cache('lookup')
% psyrat_startgate_cache('store', key, depvercheck, mincmdstan_ok, ...
%     psyratvercheck, release_status)
% psyrat_startgate_cache('clear')
% rec = psyrat_startgate_cache('get')
%
%Background:
% psyrat_start checks the dependents on every launch: psyrat_checkversionsofdeps
% spawns `stanc --version` (behind a whole-second processManager wait, and a
% second time when the captured output is empty), psyrat_require_min_cmdstan
% spawns it again, and psyrat_checkgithubreleases asks the GitHub releases API
% with an 8 s timeout. The launcher is relaunched from 23 bare `psyrat_start;`
% sites (the viewer Back buttons, cancelled file dialogs, the installer and
% updater arms), so every return to the home screen paid for all of that
% again -- around ten seconds on a firewalled machine (user reports,
% 2026-09-15). This function keeps the three verdicts for the rest of the
% MATLAB session so a relaunch can reuse them.
%
% WHERE THE STORE LIVES, AND WHY HERE. The store is a `persistent` variable
% in THIS function only. It is not inside psyrat_start (its install/update
% callbacks are local functions of that file, so the launcher is on the call
% stack twice over while the installer runs and re-enters it; a `clear
% psyrat_start` from there would be a silent no-op), and it is not inside any
% of the functions whose results it keeps. Those leaf functions must stay
% stateless: TestCmdstanDiscovery and TestInstallationAndWorkflowValidation
% call psyrat_locate_cmdstan, psyrat_detect_cmdstan_version,
% psyrat_checkversionsofdeps and psyrat_checkgithubreleases live, in one
% MATLAB process, under changed environment variables, path shadows and
% injected mocks, and expect a fresh answer every time; the updater
% (psyrat_updatepsyrat) and estimation (psyrat_run_stan's minimum-version
% preflight and psyrat_find_cmdstan_home) consult them live as well. A memo
% inside any of them would make the suite order-dependent and could feed a
% fit a stale install. TestStartgateCache pins that they carry no memo.
%
% WHY `persistent` AND NOT ROOT APPDATA (owner ruling 2026-09-15). Every way
% a persistent store can be shortened produces a RECHECK, never a stale
% answer: a MATLAB restart, `clear all`, `clear functions`, `clear
% psyrat_startgate_cache`, or MATLAB reloading this file after it changes on
% disk. Root appdata would survive `clear all`, so a habitual fresh start
% would no longer recheck and a stale verdict would have no reset short of
% restarting MATLAB. The idiom is the toolbox's own precedent
% (psyrat_gauss_hermite's node table). The one deviation from "recheck only
% after a restart" is therefore that `clear all` also rechecks, which is the
% safe direction and costs one launch.
%
% THE FINGERPRINT (validity key) is what makes reuse self-invalidating. It is
% built on every lookup from CHEAP reads only -- no process spawn, no
% network, no MATLAB-path walk beyond what the launcher's own gate already
% does -- and names the things the three checks actually read:
%   toolboxroot  fileparts(which('psyrat_start')): the checkout in use.
%                psyrat_updatepsyrat rmpaths the old toolbox and addpaths the
%                new one before relaunching, and a second checkout or a
%                worktree copy first on the path changes it too.
%   psyratver    psyrat_defineversion: a newer release relaunched in-session.
%   cmdstanhome  psyrat_locate_cmdstan(): the SAME live resolver the gate,
%                the version check and a fit consult, so it changes when the
%                resolver's ANSWER changes: an environment variable naming a
%                different VALID CmdStan root, a newer managed install
%                (newest-first globs) or a removed install (''). A variable
%                pointing at a folder that is not a CmdStan root is skipped
%                and the resolver falls through to the next candidate
%                (+mstan/stan_home.m), so that change alone leaves this
%                field, and the reuse, unchanged (seen live 2026-09-15).
%   stanc        size and modification time of <home>/bin/stanc (.exe on
%                Windows, as psyrat_detect_cmdstan_version resolves it):
%                catches an in-place rebuild or a same-version reinstall
%                into the same directory, which no path-based signal sees.
%   matlabstan   which('StanModel.m'), the anchor psyrat_matlabstan_bundleinfo
%                reads, and
%   stamp        the size/time of PSYRAT_PATCHLEVEL.txt beside it: the
%                installer re-copying the bundle into a same-named folder.
%   resolver     which('mstan.stan_home'): an unpatched upstream MatlabStan
%                shadowing the bundled one would change discovery itself.
%   mpm          which('processManager.m') and
%   mpmstamp     its size/time: the file psyrat_checkversionsofdeps reads
%                for the MatlabProcessManager version.
% Deliberately NOT in the key: the whole MATLAB path string (the import GUI
% and any user addpath would then invalidate on every launch, turning the
% cache off in practice) and the raw environment variable values (redundant
% with the resolved home). If building the key throws for any reason the key
% is [] and the lookup misses; 'store' refuses an empty key, so a broken
% fingerprint degrades to "always recheck", never to a stale hit. The key
% draws nothing from the random number generator.
%
% A LOOKUP THAT MISMATCHES DROPS THE STORE before returning, so a stale
% verdict can never outlive a key change even if the fresh check that follows
% throws. A failed dependency gate is never stored at all: psyrat_start
% returns from its not-found branch before its only lookup, and its only
% store sits after every live check has returned, before the home screen is
% built. A launch whose gate fails leaves any record from an earlier,
% successful launch in place: nothing consults the store until the gate
% passes again, and then the fingerprint decides whether that record still
% describes the install that was found. Whether a verdict is stored does not depend on its value -- an
% "unable to check" release outcome (-1) is reused like any other (owner
% ruling 2026-09-15), because that is exactly the offline user who paid the
% timeout on every Back click, and no update prompt can fire from it.
%
% NEVER call this function from estimation code: a fit must see the live
% install. TestStartgateCache pins that nothing under subroutines/estimation
% names it.
%
%Actions:
% 'lookup' - build the fingerprint and compare it with the stored one.
%   hit    - true when a record exists and its key equals the fingerprint
%   cached - the record on a hit, [] otherwise
%   key    - the fingerprint just built ([] if it could not be built); pass
%            it back to 'store' so what is stored is the state observed
%            BEFORE the live checks ran
%   why    - '' on a hit; 'no earlier check this MATLAB session' when the
%            store is empty; otherwise a short phrase naming the first key
%            field that differs ('the CmdStan location changed', ...), which
%            the launcher prints so a recheck is never silent
% 'store' - keep the fresh verdicts: the key from 'lookup', the depvercheck
%   struct from psyrat_checkversionsofdeps, the logical from
%   psyrat_require_min_cmdstan, the release statusCode (or the
%   psyrat_updatecheck value) and the release status label. Refused, with no
%   error, when the key is empty.
% 'clear' - forget the record. Called by the launcher's install/update
%   callbacks before the installer or updater runs, and the documented way
%   for a user to force a recheck without restarting MATLAB
%   (psyrat_startgate_cache('clear') and then psyrat_start).
% 'get'   - the whole record for inspection ([] when empty): key,
%   depvercheck, mincmdstan_ok, psyratvercheck, release_status and
%   checked_at (HH:mm:ss of the fresh check).

% Copyright (C) 2016-2026 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see
%     <http://www.gnu.org/licenses/>.
%

%the one session store: [] or a struct with the fields listed under 'get'
persistent store

if nargin < 1 || ~(ischar(action) || (isstring(action) && isscalar(action)))
    error('psyrat_startgate_cache:badAction', ...
        ['psyrat_startgate_cache needs an action. How to fix: call it with ' ...
        '''lookup'', ''store'', ''clear'' or ''get''.']);
end
action = lower(char(action));

switch action
    case 'lookup'
        key = local_fingerprint();
        hit = false;
        cached = [];
        if isempty(store)
            why = 'no earlier check this MATLAB session';
        elseif isempty(key)
            %a fingerprint that could not be built cannot vouch for anything
            store = [];
            why = 'the fingerprint could not be built';
        else
            why = local_first_difference(store.key, key);
            if isempty(why)
                hit = true;
                cached = store;
            else
                %drop the stale record NOW, before the caller reruns the
                %checks: if one of them throws, nothing stale survives
                store = [];
            end
        end
        varargout = {hit, cached, key, why};

    case 'store'
        if numel(varargin) ~= 5
            error('psyrat_startgate_cache:badInput', ...
                ['''store'' takes exactly five inputs. How to fix: pass the ' ...
                'key from ''lookup'', the depvercheck struct, the ' ...
                'minimum-version logical, the release status code and the ' ...
                'release status label.']);
        end
        [key, depvercheck, mincmdstan_ok, psyratvercheck, release_status] = ...
            varargin{:};
        if ~isstruct(depvercheck) || ~(islogical(mincmdstan_ok) || ...
                isnumeric(mincmdstan_ok)) || ~(ischar(release_status) || ...
                (isstring(release_status) && isscalar(release_status)))
            error('psyrat_startgate_cache:badInput', ...
                ['''store'' inputs have the wrong types. How to fix: pass the ' ...
                'depvercheck struct, a logical, the status code and a char label.']);
        end
        %a new store always replaces the old one; an empty key means the
        %fingerprint failed, and nothing may be stored against that
        store = [];
        if isempty(key) || ~isstruct(key)
            return;
        end
        store = struct( ...
            'key', key, ...
            'depvercheck', depvercheck, ...
            'mincmdstan_ok', logical(mincmdstan_ok), ...
            'psyratvercheck', psyratvercheck, ...
            'release_status', char(release_status), ...
            'checked_at', char(datetime('now', 'Format', 'HH:mm:ss')));

    case 'clear'
        store = [];

    case 'get'
        varargout = {store};

    otherwise
        error('psyrat_startgate_cache:badAction', ...
            ['Unknown action ''%s''. How to fix: use ''lookup'', ''store'', ' ...
            '''clear'' or ''get''.'], action);
end

end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function key = local_fingerprint()
%LOCAL_FINGERPRINT The validity key described in the header, or [] on error.
%   Cheap reads only. The try/catch is the fail-safe: an error here must
%   disable reuse (a miss every time), never enable a stale hit.

key = [];
try
    k = struct();
    k.toolboxroot = fileparts(which('psyrat_start'));
    k.psyratver = psyrat_defineversion;
    k.cmdstanhome = psyrat_locate_cmdstan();
    stanc = '';
    if ~isempty(k.cmdstanhome)
        %the compiler the version readers run, resolved the way
        %psyrat_detect_cmdstan_version resolves it (.exe on Windows)
        stanc = fullfile(k.cmdstanhome, 'bin', 'stanc');
        if ispc && exist([stanc '.exe'], 'file') == 2
            stanc = [stanc '.exe'];
        end
    end
    k.stanc = local_filestamp(stanc);
    k.matlabstan = which('StanModel.m');
    stamp = '';
    if ~isempty(k.matlabstan)
        stamp = fullfile(fileparts(k.matlabstan), 'PSYRAT_PATCHLEVEL.txt');
    end
    k.stamp = local_filestamp(stamp);
    k.resolver = which('mstan.stan_home');
    k.mpm = which('processManager.m');
    k.mpmstamp = local_filestamp(k.mpm);
    key = k;
catch
    key = [];
end

end

function s = local_filestamp(f)
%LOCAL_FILESTAMP Size and modification time of one file as a char, or
%   'missing' when there is no such file. Both fields, so an in-place rewrite
%   is caught even when the file system's time resolution is coarse. Never
%   NaN: isequal(NaN, NaN) is false and would invalidate on every launch.

s = 'missing';
if isempty(f)
    return;
end
d = dir(f);
if isempty(d) || d(1).isdir
    return;
end
s = sprintf('%d bytes at %.9f', d(1).bytes, d(1).datenum);

end

function why = local_first_difference(stored, current)
%LOCAL_FIRST_DIFFERENCE '' when every key field matches, otherwise a phrase
%   naming the FIRST field that differs, in the order the header lists them.
%   The launcher prints it, so an automatic recheck is announced in words a
%   user can act on rather than as a bare miss.

fields = { ...
    'toolboxroot', 'the toolbox location changed'; ...
    'psyratver',   'the toolbox version changed'; ...
    'cmdstanhome', 'the CmdStan location changed'; ...
    'stanc',       'the CmdStan compiler changed'; ...
    'matlabstan',  'the MatlabStan copy changed'; ...
    'stamp',       'the MatlabStan patch-level stamp changed'; ...
    'resolver',    'the MatlabStan resolver copy changed'; ...
    'mpm',         'the MatlabProcessManager copy changed'; ...
    'mpmstamp',    'the MatlabProcessManager file changed'};

why = '';
if ~isstruct(stored) || ~isstruct(current)
    why = 'the fingerprint changed';
    return;
end
for ii = 1:size(fields, 1)
    f = fields{ii, 1};
    if ~isfield(stored, f) || ~isfield(current, f) || ...
            ~isequal(stored.(f), current.(f))
        why = fields{ii, 2};
        return;
    end
end

end
