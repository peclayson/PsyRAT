function [verstr, home, ranok] = psyrat_detect_cmdstan_version()
%Locate the CmdStan install and read its version from `stanc --version`.
%
%[verstr, home, ranok] = psyrat_detect_cmdstan_version()
%
%Background:
% Several PsyRAT components need to know whether (and which) CmdStan is
% installed: the minimum-version preflight (psyrat_require_min_cmdstan) and the
% estimation-engine availability check (psyrat_estimation_engines_available).
% This function centralizes the version READ; since G43 the discovery itself
% goes through psyrat_locate_cmdstan (which delegates to mstan.stan_home, the
% resolver estimation uses), replacing the path-then-env lookup that was
% originally extracted from psyrat_require_min_cmdstan. Env-var-only and
% managed-home installs are therefore found here the same way a fit would
% find them.
%
%Required Inputs:
% No inputs are required.
%
%Outputs:
% verstr - CmdStan version as a 'major.minor.patch' string, or '' when no
%  install is found or the version cannot be read (e.g. stanc fails to run).
% home - the discovered CmdStan home directory, or '' when none was found.
%  NOTE: home can be non-empty while verstr is '' (an install was located but
%  its version could not be read); callers that only need "is CmdStan present"
%  should test home, not verstr.
% ranok - logical; true if the stanc compiler was located AND executed
%  successfully (exit status 0). false means either no install was found or stanc
%  is present but could not be run (wrong-architecture binary, missing shared
%  libraries, no execute permission). A present-but-broken stanc cannot compile a
%  model, so callers should treat ~ranok as "CmdStan unusable" rather than
%  fail-open: distinguishing this from a merely-unreadable version (ranok true,
%  verstr '') lets the engine cascade fall back instead of dying at compile time.

% Copyright (C) 2016-2025 Peter E. Clayson
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

verstr = '';
ranok = false;

%locate the CmdStan home directory through the single toolbox resolver
%(G43): psyrat_locate_cmdstan delegates to mstan.stan_home, the resolver
%estimation actually uses, so this check and a fit see the same install.
%This replaced an inline path-then-env lookup that carried its own (reversed)
%env order and no home-glob tier. Note the semantics that came with it: a
%candidate counts only when it is a VALID CmdStan root (makefile + bin/), so
%a stale env var pointing at a non-CmdStan directory now yields home = ''
%("not found") where the old lookup returned the bogus path with ranok
%false ("present but broken"). Both read as unusable to the callers.
home = psyrat_locate_cmdstan();
if isempty(home)
    return;
end

%path to the stanc compiler (with the .exe extension on Windows)
stanc = fullfile(home,'bin','stanc');
if ispc && exist([stanc '.exe'],'file') == 2
    stanc = [stanc '.exe'];
end

%quote the path so directories with spaces still work
command = ['"' stanc '" --version'];
[status,out] = system(command);
%ranok records whether stanc actually executed (exit status 0). A non-zero exit
%means the binary is present but unusable (wrong arch, missing libs, no exec
%permission); leave ranok false so callers can treat CmdStan as unavailable. A
%zero exit with empty/unparseable output is a usable stanc whose version banner
%just could not be read (ranok true, verstr ''), which stays fail-open.
ranok = (status == 0);
if status ~= 0 || isempty(out)
    return;
end

verstr = regexp(out,'\d+\.\d+\.\d+','match','once');

end
