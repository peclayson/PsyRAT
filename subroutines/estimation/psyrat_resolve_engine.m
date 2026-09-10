function chosen = psyrat_resolve_engine(engine_code,nativespec,avail)
%Resolve which estimation engine actually runs, applying the auto cascade,
%capability gating, and clear user-facing warnings/errors.
%
%chosen = psyrat_resolve_engine(engine_code,nativespec)
%chosen = psyrat_resolve_engine(engine_code,nativespec,avail)
%
%Inputs:
% engine_code - 0 = auto (cascade), 1 = cmdstan, 2 = fitlme, 3 = hmc. A forced
%  value (1/2/3) is an explicit user opt-in and overrides the cascade.
% nativespec - struct describing the run for native engines, with at least
%  .design.key (see the psyrat_nativespec local function inside
%  psyrat_computevarcomp). May be [] for a plain CmdStan call.
% avail - (optional) the struct returned by psyrat_estimation_engines_available.
%  Passed in for testability; detected here when omitted.
%
%Output:
% chosen - 'cmdstan' | 'hmc' | 'fitlme'. Errors when no installed engine can
%  handle the request.
%
%Cascade (engine_code == 0): CmdStan -> native HMC -> native fitlme, taking the
%first engine that is both installed and implements this design. When CmdStan is
%present the cascade chooses it silently (the normal, unchanged path). Whenever
%it falls back off CmdStan it emits a clear, prominent warning. fitlme fallbacks
%additionally warn that the intervals are approximate (frequentist), not
%Bayesian credible intervals.
%
%See also: psyrat_estimation_engines_available, psyrat_native_supported

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

if nargin < 3
    avail = [];
end

%design key and which native engines implement it
key = '';
if isstruct(nativespec) && isfield(nativespec,'design') && ...
        isfield(nativespec.design,'key')
    key = char(nativespec.design.key);
end
supported = psyrat_native_supported(key);

switch engine_code
    case 1
        %forced CmdStan: keep today's behavior; the downstream preflight
        %(psyrat_require_min_cmdstan) reports a missing/old CmdStan clearly.
        %No capability detection needed - return before touching `avail` so the
        %common path adds zero overhead.
        chosen = 'cmdstan';

    case 2
        %forced native fitlme (explicit opt-in)
        avail = local_ensure_avail(avail);
        if ~avail.fitlme
            error('PsyRAT:engineUnavailable',['Cannot use the native fitlme ' ...
                'engine: %s\n\nHow to fix: install the Statistics and Machine ' ...
                'Learning Toolbox, or choose a different engine.'], ...
                avail.details.fitlme);
        end
        if ~ismember('fitlme',supported)
            error('PsyRAT:engineUnsupportedDesign',['The native fitlme engine ' ...
                'does not support this analysis (design ''%s''). It handles only ' ...
                'the plain one-facet, test-retest, and splits designs. How to ' ...
                'fix: use CmdStan (engine 1) or, when available, the native HMC ' ...
                'engine (engine 3).'],key);
        end
        local_warn_fitlme_approx(true);
        chosen = 'fitlme';

    case 3
        %forced native HMC (explicit opt-in)
        avail = local_ensure_avail(avail);
        if ~avail.hmc
            error('PsyRAT:engineUnavailable',['Cannot use the native HMC ' ...
                'engine: %s\n\nHow to fix: install the Statistics and Machine ' ...
                'Learning Toolbox, or choose a different engine.'], ...
                avail.details.hmc);
        end
        if ~ismember('hmc',supported)
            error('PsyRAT:engineUnsupportedDesign',['The native HMC engine does ' ...
                'not yet implement this analysis (design ''%s''). How to fix: use ' ...
                'CmdStan (engine 1)%s.'],key,local_fitlme_hint(supported,avail));
        end
        chosen = 'hmc';

    case 0
        %auto cascade: CmdStan -> native HMC -> native fitlme
        avail = local_ensure_avail(avail);
        if avail.cmdstan
            %the normal path: CmdStan present, run unchanged, no warning
            chosen = 'cmdstan';
            return;
        end
        %CmdStan is present but BELOW the minimum version: do not silently fall
        %back to an approximate native engine. Match the forced-CmdStan preflight
        %(psyrat_require_min_cmdstan), which hard-errors on an old CmdStan, and
        %tell the user to upgrade. A native engine is still reachable by forcing
        %it explicitly (engine 2/3), which the message points out.
        if isfield(avail,'cmdstan_status') && strcmp(avail.cmdstan_status,'tooold')
            error('PsyRAT:cmdstanVersion',['%s\n\nHow to fix: upgrade CmdStan ' ...
                '(psyrat_installdependents installs the tested version) and ' ...
                'restart MATLAB. To proceed without upgrading, explicitly choose ' ...
                'a native engine (engine 3 = HMC or engine 2 = fitlme); native ' ...
                'results are approximate / not bit-identical to CmdStan.'], ...
                avail.details.cmdstan);
        end
        %CmdStan absent/broken/unreadable: fall back, warning clearly at each step
        if avail.hmc && ismember('hmc',supported)
            local_warn_banner('PsyRAT:fallbackHMC',sprintf(['CmdStan is ' ...
                'unavailable (%s).\nFalling back to the native HMC engine ' ...
                '(hmcSampler). Results are Bayesian but are NOT bit-identical ' ...
                'to CmdStan.'],avail.details.cmdstan));
            chosen = 'hmc';
        elseif avail.fitlme && ismember('fitlme',supported)
            local_warn_banner('PsyRAT:fallbackFitlme',sprintf(['CmdStan is ' ...
                'unavailable (%s).\nFalling back to the native fitlme engine. ' ...
                'Its intervals are APPROXIMATE (frequentist parametric ' ...
                'bootstrap), NOT Bayesian credible intervals, and may differ ' ...
                'from a CmdStan run.'],avail.details.cmdstan));
            chosen = 'fitlme';
        else
            error('PsyRAT:noEngineAvailable',['No estimation engine can run ' ...
                'this analysis (design ''%s''). CmdStan is unavailable (%s) and ' ...
                'no installed native engine implements this design. How to fix: ' ...
                'install CmdStan (psyrat_installdependents), or install the ' ...
                'Statistics and Machine Learning Toolbox for the designs the ' ...
                'native engines support.'],key,avail.details.cmdstan);
        end

    otherwise
        error('PsyRAT:badEngineCode',['Invalid engine code %g (expected 0 ' ...
            'auto, 1 cmdstan, 2 fitlme, 3 hmc).'],engine_code);
end
end

function avail = local_ensure_avail(avail)
%Detect engine availability on first need (lazy), so the forced-CmdStan path
%returns without ever probing for toolboxes/CmdStan.
if isempty(avail)
    avail = psyrat_estimation_engines_available();
end
end

function local_warn_banner(id,msg)
%Emit a prominent, hard-to-miss fallback notice: a stderr banner (shows even
%when warnings are suppressed) plus a standard MATLAB warning.
bar = repmat('=',1,72);
fprintf(2,'\n%s\n[PsyRAT engine] %s\n%s\n\n',bar,msg,bar);
warning(id,'%s',msg);
end

function local_warn_fitlme_approx(forced)
%Warn that the native fitlme engine's intervals are approximate. 'forced' notes
%that this was an explicit opt-in rather than a cascade fallback.
if forced
    lead = 'Using the native fitlme engine (explicit opt-in).';
else
    lead = 'Using the native fitlme engine.';
end
msg = [lead ' Its intervals are APPROXIMATE (frequentist parametric ' ...
    'bootstrap), NOT Bayesian credible intervals, and may differ from a ' ...
    'CmdStan run.'];
local_warn_banner('PsyRAT:fitlmeApprox',msg);
end

function s = local_fitlme_hint(supported,avail)
%Suggest the fitlme engine in an HMC-unsupported error when fitlme could in
%fact handle the design on this machine.
if ismember('fitlme',supported) && avail.fitlme
    s = ', or the native fitlme engine (engine 2)';
else
    s = '';
end
end
