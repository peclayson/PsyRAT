function fit = psyrat_run_native(chosen,nativespec,varargin)
%Run a native (non-CmdStan) estimation engine and return a PsyRATNativeFit.
%
%fit = psyrat_run_native(chosen,nativespec,varargin)
%
%Inputs:
% chosen - 'hmc' or 'fitlme' (already resolved by psyrat_resolve_engine; this
%  function is never called for 'cmdstan').
% nativespec - struct describing the run for native engines, with fields:
%   .design - struct identifying the model (.key plus any extra descriptor
%             fields the fitter needs beyond the Stan data struct).
%   .opts   - struct of native-engine knobs (bootstrap_reps, hmc_gradient).
% varargin - the SAME name-value list psyrat_run_stan received (minus the
%  'native' spec), so the engine can read 'data', 'seed', 'chains', 'warmup',
%  and 'iter' (= post-warmup sampling) exactly as CmdStan would.
%
%Output:
% fit - a PsyRATNativeFit exposing .block()/.extract()/.convstats(), so the
%  per-case extraction code in psyrat_computevarcomp consumes it unchanged.
%
%See also: psyrat_resolve_engine, psyrat_native_fitlme, psyrat_native_hmc,
% PsyRATNativeFit

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

%pull the inputs the native engines need from the stan-style argument list.
%psyrat_run_stan passes the post-warmup sampling count as 'iter' (MatlabStan's
%name), so map it to nsampling here.
data    = local_nv(varargin,'data',[]);
seed    = local_nv(varargin,'seed',12345);
nchains = local_nv(varargin,'chains',4);
nsamp   = local_nv(varargin,'iter',1000);
nwarm   = local_nv(varargin,'warmup',1000);

design = struct('key','');
if isstruct(nativespec) && isfield(nativespec,'design')
    design = nativespec.design;
end
opts = struct();
if isstruct(nativespec) && isfield(nativespec,'opts')
    opts = nativespec.opts;
end

switch chosen
    case 'fitlme'
        fprintf(['\nEstimating variance components with the native fitlme ' ...
            'engine (REML + parametric bootstrap)...\n']);
        fit = psyrat_native_fitlme(design,data,seed,nchains,opts);
    case 'hmc'
        fprintf(['\nEstimating variance components with the native HMC ' ...
            'engine (hmcSampler)...\n']);
        fit = psyrat_native_hmc(design,data,seed,nchains,nwarm,nsamp,opts);
    otherwise
        error('PsyRAT:badNativeEngine',['psyrat_run_native received an ' ...
            'unexpected engine ''%s'' (expected ''fitlme'' or ''hmc'').'],chosen);
end
end

function v = local_nv(args,name,default)
%Fetch a name-value from a stan-style argument list, with a default.
v = default;
for i = 1:numel(args)-1
    if (ischar(args{i}) || (isstring(args{i}) && isscalar(args{i}))) && ...
            strcmp(char(args{i}),name)
        v = args{i+1};
        return;
    end
end
end
