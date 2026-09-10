classdef PsyRATNativeFit
%Duck-typed stand-in for a MatlabStan fit object, returned by the native
%estimation engines (psyrat_run_native) so the rest of psyrat_computevarcomp can
%consume native results through the SAME interface it uses for CmdStan.
%
%Why this exists:
% Every analysis case in psyrat_computevarcomp pulls posterior draws with the
% identical idiom
%     fit.block();
%     REL.out.<name> = fit.extract('pars','<name>').<name>;
%     REL.out.conv.data = psyrat_storeconv(fit);
% (psyrat_storeconv is a local function inside psyrat_computevarcomp, not a
% file of its own, so it has no standalone help entry.)
% Returning an object that mimics .block()/.extract() (and a convergence cell)
% lets the native engines plug in without editing those per-case extraction
% lines - the minimal-diff requirement for the two fragile cores.
%
%Contract reproduced from MatlabStan:
% - block() is a no-op here (native engines run synchronously; there is no
%   background CmdStan process to wait on). Returned for call-compatibility.
% - extract('pars',NAME) returns a struct with one field per requested
%   parameter, each holding the draw array shaped (nchains*nsampling) x dim,
%   exactly as fit.extract returns for CmdStan. NAME may be a char or cellstr.
% - the convergence cell (see psyrat_storeconv's native branch) is the same
%   {'name','n_eff','r_hat'; ...} layout psyrat_checkconv reads. For the HMC
%   engine these are genuine split-R-hat / n_eff values; for the fitlme engine
%   they encode REML convergence (r_hat = 1 and a large n_eff when the fit
%   converged, r_hat = Inf / n_eff = 0 otherwise) so psyrat_checkconv's
%   thresholds resolve to the correct converged/not decision.
%
%See also: psyrat_run_native, psyrat_checkconv

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

    properties
        %struct mapping each Stan parameter name to its draw array,
        %(nchains*nsampling) x dim. Built by the native engine.
        draws struct = struct();

        %convergence cell in psyrat_storeconv layout:
        %{'name','n_eff','r_hat'; label, neff, rhat; ...}
        conv cell = {'name','n_eff','r_hat'};

        %which engine produced this fit: 'hmc' or 'fitlme'. Used for provenance
        %and to drive engine-aware output labeling.
        engine char = '';

        %free-form provenance/diagnostic metadata (e.g. gradient method, REML
        %convergence flag, sampler tuning info). Not consumed by the pipeline.
        meta struct = struct();
    end

    methods
        function obj = PsyRATNativeFit(draws,conv,engine,meta)
            %Construct a native fit. All arguments optional so the object can be
            %built incrementally by an engine if needed.
            if nargin >= 1 && ~isempty(draws)
                obj.draws = draws;
            end
            if nargin >= 2 && ~isempty(conv)
                obj.conv = conv;
            end
            if nargin >= 3 && ~isempty(engine)
                obj.engine = engine;
            end
            if nargin >= 4 && ~isempty(meta)
                obj.meta = meta;
            end
        end

        function block(~)
            %No-op: native engines run synchronously. Present so the per-case
            %fit.block() call (which waits on the CmdStan process) is a harmless
            %no-op for native fits.
        end

        function s = extract(obj,varargin)
            %Mimic MatlabStan's fit.extract('pars',NAME). Returns a struct with
            %one field per requested parameter holding its draw array.
            pars = local_parse_pars(varargin);
            if isempty(pars)
                %extract() with no 'pars' returns all stored parameters, matching
                %MatlabStan's "return everything" behavior closely enough for
                %internal use.
                pars = fieldnames(obj.draws);
            end
            s = struct();
            for i = 1:numel(pars)
                name = pars{i};
                if ~isfield(obj.draws,name)
                    error('PsyRATNativeFit:missingPar', ...
                        ['The native engine did not produce parameter ''%s''. ' ...
                        'This is an internal error: the engine for this design ' ...
                        'must populate every parameter the case extracts.'],name);
                end
                s.(name) = obj.draws.(name);
            end
        end

        function c = convstats(obj,~)
            %Return the convergence cell. The optional second argument (trt) is
            %accepted for signature-compatibility with psyrat_storeconv but is
            %unused: native engines pre-resolve which parameters to monitor.
            c = obj.conv;
        end
    end
end

function pars = local_parse_pars(args)
%Pull the value of a 'pars' name-value out of an extract() argument list and
%normalize it to a cellstr. Returns {} if 'pars' was not supplied.
pars = {};
for i = 1:numel(args)-1
    if (ischar(args{i}) || (isstring(args{i}) && isscalar(args{i}))) && ...
            strcmpi(char(args{i}),'pars')
        val = args{i+1};
        if ischar(val)
            pars = {val};
        elseif isstring(val)
            pars = cellstr(val);
        elseif iscell(val)
            pars = val;
        end
        return;
    end
end
end
