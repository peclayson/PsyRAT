function RELout = psyrat_checkconv(REL)
%Check for convergence by examining the potential scale reduction factor
%(r_hat) and the effective sample size (n_eff)
%
%[REL, rerun] = psyrat_checkconv(REL)
%
%Lasted Updated 1/19/17
%
%Required Input:
% REL - structure array created by psyrat_computevarcomp
%
%Outputs:
% RELout - structure array with the following added fields
%  out.conv.converged - 1 data converged, 0 data did not converge

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

%nargin is a scalar, so the length(nargin) this guard used to test was always 1
%and the guard never fired. The signature takes exactly one positional input, so
%MATLAB rejects two or more at the call itself; nargin ~= 1 therefore means the
%caller passed none, which previously surfaced as an opaque "Undefined variable
%REL" from the size() call below instead of the message written here.
if nargin ~= 1
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Input not specified \n\n',... 
    'See help psyrat_checkconv for more information on  inputs'));
end

%Only model-backed deterministic rows may be excluded. Unexpected NaNs,
%nonfinite values, malformed rows and missing required diagnostics fail closed.
tables=REL.out.conv.data;
if isempty(tables); tables={};
elseif iscell(tables)&&all(cellfun(@iscell,tables(:))); tables=reshape(tables,1,[]);
else; tables={tables}; end
result=~isempty(tables); excluded={}; invalid={}; diagnostics=cell(size(tables));
for k=1:numel(tables)
    tbl=tables{k}; m=psyrat_diagnostic_metadata(tbl); diagnostics{k}=m;
    if ~iscell(tbl) || size(tbl,1)<2 || size(tbl,2)<3
        result=false; invalid{end+1}='diagnostics_unavailable'; continue; %#ok<AGROW>
    end
    names=tbl(2:end,1); assessed=0;
    validnames=cellfun(@(x)ischar(x)||(isstring(x)&&isscalar(x)),names);
    if ~all(validnames); result=false; end
    names=cellfun(@char,names(validnames),'UniformOutput',false);
    if ~isempty(setdiff(m.expected_names,names))
        result=false; missing=setdiff(m.expected_names,names); invalid=[invalid reshape(missing,1,[])]; %#ok<AGROW>
    end
    for j=2:size(tbl,1)
        name=tbl{j,1}; ess=tbl{j,2}; rhat=tbl{j,3};
        if ~(ischar(name) || (isstring(name)&&isscalar(name))) || ...
                ~isnumeric(ess) || ~isnumeric(rhat) || ~isscalar(ess) || ~isscalar(rhat)
            result=false; invalid{end+1}='malformed_diagnostic_row'; continue; %#ok<AGROW>
        end
        name=char(name);
        if isreal(ess) && isreal(rhat) && isnan(ess) && isnan(rhat) && ...
                any(strcmp(name,m.permitted_constants))
            excluded{end+1}=name; %#ok<AGROW>
            continue;
        end
        assessed=assessed+1;
        if ~isreal(ess) || ~isreal(rhat) || ~isfinite(ess) || ~isfinite(rhat) || ...
                ~(ess>=10*REL.nchains && rhat>0 && rhat<1.1)
            result=false; invalid{end+1}=name; %#ok<AGROW>
        end
    end
    if assessed==0; result=false; end
end
REL.out.conv.converged=double(result);
REL.out.conv.excluded=unique(excluded(:),'stable');
REL.out.conv.invalid=unique(invalid(:),'stable');
REL.out.conv.diagnostics=diagnostics;
REL.out.conv.assessment_scope='early-workflow Rhat/ESS check';
REL.out.conv.thresholds=struct('rhat_less_than',1.1,'ess_at_least',10*REL.nchains);
if isfield(REL,'engine') && strcmpi(REL.engine,'fitlme')
    REL.out.conv.assessment_scope='REML optimizer validity and bootstrap count floor';
    % The Rhat/ESS row is an internal compatibility sentinel, not measured
    % MCMC diagnostics. Describe the actual approved per-fit rule publicly.
    REL.out.conv.thresholds=struct('original_fit_must_be_valid',true, ...
        'per_fit',true,'bootstrap_success_fraction_at_least',.9, ...
        'bootstrap_success_count_at_least',200, ...
        'bootstrap_count_floor_applies_from_requested',200);
end
REL.out.conv.final_inference='unassessed';
REL.out.conv.final_inference_reason='Quantity-specific MCSE and uncertainty calibration are not established by this gate.';
RELout=REL;
end
