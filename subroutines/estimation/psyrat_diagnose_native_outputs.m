function fit = psyrat_diagnose_native_outputs(fit,key,data,nsamp,nchains)
%Diagnose returned quantities without changing samples or their chain order.
%The native sampler packs nsamp consecutive rows per chain. Constructor-backed
%constants are dispositions, not empirical zero-variance exemptions.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
c=fit.conv; meta=psyrat_diagnostic_metadata(c); absent={};
schema=psyrat_native_output_schema(key,data); required=fieldnames(schema); badfields={};
if isempty(required); addrow('native_output_contract_unavailable',0,Inf,0,0,NaN); end
for k=1:numel(required)
    name=required{k}; shape=[nsamp*nchains schema.(name)]; valid=isfield(fit.draws,name)&&all(isfinite(shape));
    if valid
        v=fit.draws.(name);
        valid=isnumeric(v)&&numel(v)==prod(shape)&&all(arrayfun(@(d)size(v,d),1:numel(shape))==shape);
    end
    if ~valid
        addrow(['required_shape.' name],0,Inf,0,0,NaN); badfields{end+1}=name;
    end
end
fields=fieldnames(fit.draws);
for f=1:numel(fields)
    name=fields{f}; values=fit.draws.(name);
    append_field(name,values,false);
end
if isfield(fit.draws,'er_var_ss')&&~ismember('er_var_ss',badfields)
    append_field('residual_variance',exp(2*fit.draws.er_var_ss),true);
end
if all(isfield(fit.draws,{'Intercept_sigma','ind_sd'}))&&~any(ismember({'Intercept_sigma','ind_sd'},badfields))
    append_field('residual_variance_at_dimension_zero', ...
        exp(2*(fit.draws.Intercept_sigma+fit.draws.ind_sd)),true);
end
meta.monitor_scope=['selected unconstrained parameters; all nonempty returned ', ...
    'quantities; subject residual variances at the reference dimension'];
meta.unassessed_quantities={'dynamic surfaces away from reference dimension', ...
    'DoD subject/trial conditional residual combinations','G/D/ICC ratios','interval endpoints','trial-count cutoffs'};
meta.absent_fields=absent;
meta.contract_status='native constructor rules';
meta.native_design=key;
meta.required_output_schema=schema;
meta.invalid_output_fields=badfields;
meta.native_source_hash=psyrat_hash_file(which('psyrat_native_hmc'));
meta.chain_layout=struct('order','contiguous chain blocks','iterations',nsamp,'chains',nchains);
for k=1:size(c,2)
    if isstruct(c{1,k})&&isfield(c{1,k},'schema'); c{1,k}=meta; break; end
end
fit.conv=c; fit.meta.diagnostics=meta;

    function append_field(name,values,derived)
        if ~isnumeric(values)||size(values,1)~=nsamp*nchains
            addrow(['returned.' name],0,Inf,0,0,NaN); return;
        end
        if isempty(values)
            if isfield(schema,name)&&isequal(schema.(name),0)&&~ismember(name,badfields)
                absent{end+1}=[name ': zero parameter dimension'];
            else
                addrow(['returned.' name],0,Inf,0,0,NaN);
            end
            return;
        end
        dims=size(values); dims=dims(2:end);
        flat=reshape(values,nsamp*nchains,[]);
        for j=1:size(flat,2)
            if numel(dims)==1; ix=j;
            else
                indices=cell(1,numel(dims)); [indices{:}]=ind2sub(dims,j); ix=cell2mat(indices);
            end
            label=['returned.' name];
            if size(flat,2)>1||numel(dims)>1
                label=[label '[' strjoin(arrayfun(@num2str,ix,'UniformOutput',false),',') ']']; %#ok<AGROW>
            end
            x=reshape(flat(:,j),nsamp,nchains);
            [allowed,value,jitter]=constant_entry(name,ix,derived);
            if allowed && isreal(x) && all(isfinite(x(:))) && ...
                    ((jitter && all(abs(x(:)-1)<=32*(eps/2)/(1-32*(eps/2)))) || (~jitter && all(x(:)==value)))
                addrow(label,NaN,NaN,NaN,NaN,NaN);
                meta.permitted_constants{end+1}=label;
                meta.constant_reasons{end+1}=['native ' key ' constructor structural constant'];
            elseif allowed
                addrow(label,0,Inf,0,0,NaN);
            else
                d=psyrat_rank_diagnostics(x);
                e=min(d.ess_bulk,d.ess_tail);
                if isnan(d.ess_bulk)||isnan(d.ess_tail); e=NaN; end
                addrow(label,e,d.rhat,d.ess_bulk,d.ess_tail,d.mcse_mean);
            end
        end
    end
    function addrow(name,e,r,b,t,m)
        c(end+1,1:6)={name,e,r,b,t,m};
        meta.expected_names{end+1}=name;
    end
    function [allowed,value,jitter]=constant_entry(name,ix,derived)
        allowed=false; value=NaN; jitter=false;
        if derived; return; end
        if numel(ix)==2
            i=ix(1); j=ix(2);
            loc=ismember(key,{'sserr','dynrel','dynrel_trt'})&&strcmp(name,'chol_corrmat')&&all(ix<=2);
            subj=ismember(key,{'diff_dynrel_sserr','diff_dynrel_sserr_trt', ...
                'diff_dynrel_sserr_rescor','diff_dynrel_sserr_trt_rescor', ...
                'diff_dynrel_sserr_rho','diff_dynrel_sserr_trt_rho'})&&strcmp(name,'L_id')&&all(ix<=4);
            if (loc||subj)&&(i<j||(i==1&&j==1)); allowed=true; value=double(i==j); end
            if strcmp(key,'dod')&&isfield(data,'C')&&isequal(double(data.C),4)&& ...
                    ismember(name,{'cor_ID_Cell','cor_TrialNumber_Cell'})&&i==j&&i<=4
                %L*L' is mathematically unit diagonal; its roundoff is not a
                %stochastic parameter. Preserve raw values; assert the independently
                %derived gamma_32 forward-error bound on their unit value.
                %See review-science/DOD_DIAGONAL_BOUND.md.
                allowed=true; value=1; jitter=true;
            end
        end
        if strcmp(key,'diff_sserr')&&isfield(data,'estimate_rescor')&&isequal(double(data.estimate_rescor),0)
            if ismember(name,{'rescor','wp_cov_ss'})|| ...
                    (strcmp(name,'err_varcov')&&numel(ix)==2&&ix(1)~=ix(2)&&all(ix<=2))
                allowed=true; value=0;
            end
        end
    end
end
