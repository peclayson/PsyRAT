function schema = psyrat_native_output_schema(key,data)
%Expected returned fields independently transcribed from native constructors.
%Every shape follows the draw dimension. Missing/invalid dimensions are NaN.
%Source inventory: a development-side returned-field inventory (not shipped);
%see documentation/diagnostic_output_contract.md. No runtime dependency on it.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
schema=struct(); K=NaN; NSUB=NaN; C=NaN;
if isfield(data,'KDIM'); K=data.KDIM; end
if isfield(data,'NSUB'); NSUB=data.NSUB; end
if isfield(data,'C'); C=data.C; end
if ismember(key,{'sserr','dynrel','diff_trt','diff_dynrel_trt','diff_trt_rescor','diff_dynrel_trt_rescor'}) && ...
        (~isfield(data,'KDIM')||isempty(data.KDIM)); K=0; end
K=dimension(K,0); NSUB=dimension(NSUB,1); C=dimension(C,1);
switch key
    case {'sserr','dynrel'}
        for name={'Intercept','Intercept_sigma','sig_trl'}; schema.(name{1})=[1]; end
        for name={'b','b_sigma'}; schema.(name{1})=[K]; end
        for name={'chol_corrmat'}; schema.(name{1})=[2 2]; end
        for name={'gro_sds'}; schema.(name{1})=[2]; end
        for name={'ind_bs','ind_sd'}; schema.(name{1})=[NSUB]; end
    case {'trt'}
        for name={'pop_int','sig_err','sig_id','sig_occ','sig_occxid','sig_trl','sig_trlxid','sig_trlxocc'}; schema.(name{1})=[1]; end
    case {'diff'}
        for name={'b','b_sigma','sd_id','sd_trl'}; schema.(name{1})=[2]; end
        for name={'id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
    case {'diff_rescor'}
        for name={'b','b_sigma','sd_id','sd_trl'}; schema.(name{1})=[2]; end
        for name={'err_varcov','id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
    case {'diff_dynrel'}
        for name={'b','b_sigma','sd_id','sd_trl'}; schema.(name{1})=[2]; end
        for name={'b_dim','b_sigma_dim'}; schema.(name{1})=[K]; end
        for name={'id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
    case {'diff_dynrel_sserr'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma','sd_trl'}; schema.(name{1})=[2]; end
        for name={'b_dim','b_sigma_dim'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'dynrel_trt'}
        for name={'Intercept','Intercept_sigma','sig_occ','sig_occxid','sig_trl','sig_trlxid','sig_trlxocc'}; schema.(name{1})=[1]; end
        for name={'b','b_sigma'}; schema.(name{1})=[K]; end
        for name={'chol_corrmat'}; schema.(name{1})=[2 2]; end
        for name={'gro_sds'}; schema.(name{1})=[2]; end
        for name={'ind_sd'}; schema.(name{1})=[NSUB]; end
    case {'dod'}
        for name={'b_cell','b_sigma_cell','err_var_Cell'}; schema.(name{1})=[C]; end
        for name={'cor_ID_Cell','cor_TrialNumber_Cell','id_varcov','trl_varcov'}; schema.(name{1})=[C C]; end
    case {'diff_trt','diff_dynrel_trt'}
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim','b_sigma_dim'}; schema.(name{1})=[K]; end
        for name={'id_varcov','occ_varcov','oid_varcov','tid_varcov','to_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
    case {'diff_dynrel_sserr_trt'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim','b_sigma_dim'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'id_varcov','occ_varcov','oid_varcov','tid_varcov','to_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'diff_dynrel_rescor'}
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'err_varcov','id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
    case {'diff_trt_rescor','diff_dynrel_trt_rescor'}
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'err_varcov','id_varcov','occ_varcov','oid_varcov','tid_varcov','to_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
    case {'diff_dynrel_sserr_rescor'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'err_varcov','id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'diff_dynrel_sserr_trt_rescor'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'err_varcov','id_varcov','occ_varcov','oid_varcov','tid_varcov','to_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'diff_dynrel_sserr_rho'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'err_varcov','id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor','rescor_mu','sd_rescor'}; schema.(name{1})=[1]; end
        for name={'rho'}; schema.(name{1})=[NSUB]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'diff_dynrel_sserr_trt_rho'}
        for name={'L_id'}; schema.(name{1})=[4 4]; end
        for name={'b','b_sigma'}; schema.(name{1})=[2]; end
        for name={'b_dim1','b_dim2','b_sigma_dim1','b_sigma_dim2'}; schema.(name{1})=[K]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'err_varcov','id_varcov','occ_varcov','oid_varcov','tid_varcov','to_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor','rescor_mu','sd_rescor'}; schema.(name{1})=[1]; end
        for name={'rho'}; schema.(name{1})=[NSUB]; end
        for name={'sd_id'}; schema.(name{1})=[4]; end
    case {'diff_sserr'}
        for name={'b','b_sigma','sd_id','sd_trl'}; schema.(name{1})=[2]; end
        for name={'er_var_ss'}; schema.(name{1})=[NSUB 2]; end
        for name={'err_varcov','id_varcov','trl_varcov'}; schema.(name{1})=[2 2]; end
        for name={'rescor'}; schema.(name{1})=[1]; end
        for name={'wp_cov_ss'}; schema.(name{1})=[NSUB]; end
end
end
function x=dimension(x,minimum)
if ~isnumeric(x)||~isscalar(x)||~isreal(x)||~isfinite(x)||x<minimum||x~=fix(x); x=NaN; end
end
