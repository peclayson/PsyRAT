function conv = psyrat_diagnostic_table(rows,info,source,data,expected)
%Assemble diagnostics with bounded model-backed constants and source evidence.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if nargin<5; expected=rows(:,1); end
contract=psyrat_constant_contract(source,data,rows(:,1));
info.permitted_constants=contract.names;
info.constant_reasons=contract.reasons;
info.contract_status=contract.status;
info.expected_names=reshape(expected,1,[]);
info.model_source_hash=psyrat_hash_struct(source);
info.model_data_hash=psyrat_hash_struct(data);
info.model_flags=struct();
for name={'estimate_rescor','NSUB','C'}
    if isstruct(data)&&isfield(data,name{1}); info.model_flags.(name{1})=data.(name{1}); end
end
conv=cell(1,8);
conv(1,1:3)={'name','n_eff','r_hat'};
conv{1,4}=NaN;
conv(1,5:7)={'ess_bulk','ess_tail','mcse_mean'};
conv{1,8}=info;
if isempty(rows); rows={'convergence_unavailable',0,Inf,NaN,NaN,NaN}; end
conv(2:size(rows,1)+1,1:6)=rows;
end
