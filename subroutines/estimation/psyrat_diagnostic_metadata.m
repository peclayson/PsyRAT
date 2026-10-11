function m = psyrat_diagnostic_metadata(table)
%Read the schema record in a convergence table's header; absent is unavailable.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
m=struct('schema','psyrat_diagnostics/1','format','unrecorded', ...
    'precision_significant_digits',NaN,'degraded',true, ...
    'degradation_reason','Diagnostic source/precision not recorded.', ...
    'permitted_constants',{{}},'constant_reasons',{{}}, ...
    'contract_status','unavailable','expected_names',{{}});
if ~iscell(table) || isempty(table); return; end
for k=1:size(table,2)
    v=table{1,k};
    if isstruct(v) && isscalar(v) && isfield(v,'schema') && strcmp(v.schema,m.schema)
        names=fieldnames(v);
        for j=1:numel(names); m.(names{j})=v.(names{j}); end
        return;
    end
end
end
