function out = psyrat_append_draws(out,values)
%Append an independent fit's column, preserving every successful draw.
%Equal lengths retain the historical numeric matrix. Unequal lengths use one
%numeric vector per cell. Rows in different columns do not imply joint draws.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if ~isnumeric(values) || ~isvector(values) || isempty(values)
    error('PsyRAT:drawColumn','Expected a nonempty numeric draw vector.');
end
values = values(:);
if isempty(out)
    out = values;
elseif iscell(out)
    out{end+1} = values;
elseif size(out,1) == numel(values)
    out(:,end+1) = values;
else
    out = num2cell(out,1);
    out{end+1} = values;
end
end
