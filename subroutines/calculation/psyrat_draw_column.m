function values = psyrat_draw_column(draws,column)
%Read one independent fit's complete numeric vector from a draw carrier.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if iscell(draws)
    values = draws{column};
else
    values = draws(:,column);
end
if ~isnumeric(values) || ~isvector(values)
    error('PsyRAT:drawColumn','Expected one numeric draw vector per fit.');
end
values = values(:);
end
