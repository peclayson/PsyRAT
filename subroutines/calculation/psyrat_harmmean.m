function hm = psyrat_harmmean(x)
%Compute harmonic mean without requiring toolbox functions.
%
% hm = psyrat_harmmean(x)
%
% Input:
% x - numeric vector of strictly positive values
%
% Output:
% hm - harmonic mean of x

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

x = x(:);
if isempty(x) || any(~isfinite(x)) || any(x <= 0)
    error('varargin:harmmean',... %Error code and associated error
        'Harmonic mean input must be finite values greater than zero.');
end

hm = numel(x) / sum(1 ./ x);

end
