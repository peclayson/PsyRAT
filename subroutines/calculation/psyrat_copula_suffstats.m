function sufficient = psyrat_copula_suffstats(z1, z2)
%PSYRAT_COPULA_SUFFSTATS Sufficient statistics of paired normal scores for the
% Gaussian-copula likelihood.
%
%   sufficient = psyrat_copula_suffstats(z1, z2)
%
% Inputs
%   z1, z2 - paired normal scores (equal-length vectors), typically the output
%            of PSYRAT_GAMMA_PIT for the two concurrently measured components
%
% Output struct sufficient:
%   n         - number of pairs
%   sum_sq_1  - sum(z1.^2)
%   sum_sq_2  - sum(z2.^2)
%   sum_cross - sum(z1 .* z2)
%
% These four numbers are everything the Gaussian-copula likelihood needs about
% the data (see PSYRAT_COPULA_LOGLIK), so forming them once collapses a
% grid evaluation over n observations into a grid evaluation over four scalars.
%
% Port of copula_sufficient_statistics() in the owner's reference bundle
% (R/copula_stage.R).
%
% See also PSYRAT_COPULA_LOGLIK, PSYRAT_COPULA_RHO_GRID, PSYRAT_GAMMA_PIT.

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

z1 = z1(:);
z2 = z2(:);

if isempty(z1) || numel(z1) ~= numel(z2)
    error('psyrat_copula_suffstats:pairs', ...
        ['The two normal-score vectors must be non-empty and the same length: '...
        'the copula is defined on PAIRED observations. How to fix: check that '...
        'the two components were measured concurrently and share a row order.']);
end
if ~all(isfinite(z1)) || ~all(isfinite(z2))
    error('psyrat_copula_suffstats:nonfinite', ...
        'Normal scores must be finite; check the probability transform upstream.');
end

sufficient = struct( ...
    'n', numel(z1), ...
    'sum_sq_1', sum(z1.^2), ...
    'sum_sq_2', sum(z2.^2), ...
    'sum_cross', sum(z1 .* z2));

end
