function loglik = psyrat_copula_loglik(rho, sufficient)
%PSYRAT_COPULA_LOGLIK Gaussian-copula log likelihood for a bivariate sample,
% evaluated from SUFFICIENT STATISTICS of the normal scores.
%
%   loglik = psyrat_copula_loglik(rho, sufficient)
%
% Inputs
%   rho        - copula correlation(s) to evaluate (scalar or vector)
%   sufficient - struct from PSYRAT_COPULA_SUFFSTATS, or any struct with fields
%                n, sum_sq_1, sum_sq_2, sum_cross
%
% Output
%   loglik - log likelihood at each rho, same size as rho; -Inf where |rho| >= 1
%
% For paired normal scores (z1_j, z2_j), j = 1..n, with a bivariate normal copula
% of correlation rho, the copula log density summed over observations is
%
%   -n/2 * log(1 - rho^2)
%   - [ rho^2 * sum(z1^2 + z2^2) - 2*rho*sum(z1.*z2) ] / [ 2*(1 - rho^2) ].
%
% Everything the likelihood needs about the data therefore enters through three
% scalars plus the count, which is what makes evaluating this on a dense grid of
% rho cheap regardless of sample size: the cost is independent of n once the
% sufficient statistics are formed. That is the property the modular cut-copula
% stage relies on, since it re-evaluates the posterior of rho on a 20,001-point
% grid for every retained marginal posterior draw.
%
% Note this is the copula log likelihood ONLY - the marginal densities are
% omitted because they do not involve rho and so cancel from the conditional
% posterior. Do not read the returned value as a full joint log likelihood.
%
% Port of gaussian_copula_log_likelihood() in the owner's reference bundle
% (R/copula_stage.R).
%
% See also PSYRAT_COPULA_SUFFSTATS, PSYRAT_COPULA_RHO_GRID, PSYRAT_GAMMA_PIT.

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

required = {'n', 'sum_sq_1', 'sum_sq_2', 'sum_cross'};
missing = required(~isfield(sufficient, required));
if ~isempty(missing)
    error('psyrat_copula_loglik:sufficient', ...
        ['The sufficient-statistic struct is missing the field(s): %s. How to '...
        'fix: build it with psyrat_copula_suffstats.'], strjoin(missing, ', '));
end

rho_in = rho;
rho = rho(:);
loglik = -Inf(size(rho));

% The copula density is defined only on the open interval; the grid caller keeps
% its support strictly inside it, but guard rather than return NaN.
valid = isfinite(rho) & abs(rho) < 1;
if ~any(valid)
    loglik = reshape(loglik, size(rho_in));
    return
end

r = rho(valid);
denominator = 1 - r.^2;
sum_squares = sufficient.sum_sq_1 + sufficient.sum_sq_2;

loglik(valid) = -0.5 * sufficient.n * log(denominator) ...
    - 0.5 * (r.^2 * sum_squares - 2 * r * sufficient.sum_cross) ./ denominator;

loglik = reshape(loglik, size(rho_in));

end
