function out = psyrat_gamma_copula_rescov(nu_1, nu_2, rho, tolerance, node_pairs, strict)
%PSYRAT_GAMMA_COPULA_RESCOV Observed-scale residual covariance implied by a
% Gaussian-copula correlation between two scaled-chi-square residuals, by
% converged Gauss-Hermite quadrature.
%
%   out = psyrat_gamma_copula_rescov(nu_1, nu_2, rho)
%   out = psyrat_gamma_copula_rescov(nu_1, nu_2, rho, tolerance, node_pairs)
%
% Inputs
%   nu_1, nu_2 - scaled-chi-square degrees of freedom of the two components,
%                evaluated at the point the estimand is defined at (see the
%                ESTIMAND note below)
%   rho        - Gaussian-copula correlation, |rho| < 1
%   tolerance  - absolute agreement required between two node counts, default 1e-6
%   node_pairs - cell array of [low high] node counts to escalate through,
%                default {[24 32], [32 48], [48 64]}
%   strict     - true (default): exhausting the node pairs is an ERROR.
%                false: return the largest-node estimate with converged = false
%                and the achieved error, leaving the decision to the caller.
%
% Output struct out:
%   value     - covariance of the two STANDARDIZED (mean 1) components
%   nodes     - node count used (0 when rho is numerically zero)
%   error     - absolute difference between the last two node counts
%   converged - whether `error` met `tolerance`
%
% WHEN TO PASS strict = false, AND WHAT IT DOES NOT LICENSE. A single posterior
% draw can put the two components in wildly different regimes - one residual SD
% far below its mean, the other far above - and there the integrand is genuinely
% hard: one margin is nearly a point mass while the other is extremely
% heavy-tailed. Aborting the whole analysis for a handful of tail draws is the
% wrong failure mode when the run took hours. Non-strict mode exists for that,
% and ONLY for that: it still reports the achieved error and marks the draw, so
% the caller can count and surface them. It never silently substitutes zero,
% which would bias the result toward independence. A HIGH non-converged fraction
% is a statement about the fit, not about the quadrature - investigate the
% margins rather than raising the tolerance.
%
% WHY THIS FUNCTION EXISTS. A Gaussian-copula correlation is NOT the
% observed-score residual correlation when the margins are non-Gaussian. The
% copula couples the two components on the latent normal scale; the quantity
% reliability needs is the covariance of the two residuals on the observed
% score scale. With Gamma margins there is no closed form for that map, which is
% precisely the objection recorded in psyrat_computevarcomp.m (ruling E) against
% enabling location-scale parameterizations for the concurrent designs. This
% function supplies the missing map by direct numerical integration.
%
% Writing S_m for the standardized component m (scaled chi-square with nu_m
% degrees of freedom, rescaled to mean 1, so Var(S_m) = 2/nu_m), and coupling
% the two through a bivariate normal (Z1, Z2) with correlation rho,
%
%   S_m = F_m^-1( Phi(Z_m) ),   Cov(S_1, S_2) = E[S_1 * S_2] - 1,
%
% because E[S_m] = 1 by construction. The expectation is a two-dimensional
% Gaussian integral, evaluated here on a tensor grid of Gauss-Hermite nodes with
% Z2 written as rho*Z1 + sqrt(1 - rho^2)*eps for independent standard normal eps.
% The caller rescales by the two components' expected means to reach the
% observed-score covariance.
%
% ESTIMAND (read before interpreting the result). The degrees of freedom passed
% in fix the point at which the map is evaluated. In the person-specific
% workflow they are computed at each participant's posterior expected component
% means, which makes the resulting covariance a MEAN-MATCHED DIAGNOSTIC quantity,
% not a covariance integrated over every realized location-facet effect. It must
% not be described as a fully facet-marginalized residual covariance. This
% labeling requirement comes from the owner's reference bundle and is carried
% through every output surface.
%
% CONVERGENCE. The integrand is smooth but its tails are stretched by the
% quantile transform, so a fixed node count cannot be trusted across the range of
% nu this family produces. The rule escalates through node pairs and accepts the
% first pair whose two estimates agree to `tolerance`, reporting which pair
% succeeded; exhausting the pairs is an error rather than a silently returned
% approximation. In practice the first pair (24 vs 32) converges for realistic
% nu, so the usual cost is two evaluations.
%
% The result is checked against the Cauchy-Schwarz bound
% sqrt(Var(S_1)*Var(S_2)) = sqrt((2/nu_1)*(2/nu_2)); a violation means the
% quadrature has failed rather than that the covariance is large.
%
% BASE MATLAB ONLY (gammaincinv, erfc). No Statistics Toolbox dependency: the
% Gamma quantile is taken through gammaincinv, and the normal CDF through erfc.
% Both tails are handled explicitly, because Phi(z) saturates at exactly 1 in
% double precision beyond about z = 8.3 while the 64-node rule reaches |z| ~ 14.9;
% inverting a saturated 1 would return Inf.
%
% Port of copula_standard_covariance() in the owner's reference bundle
% (R/reliability.R), which reaches the same accuracy via R's log-scale
% pnorm/qgamma.
%
% See also PSYRAT_GAUSS_HERMITE, PSYRAT_GAMMA_CUT_COPULA_STAGE.

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

if nargin < 4 || isempty(tolerance), tolerance = 1e-6; end
if nargin < 5 || isempty(node_pairs), node_pairs = {[24 32], [32 48], [48 64]}; end
if nargin < 6 || isempty(strict), strict = true; end

if ~isscalar(nu_1) || ~isscalar(nu_2) || ~isfinite(nu_1) || ~isfinite(nu_2) ...
        || nu_1 <= 0 || nu_2 <= 0
    error('psyrat_gamma_copula_rescov:nu', ...
        ['Both degrees of freedom must be finite and positive (received '...
        '%.6g and %.6g). A nonpositive value means the component mean or '...
        'residual SD reaching this point is invalid.'], nu_1, nu_2);
end
if ~isscalar(rho) || ~isfinite(rho) || abs(rho) >= 1
    error('psyrat_gamma_copula_rescov:rho', ...
        'The copula correlation must be finite with magnitude below one.');
end

% Exactly independent residuals: the covariance is zero with no integration, and
% reporting nodes = 0 records that no quadrature was needed.
if abs(rho) < 1e-14
    out = struct('value', 0, 'nodes', 0, 'error', 0, 'converged', true);
    return
end

bound = sqrt((2/nu_1) * (2/nu_2));

% Evaluations are cached across the escalation because consecutive pairs share a
% node count (32 appears in pairs 1 and 2, 48 in pairs 2 and 3).
evaluated = containers.Map('KeyType', 'double', 'ValueType', 'double');

lastError = NaN;
lastHigh = NaN;
lastNodes = NaN;
for p = 1:numel(node_pairs)
    pair = node_pairs{p};
    low = localCovarianceOnce(pair(1), nu_1, nu_2, rho, bound, evaluated);
    high = localCovarianceOnce(pair(2), nu_1, nu_2, rho, bound, evaluated);
    lastError = abs(high - low);
    lastHigh = high;
    lastNodes = pair(2);
    if lastError <= tolerance
        out = struct('value', high, 'nodes', pair(2), 'error', lastError, ...
            'converged', true);
        return
    end
end

if ~strict
    %Best available estimate, explicitly marked as not meeting tolerance. It is
    %still bound-checked (localCovarianceOnce would have thrown otherwise), so
    %this is an imprecise number rather than an impossible one.
    out = struct('value', lastHigh, 'nodes', lastNodes, 'error', lastError, ...
        'converged', false);
    return
end

error('psyrat_gamma_copula_rescov:convergence', ...
    ['Gauss-Hermite quadrature for the copula residual covariance did not '...
    'converge for nu = (%.6g, %.6g), rho = %.6g: the largest node pair still '...
    'disagreed by %.6g against a tolerance of %.6g. How to fix: widen '...
    '''node_pairs'' or loosen ''tolerance'' deliberately, and record that the '...
    'estimand was evaluated at reduced accuracy - do not accept the '...
    'unconverged value silently.'], nu_1, nu_2, rho, lastError, tolerance);

end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function value = localCovarianceOnce(n_nodes, nu_1, nu_2, rho, bound, evaluated)
%LOCALCOVARIANCEONCE Tensor-grid Gauss-Hermite estimate at one node count.

if isKey(evaluated, n_nodes)
    value = evaluated(n_nodes);
    return
end

[nodes, weights] = psyrat_gauss_hermite(n_nodes);

% Tensor grid: z1 varies slowly (each node repeated n times), eps varies fast.
z1 = repelem(nodes, n_nodes);
epsilon = repmat(nodes, n_nodes, 1);
grid_weights = repelem(weights, n_nodes) .* repmat(weights, n_nodes, 1);

% Conditional construction of the second latent normal at correlation rho.
z2 = rho * z1 + sqrt(1 - rho^2) * epsilon;

% S_1 depends only on the n distinct z1 values, so transform the nodes once and
% expand. Elementwise identical to transforming the full n^2 vector.
s1 = repelem(localStandardizedQuantile(nodes, nu_1), n_nodes);
s2 = localStandardizedQuantile(z2, nu_2);

% E[S_1*S_2] - 1, since both standardized components have unit mean.
value = sum(grid_weights .* s1 .* s2) - 1;

if ~isfinite(value) || abs(value) > bound + 1e-8
    error('psyrat_gamma_copula_rescov:bound', ...
        ['The quadrature returned %.8g for nu = (%.6g, %.6g), rho = %.6g, '...
        'which violates the Cauchy-Schwarz bound %.8g. This indicates a '...
        'quadrature failure, not a large covariance.'], ...
        value, nu_1, nu_2, rho, bound);
end

evaluated(n_nodes) = value; %#ok<NASGU>

end

function s = localStandardizedQuantile(z, nu)
%LOCALSTANDARDIZEDQUANTILE Scaled-chi-square quantile at a normal score.
%   Returns F^-1(Phi(z)) for the standardized component: shape = nu/2,
%   rate = nu/2, hence mean 1 and variance 2/nu.
%
%   Each tail is taken from its OWN side so no probability near 1 is ever
%   formed: for z <= 0 the lower tail Phi(z) is small and accurate, and for
%   z > 0 the upper tail Phi(-z) is small and accurate, inverted through
%   gammaincinv's 'upper' branch. Using the lower tail throughout would return
%   Inf for the large nodes, since Phi(z) rounds to exactly 1 past z ~ 8.3.

shape = nu/2;
rate = nu/2;

s = zeros(size(z));
lower_side = z <= 0;

% Phi(z) = 0.5*erfc(-z/sqrt(2)); accurate for the small (lower) tail.
p_lower = 0.5 * erfc(-z(lower_side)/sqrt(2));
s(lower_side) = gammaincinv(p_lower, shape) ./ rate;

% Upper tail: Phi(-z) = 0.5*erfc(z/sqrt(2)), inverted on the upper branch.
q_upper = 0.5 * erfc(z(~lower_side)/sqrt(2));
s(~lower_side) = gammaincinv(q_upper, shape, 'upper') ./ rate;

end
