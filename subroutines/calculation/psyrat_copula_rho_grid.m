function post = psyrat_copula_rho_grid(sufficient, uniform_draw, grid_points, ...
    prior_sd, support)
%PSYRAT_COPULA_RHO_GRID Conditional posterior of the Gaussian-copula correlation
% on a dense one-dimensional grid, plus one draw from it.
%
%   post = psyrat_copula_rho_grid(sufficient, uniform_draw)
%   post = psyrat_copula_rho_grid(sufficient, uniform_draw, grid_points, ...
%                                 prior_sd, support)
%
% Inputs
%   sufficient   - struct from PSYRAT_COPULA_SUFFSTATS (normal-score statistics
%                  for ONE retained marginal posterior draw)
%   uniform_draw - a single uniform(0,1) variate; the returned rho_e is the
%                  inverse-CDF draw at this quantile. Supplied by the caller so
%                  that all random number generation stays with the caller's
%                  seeded stream and this function is deterministic.
%   grid_points  - number of grid points, default 20001 (must be odd, >= 2001)
%   prior_sd     - SD of the normal(0, prior_sd) prior on rho, default 0.5
%   support      - [lower upper] bounds, default [-0.95 0.95]
%
% Output struct post:
%   rho_e          - inverse-CDF draw at uniform_draw (a grid point)
%   mean, sd       - conditional posterior mean and SD
%   median         - conditional posterior median
%   q025, q975     - conditional 2.5th and 97.5th percentiles
%   mode           - grid point maximizing the log posterior
%   boundary_mass  - posterior mass at |rho| >= 0.94, i.e. pressed against the
%                    imposed support. A material value means the bound, not the
%                    data, is determining the answer, so the caller records it.
%
% WHY A GRID RATHER THAN HMC. Conditional on a retained draw of the marginal
% parameters, the copula correlation is the ONLY unknown, so its posterior is
% one dimensional and bounded. On that domain a dense deterministic grid is both
% more accurate and far cheaper than sampling, and it removes Monte Carlo noise
% from the conditional summaries entirely. This is what makes the two-stage
% (cut/modular) construction computationally feasible; see
% documentation/modular_cut_copula_workflow.md.
%
% The posterior is evaluated in logs and shifted by its maximum before
% exponentiating, so the normalization is stable even when the likelihood is
% sharply peaked (with tens of thousands of pairs, the unnormalized log density
% ranges over thousands of nats). Integration uses the trapezoid rule; the equal
% grid spacing cancels in the normalization, so only the halved end weights
% survive.
%
% Port of conditional_rho_posterior() in the owner's reference bundle
% (R/copula_stage.R), including the step-function quantile convention (the
% smallest grid point whose cumulative mass reaches the target).
%
% See also PSYRAT_COPULA_SUFFSTATS, PSYRAT_COPULA_LOGLIK,
% PSYRAT_GAMMA_CUT_COPULA_STAGE.

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

if nargin < 3 || isempty(grid_points), grid_points = 20001; end
if nargin < 4 || isempty(prior_sd), prior_sd = 0.5; end
if nargin < 5 || isempty(support), support = [-0.95 0.95]; end

if ~isscalar(uniform_draw) || ~isfinite(uniform_draw) || ...
        uniform_draw <= 0 || uniform_draw >= 1
    error('psyrat_copula_rho_grid:uniform', ...
        ['The inverse-CDF quantile must be strictly between zero and one. '...
        'How to fix: pass a single uniform(0,1) variate from the caller''s '...
        'seeded random stream.']);
end
% Odd so the grid contains rho = 0 exactly when the support is symmetric, which
% keeps the null value representable rather than straddled.
if ~isscalar(grid_points) || ~isfinite(grid_points) ...
        || grid_points ~= fix(grid_points) ...
        || grid_points < 2001 || mod(grid_points, 2) == 0
    %The finite/integer tests are load-bearing: mod(2002.5,2) is 0.5, so a
    %non-integer count would pass the parity test and linspace would silently
    %truncate to an EVEN grid that no longer contains rho = 0.
    error('psyrat_copula_rho_grid:grid', ...
        'grid_points must be an odd integer scalar of at least 2001.');
end
if ~isscalar(prior_sd) || ~isfinite(prior_sd) || prior_sd <= 0
    error('psyrat_copula_rho_grid:prior', 'prior_sd must be positive and finite.');
end
if numel(support) ~= 2 || ~all(isfinite(support)) || support(1) >= support(2) ...
        || support(1) <= -1 || support(2) >= 1
    error('psyrat_copula_rho_grid:support', ...
        'support must be [lower upper] strictly inside (-1, 1).');
end

grid = linspace(support(1), support(2), grid_points)';

% Log posterior up to a constant: copula log likelihood + normal log prior.
% The prior's normalizing constant is retained so the log density matches the
% reference term for term; it cancels in the normalization either way.
log_prior = -0.5*log(2*pi) - log(prior_sd) - 0.5*(grid./prior_sd).^2;
log_density = psyrat_copula_loglik(grid, sufficient) + log_prior;

if ~any(isfinite(log_density))
    error('psyrat_copula_rho_grid:nonfinite', ...
        ['The conditional posterior of the copula correlation is non-finite '...
        'across its whole grid. This means the supplied sufficient statistics '...
        'are corrupt rather than that the data are uninformative.']);
end

% Max-shift before exponentiating: the raw log density spans thousands of nats.
weights = exp(log_density - max(log_density));
% Trapezoid rule. Equal spacing cancels during normalization, so only the halved
% endpoints remain.
weights([1 end]) = 0.5 * weights([1 end]);
weights = weights ./ sum(weights);

cdf = cumsum(weights);

% Step-function quantiles: the first grid point whose cumulative mass reaches
% the target. Matches the reference convention exactly. The target is clamped
% to the accumulated endpoint because cumsum over 20k normalized weights can
% land ~1e-14 BELOW 1, and a uniform draw in that gap would make find return
% empty - a bare mid-loop assignment failure in the stage, after the fit
% (probability ~1e-15 per draw; 2026-08-11 code review). min() is a no-op
% whenever the target is inside the accumulated range, so every other call is
% bit-identical.
qidx = @(p) find(cdf >= min(p, cdf(end)), 1);

posterior_mean = sum(grid .* weights);
posterior_sd = sqrt(sum((grid - posterior_mean).^2 .* weights));
[~, mode_index] = max(log_density);

post = struct( ...
    'rho_e', grid(qidx(uniform_draw)), ...
    'mean', posterior_mean, ...
    'sd', posterior_sd, ...
    'median', grid(qidx(0.5)), ...
    'q025', grid(qidx(0.025)), ...
    'q975', grid(qidx(0.975)), ...
    'mode', grid(mode_index), ...
    'boundary_mass', sum(weights(abs(grid) >= 0.94)));

end
