function pit = psyrat_gamma_pit(y, mu, sigma, clip)
%PSYRAT_GAMMA_PIT Numerically stable probability integral transform (PIT) and
% normal scores for the scaled-chi-square / Gamma family in its MEAN-SD
% (location-scale) parameterization.
%
%   pit = psyrat_gamma_pit(y, mu, sigma)
%   pit = psyrat_gamma_pit(y, mu, sigma, clip)
%
% Inputs
%   y     - observed positive scores (vector)
%   mu    - conditional mean of each observation (vector or scalar)
%   sigma - conditional residual SD of each observation (vector or scalar)
%   clip  - tail probability floor, default 1e-12 (matches the reference)
%
% Output struct pit:
%   z             - normal scores, Phi^-1(F(y)), same size/order as y
%   lower_clipped - count of observations whose LOWER tail hit the floor
%   upper_clipped - count of observations whose UPPER tail hit the floor
%
% PARAMETERIZATION. shape = (mu/sigma)^2 and rate = mu/sigma^2, so that
% E(Y) = shape/rate = mu and Var(Y) = shape/rate^2 = sigma^2 exactly. This is
% the same moment matching the location-scale Stan models use (nu = 2*(mu/sigma)^2
% for the scaled chi-square); see documentation/gamma_scale_submodel_decision.md.
%
% WHY THE TWO-TAIL FORM IS REQUIRED, NOT AN OPTIMIZATION. The normal score of an
% observation deep in the UPPER tail cannot be recovered from the lower-tail CDF:
% for z beyond about +8.3, Phi(z) rounds to exactly 1 in double precision, and
% inverting that returns Inf. The 48- and 64-node quadrature rules this family
% uses reach |z| ~ 14.9, well past that point. This function therefore always
% works with the SMALLER of the two tails, which is computed to full relative
% accuracy by gammainc's 'upper' option, and reflects the sign afterwards using
% Phi^-1(1-q) = -Phi^-1(q). The quantity 1 - (small number) is never formed.
%
% The reference implementation (stable_gamma_pit in the owner's bundle,
% R/copula_stage.R) reaches the same accuracy through R's log-scale pgamma/qnorm
% (log.p = TRUE); R's -expm1(log_upper) and the reflection used here agree to
% within double-precision rounding, and the ported unit tests check that against
% frozen reference values.
%
% CLIPPING. Tail probabilities below `clip` are raised to `clip`, and the number
% of observations affected is reported per tail. Clipping is a numerical floor,
% not a model choice: a high clip fraction means the fitted margins do not
% describe the extremes of the data, which is why the caller records it as a
% diagnostic rather than silently proceeding.
%
% BASE MATLAB ONLY (gammainc, erfcinv). No Statistics Toolbox dependency.
%
% See also PSYRAT_GAMMA_CUT_COPULA_STAGE, PSYRAT_COPULA_LOGLIK.

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

if nargin < 4 || isempty(clip)
    clip = 1e-12;
end
%The floor is applied to the SMALLER tail before the sign reflection below, so
%a clip at or above 0.5 would replace every normal score with +/-Phi^-1(clip)
%and silently REVERSE the transform's ordering - a wrong-sign rho_e with no
%error. Reject rather than trust the caller.
if ~isscalar(clip) || ~isfinite(clip) || clip <= 0 || clip >= 0.5
    error('psyrat_gamma_pit:clip',...
        ['The probability floor ''clip'' must be a finite scalar strictly '...
        'between 0 and 0.5 (default 1e-12).']);
end

y = y(:);
% Scalar mu/sigma broadcast to the length of y so callers can pass either a
% per-observation vector (the estimation path) or a common value (tests).
if isscalar(mu), mu = repmat(mu, size(y)); else, mu = mu(:); end
if isscalar(sigma), sigma = repmat(sigma, size(y)); else, sigma = sigma(:); end

if numel(mu) ~= numel(y) || numel(sigma) ~= numel(y)
    error('psyrat_gamma_pit:size', ...
        ['mu and sigma must be scalars or match the number of observations '...
        'in y. How to fix: pass one mean and one residual SD per observation.']);
end
if ~all(isfinite(y)) || any(y <= 0)
    error('psyrat_gamma_pit:support', ...
        ['The scaled-chi-square family has positive support, so every score '...
        'must be finite and greater than zero.']);
end
if ~all(isfinite(mu)) || any(mu <= 0) || ~all(isfinite(sigma)) || any(sigma <= 0)
    error('psyrat_gamma_pit:parameters', ...
        ['Every conditional mean and residual SD must be finite and positive. '...
        'A nonpositive value here means the reconstructed linear predictor '...
        'overflowed or the wrong draw slot was supplied.']);
end

% Mean-SD to shape-rate. Both are strictly positive given the checks above.
shape = (mu ./ sigma).^2;
rate = mu ./ sigma.^2;

% Regularized incomplete gamma gives both tails of the CDF directly. gammainc
% takes the scaled argument rate*y because it is defined for unit rate.
scaled = rate .* y;
lower = gammainc(scaled, shape);
upper = gammainc(scaled, shape, 'upper');

% Work with the smaller tail throughout: it carries full relative precision,
% and its complement never has to be formed.
use_lower = lower < 0.5;
tail = upper;
tail(use_lower) = lower(use_lower);

% Count clipping per tail BEFORE flooring, so the diagnostic reflects the data.
atfloor = tail < clip;
lower_clipped = sum(atfloor & use_lower);
upper_clipped = sum(atfloor & ~use_lower);
tail = max(tail, clip);

% Normal score from the small tail: Phi^-1(q) = -sqrt(2)*erfcinv(2q) is the
% lower-tail form; the upper tail reflects through Phi^-1(1-q) = -Phi^-1(q).
z = -sqrt(2) * erfcinv(2 * tail);
z(~use_lower) = -z(~use_lower);

if ~all(isfinite(z))
    error('psyrat_gamma_pit:nonfinite', ...
        ['A normal score is non-finite after tail clipping, which should be '...
        'unreachable: the clip floor bounds the transform. This indicates a '...
        'corrupted mean or residual SD reaching the transform.']);
end

pit = struct( ...
    'z', z, ...
    'lower_clipped', lower_clipped, ...
    'upper_clipped', upper_clipped);

end
