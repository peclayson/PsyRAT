function cov = psyrat_gamma_crosscov_rescor(alpha_f, s_p_f, s_i_f, ...
    alpha_s, s_p_s, s_i_s, cor_p, cor_i, cov_e_obs)
%PSYRAT_GAMMA_CROSSCOV_RESCOR Cross-cell observed-score covariance components for
% the CONCURRENT one-facet gamma difference designs, including the two terms the
% non-concurrent helper does not supply.
%
%   cov = psyrat_gamma_crosscov_rescor(alpha_f, s_p_f, s_i_f, ...
%                                      alpha_s, s_p_s, s_i_s, cor_p, cor_i, ...
%                                      cov_e_obs)
%
% Inputs are the eight arguments of PSYRAT_GAMMA_CROSSCOV, plus:
%   cov_e_obs - the RESIDUAL cross-cell covariance on the observed score scale,
%               supplied by the caller (draws vector or scalar). For the modular
%               cut-copula designs this is the Gauss-Hermite conversion of the
%               copula correlation, rescaled by the two components' expected
%               means; see PSYRAT_GAMMA_COPULA_RESCOV.
%
% Output struct cov:
%   cov_p_obs      - observed-score person cross-cell covariance
%   cov_i_obs      - observed-score trial  cross-cell covariance
%   cov_mean_total - total cross-cell covariance of the two mean surfaces
%   cov_pimean_obs - the MEAN-SURFACE person x trial cross term
%   cov_e_obs      - the residual cross term, as supplied
%   cov_pi_e_obs   - cov_pimean_obs + cov_e_obs, the pooled residual channel
%
% WHY THIS EXISTS RATHER THAN AN EDIT TO PSYRAT_GAMMA_CROSSCOV. That helper
% hard-returns a zero residual cross-covariance, and its header states that both
% halves of its contract are load-bearing for the NON-concurrent designs it
% serves (analyses 12/13 and the static difference family). Changing it would
% move shipped numbers. This file leaves it untouched, calls it for the two
% terms it already computes - so the person and trial cross-covariances are
% bit-identical to the shipped path - and adds only what is new.
%
% THE POOLED RESIDUAL CHANNEL, and why it is pooled. Rast & Clayson define the
% residual as a COMBINED quantity: page 4 describes it as the person-by-task
% interaction together with unsystematic measurement error, and page 7 repeats
% that the condition-specific residual may include a person-by-condition
% interaction. In a Gaussian model there is nothing to pool - with no log link
% the interaction IS the residual. Under the log link the observed-scale
% decomposition splits the interaction out as a mean-surface term, so
% reconstructing the paper's estimand requires adding it back. (Cite that
% coefficient by PAGE: its display is unnumbered.) The diagonal counterpart is
% PSYRAT_GAMMA_VARCOMPS_LS's sigma_pi_e2, which is already the pooled quantity.
%
% NOTE the one-facet subject-level designs 6/8/13/26 pool the DIAGONAL the same
% way since S17-FIX (owner-directed 2026-08-09; the pre-fix conditional-only
% behavior was finding S17, now resolved). Their cross-event slot stays at zero
% per the S14 owner decision; the pooled cross term computed HERE ships only on
% the concurrent modular designs, which estimate the residual coupling.
%
% NUMERICS. The mean-surface cross term is the total minus the two named terms,
% but it is NOT computed that way. Writing P for the shared prefactor and a, b
% for the two log-scale covariances,
%
%   P*(e^(a+b) - 1) - P*(e^a - 1) - P*(e^b - 1) = P*(e^a - 1)*(e^b - 1)
%
% identically. The right-hand form is a product of quantities this function
% already holds, so it is exact where the subtraction would cancel three nearly
% equal numbers - and log-scale covariances are routinely small enough
% (1e-3 and below) for that cancellation to cost most of the significant digits.
%
% See also PSYRAT_GAMMA_CROSSCOV, PSYRAT_GAMMA_COPULA_RESCOV,
% PSYRAT_GAMMA_VARCOMPS_LS.

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

% The two shared terms come from the shipped helper, unchanged, so the concurrent
% and non-concurrent paths cannot drift apart on them.
base = psyrat_gamma_crosscov(alpha_f, s_p_f, s_i_f, ...
    alpha_s, s_p_s, s_i_s, cor_p, cor_i);

alpha_f = alpha_f(:); s_p_f = s_p_f(:); s_i_f = s_i_f(:);
alpha_s = alpha_s(:); s_p_s = s_p_s(:); s_i_s = s_i_s(:);
cov_e_obs = cov_e_obs(:);

% The prefactor is the product of the two cells' marginal expected scores; it is
% the same quantity psyrat_gamma_crosscov forms internally.
EY_f = exp(alpha_f + 0.5.*(s_p_f.^2 + s_i_f.^2));
EY_s = exp(alpha_s + 0.5.*(s_p_s.^2 + s_i_s.^2));
prefactor = EY_f .* EY_s;

cov_p_log = cor_p(:) .* s_p_f .* s_p_s;
cov_i_log = cor_i(:) .* s_i_f .* s_i_s;

% Total cross-covariance of the two lognormal mean surfaces.
cov_mean_total = prefactor .* expm1(cov_p_log + cov_i_log);

% The mean-surface person x trial cross term, by the product identity above.
cov_pimean_obs = prefactor .* expm1(cov_p_log) .* expm1(cov_i_log);

if isscalar(cov_e_obs)
    cov_e_obs = repmat(cov_e_obs, size(cov_pimean_obs));
end
if numel(cov_e_obs) ~= numel(cov_pimean_obs)
    error('psyrat_gamma_crosscov_rescor:size', ...
        ['The residual cross-covariance must be a scalar or match the number '...
        'of draws (%d); got %d.'], numel(cov_pimean_obs), numel(cov_e_obs));
end

cov = struct( ...
    'cov_p_obs', base.cov_p_obs, ...
    'cov_i_obs', base.cov_i_obs, ...
    'cov_mean_total', cov_mean_total, ...
    'cov_pimean_obs', cov_pimean_obs, ...
    'cov_e_obs', cov_e_obs, ...
    'cov_pi_e_obs', cov_pimean_obs + cov_e_obs);

end
