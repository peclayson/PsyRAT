function cov = psyrat_gamma_crosscov(alpha_f, s_p_f, s_i_f, ...
    alpha_s, s_p_s, s_i_s, cor_p, cor_i)
%PSYRAT_GAMMA_CROSSCOV Cross-cell observed-score covariance components for the
% scaled-chi-square / Gamma(log link) family (one-facet difference scores).
%
%   cov = psyrat_gamma_crosscov(alpha_f, s_p_f, s_i_f, ...
%                               alpha_s, s_p_s, s_i_s, cor_p, cor_i)
%
% For a two-cell (0 + Cell) joint model, the two cells' person effects are
% bivariate normal on the log expected-score scale with covariance
% cor_p*s_p_f*s_p_s (and analogously for trials). The observed-score covariance
% of the two cells' marginal means is the covariance of two JOINTLY LOGNORMAL
% variables, obtained from lognormal moments exactly as the reference does
% (extract_joint_results in model_fits_stancode_standata.R):
%
%   cov_p_log = cor_p * s_p_f * s_p_s          % log-scale person covariance
%   cov_i_log = cor_i * s_i_f * s_i_s          % log-scale trial  covariance
%   EY_k      = exp(alpha_k + 0.5*(s_p_k^2 + s_i_k^2))   % marginal E[mu], cell k
%   cov_p_obs = EY_f * EY_s * (exp(cov_p_log) - 1)       % observed person cov
%   cov_i_obs = EY_f * EY_s * (exp(cov_i_log) - 1)       % observed trial  cov
%
% Inputs are scalars or equal-length column vectors of posterior draws (scalars
% broadcast). All SD/alpha inputs are on the log expected-score scale; cor_p and
% cor_i are the cross-cell person and trial correlations on that scale.
%
% Output struct cov:
%   cov_p_obs - observed-score person cross-cell covariance
%   cov_i_obs - observed-score trial  cross-cell covariance
%   cov_e_obs - residual cross-cell covariance, defined 0 for non-concurrently
%               measured events (the reference placeholder). A concurrent
%               (correlated-residual) bivariate model is intentionally out of
%               scope for this family.
%
% These feed the off-diagonals of the person and trial 2x2 (co)variance arrays
% consumed by PSYRAT_DIFFREL. The diagonals come from PSYRAT_GAMMA_VARCOMPS under
% the log-nu parameterization and from PSYRAT_GAMMA_VARCOMPS_LS under
% location-scale (gammascale = 2).
%
% SHARED ACROSS BOTH PARAMETERIZATIONS, deliberately. Every argument above is a
% MEAN-submodel quantity, and the cross-cell covariance is a covariance of two
% conditional MEANS: E[Y | u, w] = mu under either parameterization, because the
% scale parameter governs the second moment of Y given mu, not the first. Combined
% with conditional independence of the residuals across cells (non-concurrent -
% see cov_e_obs above), that makes this function parameterization-invariant, so it
% is reused unchanged rather than duplicated. Both halves are load-bearing: under
% a CONCURRENT variant the arguments would still all be mean-submodel, but
% cov_e_obs = 0 would be wrong. See gamma_scale_submodel_decision.md section I.3.
%
% See also PSYRAT_GAMMA_VARCOMPS, PSYRAT_GAMMA_VARCOMPS_LS, PSYRAT_DIFFREL.

% Copyright (C) 2016-2025 Peter E. Clayson
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

% Column-vector form so scalars and draw vectors mix cleanly.
alpha_f = alpha_f(:); s_p_f = s_p_f(:); s_i_f = s_i_f(:);
alpha_s = alpha_s(:); s_p_s = s_p_s(:); s_i_s = s_i_s(:);
cor_p = cor_p(:); cor_i = cor_i(:);

% Log-scale cross-cell covariances of the person and trial effects.
cov_p_log = cor_p .* s_p_f .* s_p_s;
cov_i_log = cor_i .* s_i_f .* s_i_s;

% Marginal expected observed scores per cell (E[mu] for a lognormal surface).
EY_f = exp(alpha_f + 0.5.*(s_p_f.^2 + s_i_f.^2));
EY_s = exp(alpha_s + 0.5.*(s_p_s.^2 + s_i_s.^2));
ey_prod = EY_f .* EY_s;

% Observed-score cross-cell covariances via the jointly-lognormal identity.
cov_p_obs = ey_prod .* (exp(cov_p_log) - 1);
cov_i_obs = ey_prod .* (exp(cov_i_log) - 1);

% Non-concurrent difference score: residual covariance is 0 by definition.
cov_e_obs = zeros(size(cov_p_obs));

cov = struct( ...
    'cov_p_obs', cov_p_obs, ...
    'cov_i_obs', cov_i_obs, ...
    'cov_e_obs', cov_e_obs);

end
