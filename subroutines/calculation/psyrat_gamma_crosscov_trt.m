function cov = psyrat_gamma_crosscov_trt(alpha_x, s_p_x, s_o_x, s_t_x, ...
    s_po_x, s_pt_x, s_ot_x, alpha_y, s_p_y, s_o_y, s_t_y, ...
    s_po_y, s_pt_y, s_ot_y, cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot)
%PSYRAT_GAMMA_CROSSCOV_TRT Cross-cell observed-score covariance components for
% the scaled-chi-square / Gamma(log link) family in the crossed Persons x
% Occasions x Trials (test-retest) design (two-facet difference scores).
%
%   cov = psyrat_gamma_crosscov_trt(alpha_x, s_p_x, s_o_x, s_t_x, s_po_x, ...
%             s_pt_x, s_ot_x, alpha_y, s_p_y, s_o_y, s_t_y, s_po_y, s_pt_y, ...
%             s_ot_y, cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot)
%
% This is the two-facet analog of PSYRAT_GAMMA_CROSSCOV (one-facet). For a two-
% cell (0 + Cell) joint model, each grouping factor's two cells are bivariate
% normal on the log expected-score scale with covariance cor_k * s_k_x * s_k_y.
% The observed-score covariance of the two cells' marginal means is the
% covariance of jointly lognormal variables, obtained from lognormal moments
% exactly as the reference implementation does (twofacet_cov_components_from_params
% in model_fits_twofacet_chisq_stancode_standata_v2.R):
%
%   cov_k_log = cor_k * s_k_x * s_k_y                    % log-scale factor-k cov
%   EY_c      = exp(alpha_c + 0.5 * V_c)                 % marginal E[mu], cell c
%   cfun(c)   = EY_x * EY_y * (exp(c) - 1)               % jointly-lognormal cov
% with V_c the sum of that cell's six log-scale variances. Main-effect cross-
% covariances are cfun of the single factor's log-cov; two-way interaction cross-
% covariances follow by inclusion-exclusion (cfun of the two parents' + the
% interaction's log-covs, minus the two parent cross-covs); and the residual
% (three-way + dispersion) cross-covariance is recovered by subtraction from the
% total mean-surface cross-covariance.
%
% All computed cross-covariances are returned. The non-concurrent difference-
% score policy (whether trial-indexed / occasion / residual cross-covariances are
% used) is applied by the CALLER, mirroring the reference's cov_est vs cov_used
% split.
%
% Inputs are scalars or equal-length column vectors of posterior draws (scalars
% broadcast). All alpha/SD inputs are on the log expected-score scale; cor_* are
% the cross-cell correlations of each factor on that scale. Argument naming
% matches PSYRAT_GAMMA_VARCOMPS_TRT: p person, o occasion, t trial, po occasion x
% person, pt trial x person, ot trial x occasion.
%
% Output struct cov, each field a column vector matching the draw length:
%   cov_p, cov_o, cov_t   - observed-score main-effect cross-cell covariances
%   cov_po, cov_pt, cov_ot- observed-score two-way interaction cross-cell covs
%   cov_res               - observed-score residual (three-way + dispersion) cross
%                           covariance, by subtraction from cov_mu_total
%   cov_mu_total          - total mean-surface cross-cell covariance
%
% See also PSYRAT_GAMMA_CROSSCOV, PSYRAT_GAMMA_VARCOMPS_TRT, PSYRAT_DIFFREL_TRT.

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
alpha_x = alpha_x(:); s_p_x = s_p_x(:); s_o_x = s_o_x(:); s_t_x = s_t_x(:);
s_po_x = s_po_x(:); s_pt_x = s_pt_x(:); s_ot_x = s_ot_x(:);
alpha_y = alpha_y(:); s_p_y = s_p_y(:); s_o_y = s_o_y(:); s_t_y = s_t_y(:);
s_po_y = s_po_y(:); s_pt_y = s_pt_y(:); s_ot_y = s_ot_y(:);
cor_p = cor_p(:); cor_o = cor_o(:); cor_t = cor_t(:);
cor_po = cor_po(:); cor_pt = cor_pt(:); cor_ot = cor_ot(:);

% Per-cell total log-scale variance (sum of the six factors), for E[mu].
V_x = s_p_x.^2 + s_o_x.^2 + s_t_x.^2 + s_po_x.^2 + s_pt_x.^2 + s_ot_x.^2;
V_y = s_p_y.^2 + s_o_y.^2 + s_t_y.^2 + s_po_y.^2 + s_pt_y.^2 + s_ot_y.^2;

% Marginal expected observed scores per cell (E[mu] for a lognormal surface).
EY_x = exp(alpha_x + 0.5.*V_x);
EY_y = exp(alpha_y + 0.5.*V_y);
ey_prod = EY_x .* EY_y;

% Log-scale cross-cell covariances per factor.
cov_p_log  = cor_p  .* s_p_x  .* s_p_y;
cov_o_log  = cor_o  .* s_o_x  .* s_o_y;
cov_t_log  = cor_t  .* s_t_x  .* s_t_y;
cov_po_log = cor_po .* s_po_x .* s_po_y;
cov_pt_log = cor_pt .* s_pt_x .* s_pt_y;
cov_ot_log = cor_ot .* s_ot_x .* s_ot_y;

% Jointly-lognormal covariance of two cells whose shared log-effects have summed
% covariance c: cov = E[mu_x]*E[mu_y]*(exp(c) - 1).
cfun = @(c) ey_prod .* (exp(c) - 1);

% Main-effect cross-covariances.
cov_p = cfun(cov_p_log);
cov_o = cfun(cov_o_log);
cov_t = cfun(cov_t_log);

% Two-way interaction cross-covariances by inclusion-exclusion (parents + the
% interaction, minus the two parent cross-covariances).
cov_po = cfun(cov_p_log + cov_o_log + cov_po_log) - cov_p - cov_o;
cov_pt = cfun(cov_p_log + cov_t_log + cov_pt_log) - cov_p - cov_t;
cov_ot = cfun(cov_o_log + cov_t_log + cov_ot_log) - cov_o - cov_t;

% Total mean-surface cross-cell covariance and the residual (three-way P x O x T
% + observation dispersion) cross-covariance by subtraction. The scaled-chi-square
% observations of two DISTINCT cells are independent, so the total mean-surface
% cross-covariance is the full cross-covariance of the observed scores.
cov_mu_total = cfun(cov_p_log + cov_o_log + cov_t_log ...
    + cov_po_log + cov_pt_log + cov_ot_log);
cov_res = cov_mu_total - cov_p - cov_o - cov_t - cov_po - cov_pt - cov_ot;

cov = struct( ...
    'cov_p',        cov_p, ...
    'cov_o',        cov_o, ...
    'cov_t',        cov_t, ...
    'cov_po',       cov_po, ...
    'cov_pt',       cov_pt, ...
    'cov_ot',       cov_ot, ...
    'cov_res',      cov_res, ...
    'cov_mu_total', cov_mu_total);

end
