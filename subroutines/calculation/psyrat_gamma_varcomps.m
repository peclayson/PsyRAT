function vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu, s_v, cov_pv)
%PSYRAT_GAMMA_VARCOMPS One-facet observed-score-scale variance components for
% the scaled-chi-square / Gamma(log link) family.
%
%   vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu)
%   vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu, s_v, cov_pv)
%
% The Gamma family models the single-trial score with a LOG-LINKED mean,
%   log mu_{pit} = alpha + u_p + u_t,   u_p ~ N(0, s_p^2),  u_t ~ N(0, s_i^2),
% and a scaled-chi-square (equivalently Gamma) observation with a person-varying
% dispersion,
%   log nu_p = log(nu) + v_p,   v_p ~ N(0, s_v^2),   Cov(u_p, v_p) = cov_pv,
%   E(Y | mu, nu_p) = mu,   Var(Y | mu, nu_p) = 2*mu^2/nu_p.
%
% WHAT nu IS, PSYCHOMETRICALLY. Var = 2*mu^2/nu makes this a CONSTANT
% COEFFICIENT-OF-VARIATION (proportional-error) model:
%   CV(Y | mu, nu) = sqrt(2/nu),   log CV = 0.5*log(2) - 0.5*log(nu),
% which is 1/3 at the default nu = 18 (psyrat_default_priors, grep nu_center).
% nu is therefore NOT a nuisance shape parameter -- it is the whole error model,
% expressed in RELATIVE units. Two consequences follow, and both matter:
%   (1) A log-nu contrast is a statement about RELATIVE dispersion, not about
%       absolute precision. If b_nu is a log-nu slope the implied log-CV slope is
%       -b_nu/2, so a POSITIVE b_nu means DECREASING relative error. Do not call
%       b_nu "the CV slope" without the sign and the factor, and do not confuse
%       -b_nu/2 with b_sigma: the identity is b_nu = 2*b - 2*b_sigma, i.e.
%       b_sigma = b - b_nu/2, the log RESIDUAL-SD slope, which includes the mean
%       slope b. The two differ by exactly b.
%   (2) Equivalently log sigma = log mu + 0.5*log(2) - 0.5*log(nu), so between
%       any two levels d(log nu) = 2*d(log mu) - 2*d(log sigma): EVERY log-nu
%       contrast blends amplitude with absolute residual SD, whether the contrast
%       is indexed by a covariate, by an event/cell, or by a person.
% This model is not invalid -- it is a coherent constant-CV model. It simply does
% not support an absolute residual-precision reading. PSYRAT_GAMMA_VARCOMPS_LS is
% the direct-residual-SD alternative; see
% documentation/gamma_scale_submodel_decision.md.
%
% Because mu is on the log scale, the Stan/brms group-level SDs (s_p, s_i) are
% on the LOG expected-score scale and are NOT the G-theory variance components.
% This function maps them to OBSERVED-score-scale variance components using
% lognormal moments.
%
% ESTIMAND. The observation-level term is the MARGINAL expected conditional
% variance E[2*mu^2/nu_p], averaged over the population of person dispersions
% (estimand #2), so the components sum to the true observed-score variance by the
% law of total variance - the same convention the Gaussian reliability path uses.
% Passing s_v = cov_pv = 0 (the default) recovers the population-nu plug-in
% (estimand #1), which matches the reference single_cell_components_from_params
% only when dispersion is homogeneous.
%
% Inputs (scalars or equal-length column vectors of posterior draws; scalars
% broadcast against vectors):
%   alpha  - log expected-score-scale grand intercept (b_Intercept)
%   s_p    - person SD on the log expected-score scale (sd_ID__Intercept)
%   s_i    - trial  SD on the log expected-score scale (sd_TrialNumber__Intercept)
%   nu     - population scaled-chi-square degrees of freedom (exp of the log-nu
%            intercept)
%   s_v    - person log-nu SD (optional; default 0). Person heterogeneity in the
%            dispersion.
%   cov_pv - Cov(u_p, v_p), the person mean<->log-nu covariance on the log scale
%            (optional; default 0). Pass the covariance directly, not the
%            correlation, so no division by a possibly-zero cell SD is needed.
%
% Output struct vc, each field a column vector matching the draw length:
%   mu_bar      - marginal expected observed score, E[mu] = exp(alpha + 0.5*ssum)
%   sigma_p2    - observed-score person main-effect variance component
%   sigma_i2    - observed-score trial  main-effect variance component
%   var_mu      - variance of the lognormal mean surface, Var(mu)
%   e_cond_var  - expected conditional (observation) variance, E[Var(Y|mu,nu)]
%   total_var   - total observed variance, Var(mu) + E[Var(Y|mu,nu)]
%   sigma_pi_e2 - residual = person x trial interaction + observation dispersion,
%                 recovered by subtraction: total_var - sigma_p2 - sigma_i2
%
% The residual sigma_pi_e2 deliberately confounds the person x trial interaction
% with the mean-linked observation-level dispersion, exactly as a one-facet
% person x trial design's residual does. Downstream, the observed-score
% components are treated as a standard random-effects ANOVA decomposition and
% projected across trials with the usual linear 1/n' D-study law.
%
% See also PSYRAT_GAMMA_CROSSCOV, PSYRAT_REL_SING.

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

% Work in column-vector form so scalars and draw vectors mix cleanly.
if nargin < 5 || isempty(s_v),    s_v    = 0; end
if nargin < 6 || isempty(cov_pv), cov_pv = 0; end
alpha  = alpha(:);
s_p    = s_p(:);
s_i    = s_i(:);
nu     = nu(:);
s_v    = s_v(:);
cov_pv = cov_pv(:);

% Squared log-scale SDs and their sum (variance of log mu at a single cell).
sp2  = s_p.^2;
si2  = s_i.^2;
ssum = sp2 + si2;

% Marginal expected observed score: E[mu] for a lognormal mean surface.
mu_bar = exp(alpha + 0.5.*ssum);

% Observed-score person main effect: Var_p( E_t[mu | u_p] ). The inner
% expectation E_t[mu|u_p] = exp(alpha + u_p + si2/2) is lognormal in u_p, whose
% variance is exp(2*alpha + si2) * exp(sp2) * (exp(sp2) - 1).
sigma_p2 = exp(2.*alpha + si2) .* exp(sp2) .* (exp(sp2) - 1);

% Observed-score trial main effect: the symmetric expression in the trial SD.
sigma_i2 = exp(2.*alpha + sp2) .* exp(si2) .* (exp(si2) - 1);

% Variance of the full lognormal mean surface mu = exp(alpha + u_p + u_t):
% Var(mu) = E[mu^2] - E[mu]^2 = exp(2*alpha)*(exp(2*ssum) - exp(ssum)).
var_mu = exp(2.*alpha) .* (exp(2.*ssum) - exp(ssum));

% Expected conditional (observation-level) variance, MARGINALIZED over the
% person population's dispersion heterogeneity (estimand #2):
%   E[Var(Y|mu,nu_p)] = E[2*mu^2/nu_p]
%     = (2/nu)*exp(2*alpha + 2*ssum) * exp(0.5*s_v^2 - 2*cov_pv).
% The trailing factor follows from the MGF of the jointly-normal exponent
% 2*u_p + 2*u_t - v_p (only u_p covaries with v_p); it is 1 when s_v = 0,
% recovering the population-nu plug-in. It depends only on the person mean<->log-
% nu coupling (s_v, cov_pv), never on the trial facet - the identical factor
% applies in the multi-facet case (see PSYRAT_GAMMA_VARCOMPS_TRT).
e_cond_var = (2./nu) .* exp(2.*alpha + 2.*ssum) .* exp(0.5.*s_v.^2 - 2.*cov_pv);

% Total observed variance by the law of total variance.
total_var = var_mu + e_cond_var;

% Highest-order residual component (P x trial interaction + scaled-chi-square
% observation dispersion), recovered by subtraction. Clip tiny negative values
% arising only from floating-point noise; a genuinely negative result would
% indicate a specification error upstream.
sigma_pi_e2 = max(total_var - sigma_p2 - sigma_i2, 0);

vc = struct( ...
    'mu_bar',      mu_bar, ...
    'sigma_p2',    sigma_p2, ...
    'sigma_i2',    sigma_i2, ...
    'var_mu',      var_mu, ...
    'e_cond_var',  e_cond_var, ...
    'total_var',   total_var, ...
    'sigma_pi_e2', sigma_pi_e2);

end
