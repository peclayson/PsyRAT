function vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma, s_sd)
%PSYRAT_GAMMA_VARCOMPS_LS One-facet observed-score-scale variance components for
% the LOCATION-SCALE scaled-chi-square / Gamma(log link) family (gammascale = 2).
%
%   vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma)
%   vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, log_sigma, s_sd)
%
% This is the sibling of PSYRAT_GAMMA_VARCOMPS for the residual-SD
% parameterization. Both describe the same scaled-chi-square likelihood; they
% differ in which quantity carries the linear predictor and the person random
% effect. Here the model is
%   log mu_{pit} = alpha + u_p + u_t,   u_p ~ N(0, s_p^2),  u_t ~ N(0, s_i^2),
%   log sigma_p  = log_sigma + w_p,     w_p ~ N(0, s_sd^2),
%   nu_{pit}     = 2*(mu_{pit}/sigma_p)^2,
%   E(Y | mu, sigma) = mu,   Var(Y | mu, sigma) = sigma_p^2.
%
% WHY THIS EXISTS. Under the log-nu sibling the conditional variance is
% 2*mu^2/nu, which is MEAN-COUPLED: it moves with the trial-specific conditional
% mean, so a dimension slope on log-nu blends a mean-level effect with a
% dispersion effect and cannot be read as a statement about precision alone.
% Here the residual is decoupled from the mean surface entirely, so b_sigma is a
% pure residual-SD slope. That separation is the reason this parameterization
% exists; see documentation/gamma_scale_submodel_decision.md.
%
% The log-nu sibling is NOT invalid -- it is a coherent constant-CV
% (proportional-error) model. The difference is what the scale parameter means:
% log-nu describes RELATIVE dispersion, direct log-sigma describes ABSOLUTE
% residual SD. Prefer this parameterization whenever the output is described as
% absolute precision, individual error variance, or Gaussian-comparable
% dependability.
%
% CV EQUIVALENCE -- note there is deliberately NO free nu here. nu is derived
% (nu = 2*(mu/sigma)^2), so the sibling's CV = sqrt(2/nu) becomes simply
%   CV(Y | mu, sigma) = sigma/mu,
% i.e. the coefficient of variation is a DERIVED quantity that varies across the
% mean surface rather than a modeled constant. The sibling's "CV = 1/3 at
% nu = 18" anchor therefore has no counterpart parameter in this file; its
% translation onto this scale is the log-sigma prior centre
% (psyrat_default_priors, grep sigma_center), which reproduces nu = 18 at the
% reference point by construction.
%
% The scope of this parameterization is NOT limited to dimension slopes. Since
% log sigma = log mu + 0.5*log(2) - 0.5*log(nu), every log-nu contrast blends
% amplitude with absolute residual SD -- whether indexed by a covariate (dynamic
% designs), by an event or cell (static difference designs), or by a person
% (subject-level designs). See gamma_scale_submodel_decision.md section H.
%
% ESTIMAND. e_cond_var is the MARGINAL expected conditional variance
% E[sigma_p^2] over the population of person residual SDs, so the components sum
% to the true observed-score variance by the law of total variance - the same
% convention PSYRAT_GAMMA_VARCOMPS and the Gaussian path use. Passing s_sd = 0
% (the default) recovers the homogeneous-residual plug-in.
%
% NOTE there is deliberately NO cov_pv input. The log-nu sibling needs one
% because its residual carries exp(-2*cov_pv): the person mean and log-nu
% effects both enter 2*mu^2/nu. Here the residual never touches mu, so the
% person mean <-> scale correlation cannot affect E[sigma_p^2]. The correlation
% is still estimated (it is the LKJ block in the Stan model) and is still
% reported; it simply does not enter this conversion. Adding a cov_pv argument
% "for symmetry" would be wrong.
%
% Inputs (scalars or equal-length column vectors of posterior draws; scalars
% broadcast against vectors):
%   alpha     - log expected-score intercept AT THE EVALUATION POINT z, i.e.
%               Intercept + Xdim(z)*b, not the bare grand intercept
%   s_p       - person SD on the log expected-score scale (gro_sds(1))
%   s_i       - trial  SD on the log expected-score scale (sig_trl)
%   log_sigma - log residual SD AT z, i.e. Intercept_sigma + Xdim(z)*b_sigma.
%               Pass the population value; the person effect enters through
%               s_sd, or is added by the caller for a per-participant read-out.
%   s_sd      - person log-sigma SD (optional; default 0). Person heterogeneity
%               in precision.
%
% Output struct vc, each field a column vector matching the draw length. Field
% names match PSYRAT_GAMMA_VARCOMPS exactly so downstream consumers are
% parameterization-agnostic:
%   mu_bar      - marginal expected observed score, E[mu] = exp(alpha + 0.5*ssum)
%   sigma_p2    - observed-score person main-effect variance component
%   sigma_i2    - observed-score trial  main-effect variance component
%   var_mu      - variance of the lognormal mean surface, Var(mu)
%   e_cond_var  - expected conditional (observation) variance, E[Var(Y|mu,sigma)]
%   total_var   - total observed variance, Var(mu) + E[Var(Y|mu,sigma)]
%   sigma_pi_e2 - residual = person x trial interaction + observation dispersion,
%                 recovered by subtraction: total_var - sigma_p2 - sigma_i2
%
% The mean-surface fields (mu_bar, sigma_p2, sigma_i2, var_mu) are IDENTICAL to
% PSYRAT_GAMMA_VARCOMPS's, because the two parameterizations share the log-mean
% submodel. Only e_cond_var, and hence total_var and sigma_pi_e2, differ.
%
% See also PSYRAT_GAMMA_VARCOMPS, PSYRAT_GAMMA_CROSSCOV, PSYRAT_REL_SING.

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
if nargin < 5 || isempty(s_sd), s_sd = 0; end
alpha     = alpha(:);
s_p       = s_p(:);
s_i       = s_i(:);
log_sigma = log_sigma(:);
s_sd      = s_sd(:);

% Broadcast every input to the common draw length UP FRONT. The log-nu sibling
% gets this for free because each of its outputs touches alpha, so implicit
% expansion sizes them all. Here e_cond_var depends only on the scale submodel,
% so a scalar log_sigma against a vector alpha would otherwise return a scalar
% residual next to vector mean-surface terms and silently break the "every field
% is a column vector matching the draw length" contract.
lengths = [numel(alpha) numel(s_p) numel(s_i) numel(log_sigma) numel(s_sd)];
ndraw   = max(lengths);
if any(lengths ~= 1 & lengths ~= ndraw)
    error('psyrat_gamma_varcomps_ls:draws',... %Error code and associated error
        ['Posterior-draw inputs must be scalars or all the same length. '...
        'How to fix: pass each of alpha, s_p, s_i, log_sigma and s_sd as '...
        'either a scalar or a column vector of the same number of draws.']);
end
grow      = ones(ndraw,1);
alpha     = alpha     .* grow;
s_p       = s_p       .* grow;
s_i       = s_i       .* grow;
log_sigma = log_sigma .* grow;
s_sd      = s_sd      .* grow;

% Squared log-scale SDs and their sum (variance of log mu at a single cell).
sp2  = s_p.^2;
si2  = s_i.^2;
ssum = sp2 + si2;

% Marginal expected observed score: E[mu] for a lognormal mean surface.
mu_bar = exp(alpha + 0.5.*ssum);

% Observed-score person main effect: Var_p( E_t[mu | u_p] ). The inner
% expectation E_t[mu|u_p] = exp(alpha + u_p + si2/2) is lognormal in u_p, whose
% variance is exp(2*alpha + si2) * exp(sp2) * (exp(sp2) - 1). Identical to the
% log-nu sibling: the mean submodel is shared.
sigma_p2 = exp(2.*alpha + si2) .* exp(sp2) .* (exp(sp2) - 1);

% Observed-score trial main effect: the symmetric expression in the trial SD.
sigma_i2 = exp(2.*alpha + sp2) .* exp(si2) .* (exp(si2) - 1);

% Variance of the full lognormal mean surface mu = exp(alpha + u_p + u_t):
% Var(mu) = E[mu^2] - E[mu]^2 = exp(2*alpha)*(exp(2*ssum) - exp(ssum)).
var_mu = exp(2.*alpha) .* (exp(2.*ssum) - exp(ssum));

% Expected conditional (observation-level) variance, MARGINALIZED over the
% person population's residual-SD heterogeneity:
%   E[Var(Y|mu,sigma_p)] = E[sigma_p^2] = E[exp(2*log_sigma + 2*w_p)]
%     = exp(2*log_sigma + 2*s_sd^2)
% by the lognormal second moment (E[exp(2*w)] = exp(2*s_sd^2) for w ~ N(0,
% s_sd^2)). It is exp(2*log_sigma) when s_sd = 0. Note what is ABSENT compared
% with the log-nu sibling: no alpha, no ssum, no cov_pv. That is the decoupling
% this parameterization buys.
e_cond_var = exp(2.*log_sigma + 2.*s_sd.^2);

% Total observed variance by the law of total variance.
total_var = var_mu + e_cond_var;

% Highest-order residual component (P x trial interaction + observation
% dispersion), recovered by subtraction. Clip tiny negative values arising only
% from floating-point noise; a genuinely negative result would indicate a
% specification error upstream. Matches PSYRAT_GAMMA_VARCOMPS's handling.
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
