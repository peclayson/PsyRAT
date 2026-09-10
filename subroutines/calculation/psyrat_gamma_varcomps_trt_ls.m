function vc = psyrat_gamma_varcomps_trt_ls(alpha, s_p, s_o, s_t, s_po, s_pt, ...
    s_ot, log_sigma, s_sd)
%PSYRAT_GAMMA_VARCOMPS_TRT_LS Crossed Persons x Occasions x Trials observed-
% score-scale variance components for the LOCATION-SCALE scaled-chi-square /
% Gamma(log link) family (gammascale = 2).
%
%   vc = psyrat_gamma_varcomps_trt_ls(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot,
%                                     log_sigma)
%   vc = psyrat_gamma_varcomps_trt_ls(..., log_sigma, s_sd)
%
% This is the two-facet sibling of PSYRAT_GAMMA_VARCOMPS_LS and the
% location-scale sibling of PSYRAT_GAMMA_VARCOMPS_TRT. All three describe the
% same scaled-chi-square likelihood; they differ in which quantity carries the
% linear predictor and the person random effect. Here the model is
%   log mu_{pot} = alpha + u_p + u_o + u_t + u_po + u_pt + u_ot,
%   log sigma_p  = log_sigma + w_p,     w_p ~ N(0, s_sd^2),
%   nu_{pot}     = 2*(mu_{pot}/sigma_p)^2,
%   E(Y | mu, sigma) = mu,   Var(Y | mu, sigma) = sigma_p^2,
% with mutually independent zero-mean normal effects (person, occasion, trial,
% and the three two-way interactions) whose SDs are on the LOG expected-score
% scale.
%
% WHY THIS EXISTS. Under the log-nu sibling the conditional variance is
% 2*mu^2/nu, which is MEAN-COUPLED: it moves with the cell's conditional mean,
% so a log-nu contrast blends amplitude with absolute residual SD and cannot be
% read as a statement about precision alone. That matters acutely for this
% design, which is used to build a DIFFERENCE score: the blended contrast is the
% EVENT contrast, and the event contrast is the construct. Since
% log sigma_e = log mu_e + 0.5*log(2) - 0.5*log(nu_e), a fitted log-nu event
% difference b_nu,2 - b_nu,1 equals 2*(b_2 - b_1) - 2*(b_sigma,2 - b_sigma,1) and
% so moves when the two events differ in AMPLITUDE even at identical residual
% precision. Here the residual is decoupled from the mean surface entirely, so
% the per-event log_sigma difference is a pure absolute residual-SD contrast.
% See documentation/gamma_scale_submodel_decision.md section H (owner ruling H).
%
% The log-nu sibling is NOT invalid -- it is a coherent constant-CV
% (proportional-error) model. The difference is what the scale parameter means:
% log-nu describes RELATIVE dispersion, direct log-sigma describes ABSOLUTE
% residual SD.
%
% CV EQUIVALENCE -- there is deliberately NO free nu here. nu is derived
% (nu = 2*(mu/sigma)^2), so the sibling's CV = sqrt(2/nu) becomes simply
%   CV(Y | mu, sigma) = sigma/mu,
% a DERIVED quantity that varies across the mean surface rather than a modeled
% constant. Note in particular what the log-nu two-facet model cannot represent
% and this one can: an occasion whose absolute error SD is unchanged while its
% mean shifts. There, one nu is shared across occasions and trials and every
% facet sits on the log mean, so the residual SD is forced to track amplitude.
%
% ESTIMAND. e_cond_var is the MARGINAL expected conditional variance
% E[sigma_p^2] over the population of person residual SDs, so the components sum
% to the true observed-score variance by the law of total variance -- the same
% convention PSYRAT_GAMMA_VARCOMPS_TRT and the Gaussian path use. Passing
% s_sd = 0 (the default) recovers the homogeneous-residual plug-in.
%
% NOTE there is deliberately NO cov_pv input, exactly as in the one-facet
% location-scale sibling. The log-nu version needs one because its residual
% carries exp(0.5*s_v^2 - 2*cov_pv): the person mean and log-nu effects both
% enter 2*mu^2/nu. Here the residual never touches mu, so the person mean <->
% scale correlation cannot affect E[sigma_p^2]. That correlation is still
% estimated (it is the LKJ person block in the Stan model) and is still reported
% in the provenance lines; it simply does not enter this conversion. Adding a
% cov_pv argument "for symmetry" would be wrong.
%
% RELATION TO THE SUBJECT-LEVEL TWO-FACET PATH. Analysis 25 under gammascale = 2
% (psyrat_gamma_extract_trt_sserr_ls) reuses PSYRAT_GAMMA_VARCOMPS_TRT unchanged,
% with NaN dispersion arguments, because it reads only the six signal components
% -- which are parameterization-invariant -- and gets its residual straight from
% the fitted per-person log-SD. This file exists because the GROUP-level two-facet
% difference (analysis 10) additionally consumes sigma_res2, which is NOT
% parameterization-invariant. The two facts are consistent; do not "simplify" one
% into the other.
%
% Inputs (scalars or equal-length column vectors of posterior draws; scalars
% broadcast against vectors), all mean-submodel inputs on the log expected-score
% scale:
%   alpha     - log expected-score grand intercept (per event, b(:,e))
%   s_p       - person            SD on the log expected-score scale
%   s_o       - occasion          SD
%   s_t       - trial             SD
%   s_po      - occasion x person SD
%   s_pt      - trial x person    SD
%   s_ot      - trial x occasion  SD
%   log_sigma - log residual SD (the per-event b_sigma intercept). Pass the
%               population value; person heterogeneity enters through s_sd, or is
%               added by the caller for a per-participant read-out.
%   s_sd      - person log-sigma SD (optional; default 0)
%
% Output struct vc, each field a column vector matching the draw length. Field
% names match PSYRAT_GAMMA_VARCOMPS_TRT exactly so downstream consumers are
% parameterization-agnostic:
%   mu_bar      - marginal expected observed score, E[mu] = exp(alpha + V/2)
%   sigma_p2    - observed-score person   main-effect variance component
%   sigma_o2    - observed-score occasion main-effect variance component
%   sigma_t2    - observed-score trial    main-effect variance component
%   sigma_po2   - observed-score occasion x person interaction component
%   sigma_pt2   - observed-score trial x person    interaction component
%   sigma_ot2   - observed-score trial x occasion  interaction component
%   var_mu      - variance of the lognormal mean surface, Var(mu)
%   e_cond_var  - expected conditional (observation) variance, E[Var(Y|mu,sigma)]
%   total_var   - total observed variance, Var(mu) + E[Var(Y|mu,sigma)]
%   sigma_res2  - residual = three-way P x O x T interaction + observation
%                 dispersion, recovered by subtraction
%
% The mean-surface fields (mu_bar, the six components, var_mu) are IDENTICAL to
% PSYRAT_GAMMA_VARCOMPS_TRT's, because the two parameterizations share the
% log-mean submodel. Only e_cond_var, and hence total_var and sigma_res2, differ.
%
% REDUCTION TO THE ONE-FACET SIBLING, stated precisely because the obvious
% version of it is wrong. With s_o = s_po = s_ot = s_pt = 0 this function agrees
% with PSYRAT_GAMMA_VARCOMPS_LS on mu_bar, sigma_p2, var_mu, e_cond_var,
% total_var and sigma_t2 <-> sigma_i2. The residuals are NOT equal, because they
% are not the same quantity: sigma_pi_e2 there is "person x trial + dispersion",
% while here person x trial is broken out as its own sigma_pt2. The identity is
%   sigma_res2 + sigma_pt2 = sigma_pi_e2,
% and sigma_pt2 is nonzero even at s_pt = 0, since by inclusion-exclusion the
% observed-scale interaction is B*(exp(s_p^2)-1)*(exp(s_t^2)-1) whenever both
% main effects are present. Note s_pt = 0 is REQUIRED and is not implied by the
% other three: the one-facet model carries no log-scale person x trial effect, so
% a nonzero s_pt inflates V and only e_cond_var still agrees, the identity
% failing by 50% at s_pt = 0.25. The log-nu pair satisfies the same identity, so
% this is a property of the observed-scale decomposition, not of this file.
% Pinned by tests/TestGammaDiffTrtLsConversion.
%
% See also PSYRAT_GAMMA_VARCOMPS_LS, PSYRAT_GAMMA_VARCOMPS_TRT,
% PSYRAT_GAMMA_CROSSCOV_TRT, PSYRAT_REL_TRT.

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
if nargin < 9 || isempty(s_sd), s_sd = 0; end
alpha     = alpha(:);
s_p       = s_p(:);
s_o       = s_o(:);
s_t       = s_t(:);
s_po      = s_po(:);
s_pt      = s_pt(:);
s_ot      = s_ot(:);
log_sigma = log_sigma(:);
s_sd      = s_sd(:);

% Broadcast every input to the common draw length UP FRONT, for the same reason
% the one-facet location-scale sibling does. The log-nu sibling gets this for
% free because each of its outputs touches alpha, so implicit expansion sizes
% them all. Here e_cond_var depends only on the scale submodel, so a scalar
% log_sigma against a vector alpha would otherwise return a scalar residual next
% to vector mean-surface terms and silently break the "every field is a column
% vector matching the draw length" contract.
lengths = [numel(alpha) numel(s_p) numel(s_o) numel(s_t) numel(s_po) ...
    numel(s_pt) numel(s_ot) numel(log_sigma) numel(s_sd)];
ndraw   = max(lengths);
if any(lengths ~= 1 & lengths ~= ndraw)
    error('psyrat_gamma_varcomps_trt_ls:draws',... %Error code and associated error
        ['Posterior-draw inputs must be scalars or all the same length. '...
        'How to fix: pass each of alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, '...
        'log_sigma and s_sd as either a scalar or a column vector of the same '...
        'number of draws.']);
end
grow      = ones(ndraw,1);
alpha     = alpha     .* grow;
s_p       = s_p       .* grow;
s_o       = s_o       .* grow;
s_t       = s_t       .* grow;
s_po      = s_po      .* grow;
s_pt      = s_pt      .* grow;
s_ot      = s_ot      .* grow;
log_sigma = log_sigma .* grow;
s_sd      = s_sd      .* grow;

% ---- Mean submodel. Identical algebra to PSYRAT_GAMMA_VARCOMPS_TRT; these
% components cannot depend on the scale parameterization.

% Squared log-scale SDs (variances of each independent effect on the log mean).
vp  = s_p.^2;
vo  = s_o.^2;
vt  = s_t.^2;
vpo = s_po.^2;
vpt = s_pt.^2;
vot = s_ot.^2;

% V = variance of log mu at a single cell = sum of all six log-scale variances.
V = vp + vo + vt + vpo + vpt + vot;

% Common prefactor B = exp(2*alpha + V) = E[mu]^2.
B = exp(2.*alpha + V);

% Marginal expected observed score: E[mu] for the lognormal mean surface.
mu_bar = exp(alpha + 0.5.*V);

% Main-effect components. Cov of two cells sharing exactly one factor equals that
% factor's main-effect component; for a single shared factor the only nested
% effect is the main effect itself, so Cov = B*(exp(v_main)-1).
sigma_p2 = B .* (exp(vp) - 1);
sigma_o2 = B .* (exp(vo) - 1);
sigma_t2 = B .* (exp(vt) - 1);

% Two-way interaction components by inclusion-exclusion. Cov of two cells sharing
% factor set {x,y} nests the two main effects plus the xy interaction, so the
% interaction component is that covariance minus the two main-effect components.
% The max(.,0) only guards float noise, as in the log-nu sibling.
sigma_po2 = max(B .* (exp(vp + vo + vpo) - exp(vp) - exp(vo) + 1), 0);
sigma_pt2 = max(B .* (exp(vp + vt + vpt) - exp(vp) - exp(vt) + 1), 0);
sigma_ot2 = max(B .* (exp(vo + vt + vot) - exp(vo) - exp(vt) + 1), 0);

% Variance of the full lognormal mean surface:
% Var(mu) = E[mu^2] - E[mu]^2 = exp(2*alpha+2V) - exp(2*alpha+V) = B*(exp(V)-1).
var_mu = B .* (exp(V) - 1);

% ---- Scale submodel. This is the ONLY place the two parameterizations differ.

% Expected conditional (observation-level) variance, MARGINALIZED over the person
% population's residual-SD heterogeneity:
%   E[Var(Y|mu,sigma_p)] = E[sigma_p^2] = E[exp(2*log_sigma + 2*w_p)]
%     = exp(2*log_sigma + 2*s_sd^2)
% by the lognormal second moment (E[exp(2*w)] = exp(2*s_sd^2) for w ~ N(0,
% s_sd^2)). It is exp(2*log_sigma) when s_sd = 0. Note what is ABSENT compared
% with the log-nu sibling: no alpha, no V, no cov_pv. In particular the five
% within-person facet variances do not appear. Under log-nu they must, because
% the residual 2*mu^2/nu is an expectation over those facets; here sigma_p is a
% single per-participant quantity carrying no facet terms, so there is no
% expectation to take. That decoupling is what this parameterization buys, and
% restoring an exp(2V) factor "for symmetry" would re-introduce exactly the mean
% coupling it exists to remove.
e_cond_var = exp(2.*log_sigma + 2.*s_sd.^2);

% Total observed variance by the law of total variance.
total_var = var_mu + e_cond_var;

% Highest-order residual component (three-way P x O x T interaction + scaled-chi-
% square observation dispersion), recovered by subtraction. Clip tiny negative
% values arising only from floating-point noise; a genuinely negative result
% would indicate a specification error upstream.
sigma_res2 = max(total_var - sigma_p2 - sigma_o2 - sigma_t2 ...
    - sigma_po2 - sigma_pt2 - sigma_ot2, 0);

vc = struct( ...
    'mu_bar',     mu_bar, ...
    'sigma_p2',   sigma_p2, ...
    'sigma_o2',   sigma_o2, ...
    'sigma_t2',   sigma_t2, ...
    'sigma_po2',  sigma_po2, ...
    'sigma_pt2',  sigma_pt2, ...
    'sigma_ot2',  sigma_ot2, ...
    'var_mu',     var_mu, ...
    'e_cond_var', e_cond_var, ...
    'total_var',  total_var, ...
    'sigma_res2', sigma_res2);

end
