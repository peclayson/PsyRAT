function vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu, s_v, cov_pv)
%PSYRAT_GAMMA_VARCOMPS_TRT Crossed Persons x Occasions x Trials observed-score-
% scale variance components for the scaled-chi-square / Gamma(log link) family
% (test-retest, analysis 5).
%
%   vc = psyrat_gamma_varcomps_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu)
%   vc = psyrat_gamma_varcomps_trt(..., nu, s_v, cov_pv)
%
% The Gamma family models the single-trial score with a LOG-LINKED mean over the
% fully crossed test-retest design,
%   log mu_{pot} = alpha + u_p + u_o + u_t + u_po + u_pt + u_ot,
% with mutually independent zero-mean normal effects (person, occasion, trial,
% and the three two-way interactions) whose SDs (s_p, s_o, s_t, s_po, s_pt, s_ot)
% are on the LOG expected-score scale, and a scaled-chi-square (equivalently
% Gamma) observation with
%   E(Y | mu, nu) = mu,   Var(Y | mu, nu) = 2*mu^2/nu.
%
% WHAT nu IS, PSYCHOMETRICALLY. As in the one-facet sibling, Var = 2*mu^2/nu makes
% this a CONSTANT COEFFICIENT-OF-VARIATION (proportional-error) model:
%   CV(Y | mu, nu) = sqrt(2/nu),   log CV = 0.5*log(2) - 0.5*log(nu),
% which is 1/3 at the default nu = 18. So a log-nu contrast describes RELATIVE
% dispersion, not absolute precision, and log sigma = log mu + 0.5*log(2) -
% 0.5*log(nu) means every log-nu contrast blends amplitude with absolute residual
% SD. Note this design shares ONE nu across occasions and trials -- every facet
% sits on the log mean -- so the model cannot represent an occasion whose absolute
% error SD is unchanged while its mean shifts. See the full note in
% PSYRAT_GAMMA_VARCOMPS and documentation/gamma_scale_submodel_decision.md.
%
% Because mu is on the log scale, the Stan/brms group-level SDs are NOT the
% G-theory variance components. This function maps them to OBSERVED-score-scale
% variance components using lognormal moments. It is the crossed-design analog of
% PSYRAT_GAMMA_VARCOMPS (one-facet); see REDUCTION TO THE ONE-FACET SIBLING below
% for the precise sense in which it reduces to it.
%
% DERIVATION (closed form; see the header block of the accompanying test).
% For the log-linked mean mu = exp(alpha + sum of independent normal effects),
% the covariance of two cells that share exactly the factor set A (and differ
% independently on the complementary facets) is
%   Cov_A = B * (exp( sum_{S subset of A} s_S^2 ) - 1),   B = exp(2*alpha + V),
% where V is the sum of ALL six log-scale variances and the inner sum runs over
% effects whose index set is a subset of A. In a random-effects ANOVA the
% covariance of two cells sharing factor set A equals the sum of the variance
% components indexed within A, so the observed-score components follow by
% inclusion-exclusion:
%   sigma_p^2  = Cov_{p}   = B*(exp(s_p^2)-1)              (and symmetric for o,t)
%   sigma_po^2 = Cov_{po} - sigma_p^2 - sigma_o^2
%              = B*( exp(s_p^2+s_o^2+s_po^2) - exp(s_p^2) - exp(s_o^2) + 1 )
% and symmetrically for person x trial and occasion x trial. The residual
% (three-way P x O x T interaction confounded with the scaled-chi-square
% observation dispersion) is recovered by subtraction from the total observed
% variance, exactly as the one-facet residual confounds P x T with dispersion.
%
% Inputs (scalars or equal-length column vectors of posterior draws; scalars
% broadcast against vectors), all SD/alpha inputs on the log expected-score scale:
%   alpha - log expected-score grand intercept (b_Intercept)
%   s_p   - person             SD (sd_ID__Intercept)
%   s_o   - occasion           SD (sd_time__Intercept)
%   s_t   - trial              SD (sd_TrialNumber__Intercept)
%   s_po  - occasion x person  SD (person-by-occasion interaction)
%   s_pt  - trial x person     SD (person-by-trial interaction)
%   s_ot  - trial x occasion   SD (occasion-by-trial interaction)
%   nu    - population scaled-chi-square degrees of freedom (exp of the log-nu
%           intercept)
%   s_v    - person log-nu SD (optional; default 0)
%   cov_pv - Cov(u_p, v_p), person mean<->log-nu covariance on the log scale
%            (optional; default 0). The person log-nu effect couples only to the
%            person mean (not occasion/trial/interactions), so the marginal
%            observation term carries the SAME factor exp(0.5*s_v^2 - 2*cov_pv)
%            as the one-facet case (estimand #2; default 0 recovers the
%            population-nu plug-in, estimand #1).
%
% Output struct vc, each field a column vector matching the draw length:
%   mu_bar      - marginal expected observed score, E[mu] = exp(alpha + V/2)
%   sigma_p2    - observed-score person   main-effect variance component
%   sigma_o2    - observed-score occasion main-effect variance component
%   sigma_t2    - observed-score trial    main-effect variance component
%   sigma_po2   - observed-score occasion x person interaction variance component
%   sigma_pt2   - observed-score trial x person    interaction variance component
%   sigma_ot2   - observed-score trial x occasion  interaction variance component
%   var_mu      - variance of the lognormal mean surface, Var(mu)
%   e_cond_var  - expected conditional (observation) variance, E[Var(Y|mu,nu)]
%   total_var   - total observed variance, Var(mu) + E[Var(Y|mu,nu)]
%   sigma_res2  - residual = three-way P x O x T interaction + observation
%                 dispersion, recovered by subtraction:
%                 total_var - (the six components above)
%
% Downstream, the observed-score components are treated as a standard
% random-effects ANOVA decomposition and projected across trials/occasions with
% the usual linear D-study law (PSYRAT_REL_TRT), which squares SD inputs
% internally; the caller therefore passes sqrt(component).
%
% REDUCTION TO THE ONE-FACET SIBLING, stated precisely because the obvious
% version of it is wrong. With s_o = s_po = s_ot = s_pt = 0 this function agrees
% with PSYRAT_GAMMA_VARCOMPS on mu_bar, sigma_p2, var_mu, e_cond_var, total_var
% and sigma_t2 <-> sigma_i2. The residuals are NOT equal, because they are not the
% same quantity: sigma_pi_e2 there is "person x trial + dispersion", while here
% person x trial is broken out as its own sigma_pt2. The identity is
%   sigma_res2 + sigma_pt2 = sigma_pi_e2,
% and sigma_pt2 is nonzero even at s_pt = 0, since by inclusion-exclusion the
% observed-scale interaction is B*(exp(s_p^2)-1)*(exp(s_t^2)-1) whenever both
% main effects are present (0.187, or 3.4% of the one-facet residual, at the
% pinned design point). Note s_pt = 0 is REQUIRED and is not implied by the other
% three: the one-facet model carries no log-scale person x trial effect, so a
% nonzero s_pt inflates V and NO field agrees, the identity included. The
% location-scale pair satisfies the same identity, so this is a property of the
% observed-scale decomposition, not of either parameterization. Pinned by
% tests/TestGammaDiffTrtLsConversion.
%
% See also PSYRAT_GAMMA_VARCOMPS, PSYRAT_GAMMA_VARCOMPS_TRT_LS,
% PSYRAT_GAMMA_CROSSCOV, PSYRAT_REL_TRT.

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
if nargin < 9  || isempty(s_v),    s_v    = 0; end
if nargin < 10 || isempty(cov_pv), cov_pv = 0; end
alpha  = alpha(:);
s_p    = s_p(:);
s_o    = s_o(:);
s_t    = s_t(:);
s_po   = s_po(:);
s_pt   = s_pt(:);
s_ot   = s_ot(:);
nu     = nu(:);
s_v    = s_v(:);
cov_pv = cov_pv(:);

% Squared log-scale SDs (variances of each independent effect on the log mean).
vp  = s_p.^2;
vo  = s_o.^2;
vt  = s_t.^2;
vpo = s_po.^2;
vpt = s_pt.^2;
vot = s_ot.^2;

% V = variance of log mu at a single cell = sum of all six log-scale variances.
V = vp + vo + vt + vpo + vpt + vot;

% Common prefactor B = exp(2*alpha + V) = E[mu]^2 (see derivation above).
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
% Each is nonnegative analytically ((exp(v_x)-1)(exp(v_y)-1) when the interaction
% SD is zero, larger otherwise); the max(.,0) only guards float noise and mirrors
% the reference implementation (twofacet_components_from_params).
sigma_po2 = max(B .* (exp(vp + vo + vpo) - exp(vp) - exp(vo) + 1), 0);
sigma_pt2 = max(B .* (exp(vp + vt + vpt) - exp(vp) - exp(vt) + 1), 0);
sigma_ot2 = max(B .* (exp(vo + vt + vot) - exp(vo) - exp(vt) + 1), 0);

% Variance of the full lognormal mean surface:
% Var(mu) = E[mu^2] - E[mu]^2 = exp(2*alpha+2V) - exp(2*alpha+V) = B*(exp(V)-1).
var_mu = B .* (exp(V) - 1);

% Expected conditional (observation-level) variance, MARGINALIZED over the
% person population's dispersion heterogeneity (estimand #2):
%   E[Var(Y|mu,nu_p)] = E[2*mu^2/nu_p]
%     = (2/nu)*exp(2*alpha + 2V) * exp(0.5*s_v^2 - 2*cov_pv).
% The multi-facet SDs enter only the base term exp(2*alpha+2V); the correction
% factor is identical to the one-facet case because the person log-nu effect v_p
% couples only to the person mean u_p (Var(2*sum_of_mean_effects - v_p) contains
% just the one cross-covariance Cov(u_p,v_p) = cov_pv). Default s_v = 0 recovers
% the population-nu plug-in, matching PSYRAT_GAMMA_VARCOMPS.
e_cond_var = (2./nu) .* exp(2.*alpha + 2.*V) .* exp(0.5.*s_v.^2 - 2.*cov_pv);

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
