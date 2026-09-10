function priors = psyrat_default_priors
%Default population-level priors for the CmdStan estimation models.
%
% priors = psyrat_default_priors
%
%Returns a struct of the population-level (location and scale) prior scales
%used by psyrat_computevarcomp. These are the user-configurable priors: the
%population intercepts/means and the variance-component standard deviations.
%The structural non-centered reparameterization priors (the raw effects'
%normal(0,1)) are not exposed here. The Gaussian-family LKJ correlation
%priors are fixed, except the four-event difference-of-differences model,
%which exposes its LKJ shape as priors.dod.lkj; the gamma family exposes its
%two LKJ shapes as priors.gamma.lkj and priors.gamma.modular_lkj (see the
%gamma block below).
%
%IMPORTANT: these defaults are weakly informative on the ERP microvolt scale.
%They are not scale-free. For other measurement scales or modalities, set the
%values to match the data scale (see documentation/priors_and_sensitivity.md).
%
%Fields (each value is the scale parameter of the named prior):
% single  - single-occasion model (also used by the group/event variants)
%   mu     : mu ~ normal(0, mu)
%   sig_u  : sig_u ~ cauchy(0, sig_u)        between-person SD
%   sig_trl: sig_trl ~ cauchy(0, sig_trl)    trial main-effect SD (sigma_i)
%   sig_e  : sig_e ~ cauchy(0, sig_e)        person x trial + residual SD
% sserr   - subject-level (single-subject error-variance) location-scale model
%   sig_trl: sig_trl ~ cauchy(0, sig_trl)    trial main-effect SD (sigma_i).
%            The population intercept, between-person, and residual log-SD
%            priors of the sserr model stay fixed brms-style constants; only
%            the added crossed trial main effect is exposed here.
% dynrel  - dynamic/conditional reliability (Rast & Clayson) one-facet
%           location-scale model with fixed dimension slopes on the mean and
%           the log-residual. As with sserr, only these scales are exposed:
%   sig_trl: sig_trl ~ cauchy(0, sig_trl)    trial main-effect SD (sigma_i)
%   b      : b ~ normal(0, b)                mean dimension slopes
%   b_sigma: b_sigma ~ normal(0, b_sigma)    log-residual dimension slopes
%            (scales are on the standardized predictor scale, not microvolts)
% trt     - test-retest crossed design (Persons x Trials x Occasions)
%   sig_id, sig_occ, sig_trl, sig_err, sig_trlxid, sig_occxid, sig_trlxocc
%   (each: param ~ cauchy(0, value))
% diff    - difference-score / location-scale models
%   b           : b ~ normal(0, b)
%   b_sigma     : b_sigma ~ student_t(3, 0, b_sigma)
%   sd_id       : sd_id ~ student_t(3, 0, sd_id)        (half)
%   sd_trl      : sd_trl ~ student_t(3, 0, sd_trl)      (half)
%   sd_id_sigma : sd_id_sigma ~ student_t(3, 0, sd_id_sigma) (half)
%   rescor_mu   : rescor_mu ~ normal(0, rescor_mu)
%   sd_rescor   : sd_rescor ~ student_t(3, 0, sd_rescor) (half)
% diff_trt - two-facet (test-retest) difference-score model (Persons x Trials
%   x Occasions, two conditions). The six crossed cross-condition components
%   each get a half-student_t SD prior; the per-condition residual log-SD and
%   the population cell means share the diff-style scales. The person x trial x
%   occasion interaction is lumped into the residual (sigma_poi,e), so it has no
%   separate component prior. Scales mirror priors.trt on the microvolt scale.
%   b      : b ~ normal(0, b)                       cell means (per condition)
%   b_sigma: b_sigma ~ student_t(3, 0, b_sigma)     residual log-SD (per cond.)
%   sd_id  : person component sigma(p)        ~ student_t(3, 0, sd_id)   (half)
%   sd_trl : trial component sigma(i)         ~ student_t(3, 0, sd_trl)  (half)
%   sd_occ : occasion component sigma(o)      ~ student_t(3, 0, sd_occ)  (half)
%   sd_tid : person x trial sigma(pi)         ~ student_t(3, 0, sd_tid)  (half)
%   sd_oid : person x occasion sigma(po)      ~ student_t(3, 0, sd_oid)  (half)
%   sd_to  : trial x occasion sigma(io)       ~ student_t(3, 0, sd_to)   (half)
% dod     - four-event difference-of-differences model (analysis 9, Gaussian;
%   psyrat_build_stan_dodiff and its native HMC twin local_hmc_dodiff). Until
%   2026-09-05 that builder hard-coded brms-style constants and a 'priors'
%   override never reached it; these fields expose the same values, so a
%   default run emits byte-identical Stan. Scales are on the microvolt scale
%   of the four cells; the student-t degrees of freedom (10) stay fixed.
%   b_cell       : b_cell ~ normal(0, b_cell)                     cell means
%   b_sigma_cell : b_sigma_cell ~ student_t(10, 0, b_sigma_cell)  cell log-residual SDs
%   sd_id        : sd_id ~ student_t(10, 0, sd_id)        (half)  the 8 person SDs
%                  (location and scale block)
%   sd_trial     : sd_trial ~ student_t(10, 0, sd_trial)  (half)  the 8 trial SDs
%   lkj          : LKJ(lkj) on the 8-D person and trial correlation matrices
% splits  - nonparallel split (observed-design) models. The extra interaction
%   SDs the weighted split-mean models add on top of their single/test-retest
%   donors; the donors' shared components reuse priors.single / priors.trt.
%   sig_splitxid : person x split sigma(ps)            ~ cauchy(0, sig_splitxid)
%                  [case 23, p x (i:s)]
%   sig_posxid   : person x occasion x split sigma(pos)~ cauchy(0, sig_posxid)
%                  [case 24, p x o x (i:s)]
% gamma   - scaled chi-square observation-family priors (means, log-nu and
%           location-scale SD scales, copula and modular-cut scales, and the
%           two LKJ shapes). Documented inline at the assignments below;
%           see also documentation/priors_and_sensitivity.md.

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

priors = struct();

% single-occasion (shared by group/event variants)
priors.single.mu = 100;
priors.single.sig_u = 40;
priors.single.sig_trl = 10;   % trial main-effect SD (sigma_i); matches priors.trt.sig_trl
priors.single.sig_e = 40;

% subject-level error-variance (sserr) model. Only the crossed trial
% main-effect SD (sigma_i) is configurable; the population intercept,
% between-person, and residual log-SD priors remain the fixed brms-style
% constants built into psyrat_build_stan_sserr. Deliberately NOT brms's tight
% student_t(10,0,2): a weakly informative half-cauchy on the microvolt scale,
% matching priors.single.sig_trl.
priors.sserr.sig_trl = 10;

% dynamic/conditional reliability (Rast & Clayson) location-scale model. The
% one-facet sserr model (person intercept + per-subject log-residual via a 2-D
% Cholesky) extended with fixed dimension slopes on BOTH the mean (b) and the
% log-residual (b_sigma). KDIM is data-driven: 1 slope for one dimension; 3 for
% two dimensions (z1, z2, z1*z2). As with sserr, the population intercept,
% between-person, and residual log-SD priors stay the fixed brms-style
% constants in psyrat_build_stan_dynrel; only the trial main-effect SD and the
% two dimension-slope blocks are configurable here. The dimension predictors
% are standardized (mean 0, SD 1), so these scales are on the standardized
% predictor scale, not the microvolt scale.
priors.dynrel.sig_trl = 10;  % trial main-effect SD (sigma_i); matches sserr
priors.dynrel.b = 10;        % mean dimension slopes:  b ~ normal(0, b)
priors.dynrel.b_sigma = 1;   % scale dimension slopes: b_sigma ~ normal(0, b_sigma)

% test-retest crossed design
priors.trt.sig_id = 10;
priors.trt.sig_occ = 5;
priors.trt.sig_trl = 10;
priors.trt.sig_err = 20;
priors.trt.sig_trlxid = 5;
priors.trt.sig_occxid = 5;
priors.trt.sig_trlxocc = 1;

% difference-score / location-scale models
priors.diff.b = 5;
priors.diff.b_sigma = 10;
priors.diff.sd_id = 10;
priors.diff.sd_trl = 10;
priors.diff.sd_id_sigma = 2.5;
priors.diff.rescor_mu = 1;
priors.diff.sd_rescor = 1;

% two-facet (test-retest) difference-score model. Six crossed cross-condition
% components plus a per-condition residual; the pto interaction is lumped into
% the residual (Rocha Table 3 / build_stan_trt convention). Component scales
% mirror priors.trt; b/b_sigma mirror priors.diff.
priors.diff_trt.b = 5;
priors.diff_trt.b_sigma = 10;
priors.diff_trt.sd_id = 10;   % person sigma(p)            (matches priors.trt.sig_id)
priors.diff_trt.sd_trl = 10;  % trial sigma(i)             (matches priors.trt.sig_trl)
priors.diff_trt.sd_occ = 5;   % occasion sigma(o)          (matches priors.trt.sig_occ)
priors.diff_trt.sd_tid = 5;   % person x trial sigma(pi)   (matches priors.trt.sig_trlxid)
priors.diff_trt.sd_oid = 5;   % person x occasion sigma(po)(matches priors.trt.sig_occxid)
priors.diff_trt.sd_to = 1;    % trial x occasion sigma(io) (matches priors.trt.sig_trlxocc)

% four-event difference-of-differences model (analysis 9, Gaussian). These
% reproduce the brms-style constants psyrat_build_stan_dodiff hard-coded until
% 2026-09-05 (owner ruling: expose before the beta), so the emitted Stan at the
% defaults is byte-identical to the former model (pinned by the dod_event
% golden). The native HMC twin (psyrat_native_hmc, local_hmc_dodiff) reads the
% same block through native_extra. The student-t degrees of freedom are 10 (a
% fixed brms-style constant, not the 3 the diff family uses) and stay fixed;
% only the scales are exposed. Scales are on the microvolt scale of the
% four cells. Override per run via the 'priors' input.
priors.dod.b_cell = 10;       % cell means:           b_cell ~ normal(0, b_cell)
priors.dod.b_sigma_cell = 2;  % cell log-residual SDs: b_sigma_cell ~ student_t(10, 0, b_sigma_cell)
priors.dod.sd_id = 2;         % person SDs (8-D block): sd_id ~ student_t(10, 0, sd_id) (half)
priors.dod.sd_trial = 2;      % trial SDs (8-D block):  sd_trial ~ student_t(10, 0, sd_trial) (half)
priors.dod.lkj = 2;           % LKJ(eta) on the person and trial 8x8 correlations

% nonparallel split (observed-design) models. These add the explicit
% person x split [case 23, p x (i:s)] and person x occasion x split
% [case 24, p x o x (i:s)] interaction SDs that the weighted split-mean models
% estimate ON TOP of their single/test-retest donors. The donors' shared
% components reuse priors.single / priors.trt unchanged; only these two extra
% interaction scales are configurable here. Both match the microvolt-scale
% interaction defaults (priors.trt.sig_trlxid / priors.trt.sig_occxid).
priors.splits.sig_splitxid = 5;  % person x split sigma(ps)        [case 23]
priors.splits.sig_posxid = 5;    % person x occasion x split sigma(pos) [case 24]

% Gamma / scaled-chi-square family (one-facet, log-linked mean + person-varying
% log-nu dispersion). IMPORTANT: unlike the Gaussian priors above, these live on
% the LOG expected-score scale (mean/SD parameters) and the LOG degrees-of-
% freedom scale (nu), NOT the microvolt scale. They are strongly scale-dependent
% and must track the measurement scale; the defaults reproduce the reference
% single-trial time-frequency-power model (mean center exp(1.7047) = 5.5, nu
% center exp(2.8904) = 18). See documentation/priors_and_sensitivity.md.
priors.gamma.mu_center = log(5.5);  % log-mean grand intercept: Intercept ~ normal(mu_center, mu_sd)
priors.gamma.mu_sd     = 0.5;       % log-mean grand intercept SD
priors.gamma.sig_u     = 0.35;      % person mean log-SD:  sig_u   ~ student_t(3, 0, sig_u)  (half)
priors.gamma.sig_trl   = 0.35;      % trial  mean log-SD:  sig_trl ~ student_t(3, 0, sig_trl)(half)
% Test-retest (case 5) crossed-design log-mean facet SDs. The reference uses a
% single student_t(3, 0, 0.35) scale for all crossed random-effect SDs; these are
% split per facet so they can be tuned independently (mirroring the Gaussian
% test-retest per-facet priors). Each is half student_t(3, 0, scale).
priors.gamma.sig_occ     = 0.35;    % occasion       mean log-SD (case 5)
priors.gamma.sig_trlxid  = 0.35;    % trial x person mean log-SD (case 5)
priors.gamma.sig_occxid  = 0.35;    % occasion x person mean log-SD (case 5)
priors.gamma.sig_trlxocc = 0.35;    % trial x occasion mean log-SD (case 5)
priors.gamma.nu_center = log(18);   % log-nu grand intercept: Intercept_nu ~ normal(nu_center, nu_sd)
priors.gamma.nu_sd     = 0.75;      % log-nu grand intercept SD
% Person log-nu SD prior scale: sig_nu ~ student_t(3, 0, sig_nu) (half). Set to
% 1.0 (not the brms response-scale default 2.5): sig_nu lives on the LOG-nu scale,
% where 2.5 implies a prior-median s_v ~ 1.9, i.e. dispersion varying ~45x across
% people - an implausible default that only became consequential once the marginal
% observation term (estimand #2, see documentation/priors_and_sensitivity.md and
% SCIENTIFIC_FORMULA_AUDIT.md sections 18/11) began using s_v. sig_nu = 1.0 is
% better calibrated (median s_v ~ 0.77), reduces low-trial residual over-inflation,
% and does not over-shrink genuine heterogeneity at typical sample sizes (validated
% by the two identifiability simulations recorded in section 11). Override per run
% via the 'priors' input if a wider/narrower dispersion prior is warranted.
priors.gamma.sig_nu    = 1.0;       % person log-nu SD:   sig_nu  ~ student_t(3, 0, sig_nu)  (half)
% Copula-scoped person log-nu SD prior. The concurrent copula difference designs
% (cases 7/10/8, diffrescor=2) add a per-person log-nu funnel that is markedly harder
% to identify than the non-copula gamma models; the owner-supplied Stan reference
% (stan_chisquare_ssrel/...subject_specific...) tightens the log-nu SD prior to
% student_t(3, 0, 0.5) there for identifiability (S11). This field is consumed ONLY by
% the copula builders (psyrat_build_stan_diff_copula_chisq / _trt_copula_chisq); every
% non-copula gamma design keeps the validated sig_nu = 1.0 above. Override per run via
% the 'priors' input if a wider/narrower copula dispersion prior is warranted.
priors.gamma.sig_nu_copula = 0.5;   % copula-only person log-nu SD: ~ student_t(3, 0, sig_nu_copula) (half)
priors.gamma.lkj       = 1;         % LKJ(eta) on the 2-D person (mean, log-nu) correlation
% Dimension slopes for the gamma DYNAMIC/conditional reliability models (case 11
% one-facet, case 14 trial x occasion). The gamma dynrel adds fixed slopes on the
% standardized dimension predictor(s) to BOTH log-linked submodels: the log mean
% (b) and the log degrees-of-freedom / dispersion (b_nu), mirroring how the
% Gaussian dynrel adds priors.dynrel.{b,b_sigma}. KDIM is data-driven (1 slope for
% one dimension; 3 for two: z1, z2, z1*z2). Because the predictors are
% standardized (mean 0, SD 1), these scales are per-z-SD movement on the LOG
% expected-score (b) and LOG-nu (b_nu) scales, NOT microvolts. Both are centered
% at 0 (null = no dimensional effect) and weakly informative: N(0, 0.5) allows a
% one-SD change in the predictor to shift log-mu (or log-nu) by up to ~+/-1 in the
% 95% prior interval, i.e. a mean (or nu) change of up to ~2.7x per z-SD. Override
% per run via the 'priors' input if a wider/narrower dimensional effect is
% warranted. (Owner-confirmable default; see SCIENTIFIC_FORMULA_AUDIT.md sec 18.)
priors.gamma.b    = 0.5;   % mean dimension slopes:  b    ~ normal(0, b)    (log expected-score scale)
priors.gamma.b_nu = 0.5;   % log-nu dimension slopes: b_nu ~ normal(0, b_nu)(log-nu scale)
% LOCATION-SCALE gamma priors (gammascale = 2). This parameterization replaces the
% log-nu scale submodel with a direct log residual-SD submodel,
%   log sigma = Intercept_sigma + Xdim*b_sigma + ind_sd[id],   nu = 2*(mu/sigma)^2,
% so that b_sigma is a PURE residual-SD slope rather than the mean/dispersion
% blend b_nu carries (documentation/gamma_scale_submodel_decision.md).
%
% SCOPE (widened 2026-08-04, owner ruling H). These are no longer dynrel-only.
% Every gammascale = 2 design draws on this block; the STATIC designs (subject-
% level 6/25, difference 7/8/9/10) use sigma_center, sigma_sd and sig_sd but have
% no dimension slope, so b_sigma is simply unused there rather than being a
% separate default. The reason the static designs need this parameterization at
% all is that log sigma = log mu + 0.5*log(2) - 0.5*log(nu) makes EVERY log-nu
% contrast a blend of amplitude and absolute residual SD -- whether the contrast
% is indexed by a covariate, by an event/cell, or by a person.
%
% DERIVATION. Priors are NOT invariant under this reparameterization, so these
% are derived rather than copied, using the exact relation
%   log sigma = log mu + 0.5*log(2) - 0.5*log(nu).
% Centering: sigma_center = mu_center + 0.5*log(2/18) = log(5.5) - log(3) =
% log(1.8333), which reproduces exactly nu = 18 at the reference point, i.e. the
% same implied dispersion as nu_center above. (Independently, the owner's Stan
% reference bundle centers its log-sigma offset at log(mean) + 0.5*log(2/18) --
% the identical convention, arrived at separately.)
% Scales: the 0.5 coefficient on log nu halves the log-nu scales, so sig_nu = 1.0
% implies ~0.5 for the person log-sigma SD, and nu_sd = 0.75 implies ~0.375 for
% the intercept. sigma_sd is set to 0.5 rather than 0.375 to stay weakly
% informative and to match mu_sd; sig_sd takes the derived 0.5 (the reference
% bundle uses a tighter 0.35 for both, so these defaults are the more
% conservative choice).
% b_sigma is NOT a translation: the slope relation b_nu = 2*b - 2*b_sigma cannot
% be inverted into a prior (b and b_nu are specified independently above, and the
% induced value would be ~0.56). It is instead chosen for interpretability on its
% own scale: N(0, 0.25) lets a one-SD change in the predictor shift log-sigma by
% up to ~+/-0.5 in the 95% prior interval, i.e. a residual-SD change of up to
% ~1.65x per z-SD - generous for a precision effect without being unbounded.
% All four warrant a prior-sensitivity check; see
% documentation/priors_and_sensitivity.md. Override per run via 'priors'.
priors.gamma.sigma_center = log(5.5) + 0.5*log(2/18); % log residual-SD intercept: Intercept_sigma ~ normal(sigma_center, sigma_sd)
priors.gamma.sigma_sd     = 0.5;   % log residual-SD intercept SD
priors.gamma.sig_sd       = 0.5;   % person log-sigma SD: sig_sd ~ student_t(3, 0, sig_sd) (half)
priors.gamma.b_sigma      = 0.25;  % log residual-SD dimension slopes: b_sigma ~ normal(0, b_sigma)
% Gaussian-copula residual correlation for the CONCURRENT gamma difference score
% (cases 7/10, diffrescor=2): rho_e ~ normal(0, rescor_sd), bounded (-0.95, 0.95).
% rho_e couples the two concurrently recorded events' scaled-chi-square residuals
% (psyrat_build_stan_diff_copula_chisq); it is the copula analog of the Gaussian
% difference model's LKJ residual correlation. Centered at 0 (null = independent
% residuals) and weakly informative on the (-1, 1) correlation scale, matching the
% owner reference (model_fits_concurrent_theta_alpha_copula_chisq). NOTE the copula
% correlation is NOT the observed-score residual correlation (the marginals are
% non-Gaussian); the observed-scale residual cross-covariance is derived by Monte
% Carlo in psyrat_gamma_extract_diff_copula. Override per run via 'priors'.
% (Owner-confirmable default; see SCIENTIFIC_FORMULA_AUDIT.md sec 22.)
priors.gamma.rescor_sd = 0.5;  % copula residual corr: rho_e ~ normal(0, rescor_sd)
% MODULAR CUT-COPULA stage-1 margin priors (analyses 19/20, family = 'gamma',
% gammascale = 2). These are a SEPARATE block rather than a reuse of the
% location-scale fields above, and the reason is structural, not stylistic: that
% model centres its intercepts on fixed literal values (sigma_center, mu_center),
% whereas the modular margins model passes DATA-DERIVED offsets into Stan and
% fits the intercepts as DEVIATIONS from them (psyrat_build_stan_diffdynrel_sserr_
% margins_chisq). A prior centred at log(5.5) is meaningless on a deviation, so
% the centred fields cannot be borrowed even though the likelihood is the same
% family.
%
% The values reproduce the owner's reference bundle
% (person_specific_dynamic_concurrent_modular/R/model_generation.R) exactly,
% because that is the configuration its completed HPC fits were validated under
% and the one the frozen external-oracle fixtures were produced from. Changing
% them silently makes the toolbox estimate a different model from the reference
% it is checked against; run a sensitivity analysis and record it first.
%
% All are on LOG scales: modular_b and modular_b_dim on the log expected-score
% scale, the rest on the log residual-SD scale. The two LKJ blocks use eta = 2
% (mildly favouring lower correlations) rather than the eta = 1 uniform default
% used elsewhere in this family - again, the reference's choice.
priors.gamma.modular_b           = 0.75; % log-mean intercept deviation:  b ~ normal(0, .)
priors.gamma.modular_b_dim       = 0.50; % log-mean dimension slopes:     b_dim* ~ normal(0, .)
priors.gamma.modular_b_sigma     = 0.50; % log residual-SD intercept dev: b_sigma ~ normal(0, .)
priors.gamma.modular_b_sigma_dim = 0.35; % log residual-SD dim slopes:    b_sigma_dim* ~ normal(0, .)
priors.gamma.modular_sd_id       = 0.50; % person LOCATION SDs (sd_id[1:2]) ~ student_t(3, 0, .)
priors.gamma.modular_sd_id_sigma = 0.35; % person SCALE SDs   (sd_id[3:4]) ~ student_t(3, 0, .)
priors.gamma.modular_sd_facet    = 0.50; % every crossed facet SD ~ student_t(3, 0, .)
priors.gamma.modular_lkj         = 2;    % LKJ(eta) on the person and facet correlations

end
