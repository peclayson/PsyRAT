function psyrat_prefs = psyrat_defaults
%Default settings for the PsyRAT Toolbox
%
%

%This script will be read by the PsyRAT Toolbox to define the default settings
%for analyzing and viewing data in the PsyRAT Toolbox. Should the user want to
%change any default settings, the values below can be changed. Upon
%restarting the toolbox, these settings will be read. Comments are provided 
%for understanding each input.
%
%Note: updating the toolbox overwrites this file, so re-apply any local
%edits after an update.
%
%Input
% No inputs required in the command line
%
%Output
% psyrat_prefs - data structure containing the default processing and viewing
%  settings for the PsyRAT Toolbox (the per-field defaults set below)
%

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

psyrat_prefs = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Processing Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of chains for processing data in CmdStan (default = 4)
psyrat_prefs.proc.nchains = 4;

%Number of warmup iterations for processing data in CmdStan
psyrat_prefs.proc.nwarmup = 5000;

%Number of sampling iterations for processing data in CmdStan
psyrat_prefs.proc.nsampling = 5000;

%Total iterations retained for backwards compatibility (warmup + sampling)
psyrat_prefs.proc.niter = psyrat_prefs.proc.nwarmup + psyrat_prefs.proc.nsampling;

%NUTS adaptation target acceptance rate (adapt_delta) passed to CmdStan.
%Empty ([]) means "use PsyRAT's per-analysis default" (each model sets its own
%tuned value, e.g. 0.90/0.95/0.98). Set a scalar in (0,1) - typically 0.90-0.999
%- to override that default; raising it reduces divergent transitions at the
%cost of slower sampling. Opt-in: leaving this empty preserves current behavior.
psyrat_prefs.proc.adapt_delta = [];

%NUTS maximum tree depth (max_treedepth) passed to CmdStan. Empty ([]) means
%"use CmdStan's default of 10". Set a positive integer (typically 10-15) to
%override; raising it avoids premature trajectory termination on hard
%posteriors. Opt-in: leaving this empty preserves current behavior.
psyrat_prefs.proc.max_treedepth = [];

%Verbose input while processing data in CmdStan (default = 2);
% 1 = No
% 2 = Yes
psyrat_prefs.proc.verbose = 2;

%View the traceplots for parameters prior to saving the stan output 
%(default= 1)
% 1 = No
% 2 = Yes
psyrat_prefs.proc.traceplots = 1;

%%%%Design-selection and estimation preferences (difference scores, %%%%%%%%%%
%%%%subject-level, family, engine, seed, priors). Legal combinations are %%%%%
%%%%decided by the analysis selector in psyrat_computevarcomp; comments %%%%%%
%%%%here deliberately do not enumerate supported designs. %%%%%%%%%%%%%%%%%%%%
%Whether subject-level reliability metrics should be estimated
% 1 = No
% 2 = Yes
psyrat_prefs.proc.ssubjrel = 1;
%legacy alias retained for backwards compatibility
psyrat_prefs.proc.sserrvar = psyrat_prefs.proc.ssubjrel;

%Observation-likelihood family for the single-trial scores
% 'gaussian' = default (ERP amplitudes)
% 'gamma'    = scaled chi-square / log-linked mean, for strictly positive,
%              right-skewed data (single-trial time-frequency power); available
%              for the design families accepted by the estimator's
%              family/gammascale guards (the guard in psyrat_computevarcomp is
%              authoritative; documentation/model_availability_matrix.md is the
%              user-facing matrix) (CmdStan engine only)
psyrat_prefs.proc.family = 'gaussian';

%Gamma scale parameterization (the 'gammascale' option). Applies under the Gamma
%family only; inert under Gaussian.
% []  = unset / Automatic (the default). The warp OMITS the argument, so
%       psyrat_computevarcomp applies its context-aware default: location-scale
%       on every design that supports it, log-nu on the rest (rulings 33-36 in
%       documentation/gamma_scale_submodel_decision.md).
% 1   = log-nu dispersion (relative dispersion / constant CV)
% 2   = log residual SD (location-scale, absolute residual SD)
%
%Left EMPTY rather than set to 1 on purpose: a concrete value here would pin the
%estimand for every run and bypass the engine's per-design resolution, which is
%the single authority on which designs accept which parameterization.
psyrat_prefs.proc.gammascale = [];

%Whether the internal consistency of difference scores should be estimated
% 1 = No
% 2 = Yes, 2-event difference score
% 3 = Yes, 4-event difference-of-differences
psyrat_prefs.proc.diffest = 1;

%DoD contrast mapping for processing mode (diffest = 3)
%Stored event labels in ERP order: {ERP_1 ERP_2 ERP_3 ERP_4}
psyrat_prefs.proc.dodmap = {};
psyrat_prefs.proc.dodmapcol = [];
psyrat_prefs.proc.dodset = 0;

%Whether correlated residuals should be estimated for difference score
%models (only relevant when diffest = 2; a diffest = 3
%difference-of-differences run rejects a concurrent request outright)
% 1 = No
% 2 = Yes
psyrat_prefs.proc.diffrescor = 1;

%Whether within-person residual covariance should be estimated for
%co-occurring difference score events
% 1 = No (set wp_cov = 0)
% 2 = Yes
psyrat_prefs.proc.diffwpcov = 1;

%Whether within-chain parallelization should be used for Stan models
%(default = 1)
% 1 = No
% 2 = Yes
%When enabled, PsyRAT will automatically pick threads_per_chain so that at
%least 2 CPU cores remain free after chains*threads_per_chain.
psyrat_prefs.proc.withinchain = 1;

%RNG seed passed to CmdStan for reproducible estimation runs
psyrat_prefs.proc.seed = 12345;

%Estimation engine selection. PsyRAT estimates variance components with CmdStan
%(the preferred Bayesian engine) or, when CmdStan is unavailable, with native
%MATLAB engines. Availability is detected by
%psyrat_estimation_engines_available and dispatched by psyrat_run_native.
% 0 = auto    : try CmdStan, then the native HMC engine, then the native fitlme
%               engine (the cascade), considering only engines installed on
%               this machine. CmdStan is therefore the default whenever it is
%               present (the run is unchanged); a clear, non-suppressible
%               warning is issued whenever the cascade falls back off CmdStan.
% 1 = cmdstan : force CmdStan (errors if it is missing or too old).
% 2 = fitlme  : force the native frequentist engine (simple designs only; its
%               intervals are approximate, not Bayesian credible intervals).
% 3 = hmc     : force the native Hamiltonian Monte Carlo engine.
%A forced value (1/2/3) is an explicit user opt-in and overrides the cascade.
psyrat_prefs.proc.engine = 0;

%Number of parametric-bootstrap resamples the native fitlme engine uses to turn
%its REML point estimates into aligned draw vectors, so the downstream
%calculation layer can still form intervals by propagating the joint
%distribution through the (nonlinear) reliability formula. Ignored by the
%CmdStan and HMC engines, which produce posterior draws directly.
psyrat_prefs.proc.bootstrap_reps = 1000;

%Gradient method for the native HMC engine (engine 3).
% 'auto'  = automatic differentiation (dlarray) when the Deep Learning Toolbox
%           is available, otherwise numerical gradients.
% 'ad'    = force automatic differentiation (requires Deep Learning Toolbox).
% 'numer' = force numerical (finite-difference) gradients.
psyrat_prefs.proc.hmc_gradient = 'auto';

%Population-level priors for the CmdStan estimation models (means/intercepts
%and variance-component standard deviations). Defaults are weakly informative
%on the ERP microvolt scale; edit the values returned by psyrat_default_priors
%to adapt to a different measurement scale. See psyrat_default_priors and
%documentation/priors_and_sensitivity.md.
psyrat_prefs.proc.priors = psyrat_default_priors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Viewing Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Value to be used for dependability cutoff/threshold
psyrat_prefs.view.depvalue = .8;

%For the inputs below a value of 1 (default) indicates that the
%figure/table should be viewed; a value of 0 indicates the figure/table
%should not be viewed
%The default value for all of these figures/tables is 1

%Plot that shows the relationship between the number of trials retained 
%for averaging and dependability
psyrat_prefs.view.plotdep = 1;

%Plot that compares the intraclass correlation coefficients for each
%group and/or condition
psyrat_prefs.view.ploticc = 1;

%Subject-level reliability/dependability plot toggle
psyrat_prefs.view.plotssrel = 1;

%Subject-level ICC plot toggle
psyrat_prefs.view.plotssicc = 1;

%Subject-level table toggle
psyrat_prefs.view.tablessrel = 1;

%Table displaying information about cutoffs based on dependability
%threshold
psyrat_prefs.view.inctrltable = 1;

%Table displaying information about overall dependability with data
%including all trials after applying the trial cutoffs
psyrat_prefs.view.overalltable = 1;

%Table displaying information about between- and within-person standard
%deviations as well as intraclass correlation coefficients
psyrat_prefs.view.showstddevt = 1;

%Figure displaying between-person standard deviations for each group and/or
%condition
psyrat_prefs.view.showstddevf = 1;


%Which dependability estimate to plot on the figure that shows the
%relationship between the number of trials retained for averaging and
%dependability (default = 2)
% 1 = lower limit of the credible interval
% 2 = point estimate of the credible interval
% 3 = upper limit of the credible interval
psyrat_prefs.view.plotdepline = 2;

%Number of trials to plot on the figure that shows the relationship between
%the number of trials retained for averaging and dependability 
%(default = 50)
psyrat_prefs.view.ntrials = 50;

%Which cutoff to use for estimating the dependability of averages for a
%trial cutoff (default = 2);
% 1 = lower limit of the credible interval
% 2 = point estimate of the credible interval
% 3 = upper limit of the credible interval
psyrat_prefs.view.meascutoff = 2;

%Measure of central tendency to use for estimating the overall score
%dependability of waveforms after applying the trial cutoffs (default = 1)
% 1 = mean
% 2 = median
psyrat_prefs.view.depcentmeas = 1;

%%%%Additional preferences relevant to only analyses on data with %%%%%%%%%
%%%% multiple occasions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Which G-theory coefficient should be calculated
% 1 = dependability
% 2 = generalizability
psyrat_prefs.view.gcoeff = 1;

%Which reliability coefficient should be calculated
% 1 = equivalence
% 2 = stability
% 3 = equivalence and stability
psyrat_prefs.view.reltype = 1;

%Test-retest occasion estimand (how the occasion facet enters the D study)
% 1 = single-occasion score (n'_o = 1), the default: reliability of a score
%     from n'_i trials at one occasion, generalizing over the occasion facet
% 2 = multi-occasion composite score (n'_o = k): reliability of a score
%     averaged across k occasions. The SEM uses the same n'_o as the
%     coefficient (author estimand convention; scientific audit sec.14).
psyrat_prefs.view.noccmode = 1;

%Number of occasions (k) for the multi-occasion composite score. Empty means
%use the number of observed occasions. Only used when noccmode == 2.
psyrat_prefs.view.nocc = [];

%Value to be used for reliability cutoff/threshold
psyrat_prefs.view.relvalue = .8;

%Difference-score coefficient for single-session difference-score analyses
% 1 = dependability (absolute decision)
% 2 = generalizability (relative decision)
psyrat_prefs.view.diffgcoeff = 1;

%DoD view defaults (task-free 4-event difference-of-differences)
%Contrast is always (ERP_1 - ERP_2) - (ERP_3 - ERP_4), with user-selected
%event mapping in the DoD view window.
psyrat_prefs.view.dodgroup = 1;
psyrat_prefs.view.dodmap = [1 2 3 4];
psyrat_prefs.view.dodtrials = [];
psyrat_prefs.view.dodgridmax = [];
psyrat_prefs.view.dodheatmetric = 1; %1 = dependability, 2 = generalizability

%Plot that shows the relationship between the number of trials retained 
%for averaging and reliability
psyrat_prefs.view.plotrel = 1;

%Which reliability estimate to plot on the figure that shows the
%relationship between the number of trials retained for averaging and
%reliability (default = 2)
% 1 = lower limit of the credible interval
% 2 = point estimate of the credible interval
% 3 = upper limit of the credible interval
psyrat_prefs.view.plotrelline = 2;

%Measure of central tendency to use for estimating the overall score
%reliability of waveforms after applying the trial cutoffs (default = 1)
% 1 = mean
% 2 = median
psyrat_prefs.view.relcentmeas = 1;

%Criterion cutoff to use for criterion-score-specific reliability outputs
psyrat_prefs.view.criterioncutoff = 0;

%Criterion output toggles
% 1 = include, 0 = do not include
psyrat_prefs.view.plotcriterion = 1;
psyrat_prefs.view.tablecriterion = 1;

end
