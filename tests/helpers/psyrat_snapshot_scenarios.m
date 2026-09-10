function scenarioMap = psyrat_snapshot_scenarios()
%Scenarios for the generated-Stan snapshot tests.
%
% Returns the CmdStan integration scenarios (PsyRATCmdStanIntegrationFactory)
% plus a difference-of-differences (diffest=3) case. The integration factory
% does not cover difference-of-differences because the shared test data has
% only two event types; the snapshot test builds a four-event table inline for
% the DoD scenario (flagged with isDoD) and passes the matching dodmap.

scenarioMap = PsyRATCmdStanIntegrationFactory.scenarios();

dod = struct();
dod.name = 'dod_event';
dod.includeGroup = false;
dod.includeEvent = true;
dod.includeTime = false;
dod.sserrvar = 1;
dod.diffest = 3;
dod.diffwpcov = 1;
dod.diffrescor = 1;
dod.expectedAnalysis = 'ic_dodiff';
dod.expectedErrorId = '';
dod.isDoD = true;
dod.dodmap = {'E1', 'E2', 'E3', 'E4'};

scenarioMap.dod_event = dod;

% Gamma / scaled-chi-square one-facet family (analysis 1, family='gamma'). The
% family swaps the observation likelihood (log-linked mean + person-varying
% log-nu dispersion + scaled chi-square) but keeps the one-facet, no-group,
% no-event routing. The generated Stan text is independent of the data values,
% so the runner supplies a strictly positive table (the family's support) via
% the isGamma branch.
gamma1 = struct();
gamma1.name = 'ic_gamma_base';
gamma1.includeGroup = false;
gamma1.includeEvent = false;
gamma1.includeTime = false;
gamma1.sserrvar = 1;
gamma1.diffest = 1;
gamma1.diffwpcov = 1;
gamma1.diffrescor = 1;
gamma1.family = 'gamma';
gamma1.isGamma = true;
gamma1.expectedAnalysis = 'ic';
gamma1.expectedErrorId = '';

scenarioMap.ic_gamma_base = gamma1;

% Gamma family, grouped (analysis 2) and event (analysis 3): the per-stratum
% gamma facet builder (psyrat_build_stan_facet_chisq). Same isGamma positive
% table, with a group or event column driving the routing.
gammaG = gamma1;
gammaG.name = 'ic_gamma_group';
gammaG.includeGroup = true;
scenarioMap.ic_gamma_group = gammaG;

gammaE = gamma1;
gammaE.name = 'ic_gamma_event';
gammaE.includeEvent = true;
scenarioMap.ic_gamma_event = gammaE;

% Gamma family, single-occasion SUBJECT-LEVEL design (analysis 6, sserr): the
% gamma subject-level builder (psyrat_build_stan_sserr_chisq). sserrvar=2 routes
% to analysis 6. The gamma sserr model IS the one-facet gamma model (the
% person-varying log-nu already gives each participant their own residual), so the
% emitted Stan matches ic_gamma_base; the distinct scenario pins the sserr builder
% + golden and the analysis-6 (ic_sserrvar) routing. Same isGamma positive table.
%
% WHY THIS SCENARIO PINS gammascale EXPLICITLY (rulings 33-36, 2026-08-07). The
% engine default is now location-scale for every design that supports it, and
% analysis 6 is one of them. A scenario that omitted gammascale would therefore
% generate the LOCATION-SCALE model and no longer match its log-nu golden. The
% explicit 1 keeps this scenario pinning the log-nu builder it was written for,
% and keeps its golden byte-identical.
%
% Only the ELEVEN movable designs are pinned this way. The gamma scenarios for
% analyses 1-5 and the concurrent copula variants are deliberately left WITHOUT a
% gammascale field, because the default cannot move for them - so they keep
% working as tripwires that would fail if a future edit widened the flip too far.
gammaSserr = gamma1;
gammaSserr.name = 'ic_gamma_sserr';
gammaSserr.sserrvar = 2;
gammaSserr.expectedAnalysis = 'ic_sserrvar';
gammaSserr.gammascale = 1;
scenarioMap.ic_gamma_sserr = gammaSserr;

% Gamma family, one-facet SUBJECT-LEVEL under the LOCATION-SCALE parameterization
% (analysis 6, gammascale = 2). The scale submodel becomes a direct log residual
% SD (log sigma = Intercept_sigma + ind_sd[id], nu derived as 2*(mu/sigma)^2), so
% unlike the log-nu scenario above this is NOT byte-identical to ic_gamma_base:
% analyses 1-4 have no location-scale variant, being single-mean designs whose
% coefficients the parameterization leaves unchanged. It is therefore a scenario
% of its own and is deliberately NOT added to the delegation pair list in
% TestStanModelSnapshots, which pins the gammascale = 1 identities only.
gammaSserrLs = gammaSserr;
gammaSserrLs.name = 'ic_gamma_ls_sserr';
gammaSserrLs.gammascale = 2;
scenarioMap.ic_gamma_ls_sserr = gammaSserrLs;

% Gamma family, two-event difference score (analysis 7): the joint bivariate
% gamma builder (psyrat_build_stan_diff_chisq). Non-concurrent (diffwpcov=1).
gammaDiff = gamma1;
gammaDiff.name = 'ic_gamma_diff_event';
gammaDiff.includeEvent = true;
gammaDiff.diffest = 2;
gammaDiff.expectedAnalysis = 'ic_diff';
gammaDiff.gammascale = 1; %movable design; pins log-nu - see ic_gamma_sserr above
scenarioMap.ic_gamma_diff_event = gammaDiff;

% Gamma family, single-occasion SUBJECT-LEVEL two-event difference score (analysis
% 8, ic_diff_sserrvar): sserrvar=2 + diffest=2 routes to case 8. The gamma
% subject-level difference builder (psyrat_build_stan_diff_sserr_chisq) emits the
% SAME joint bivariate model as the group-level gamma difference (case 7); the
% subject-level residual is a downstream extraction difference, so the golden is
% byte-identical to ic_gamma_diff_event. The distinct scenario pins the case-8
% builder + golden and the ic_diff_sserrvar routing. Non-concurrent (diffrescor=1).
gammaDiffSserr = gammaDiff;
gammaDiffSserr.name = 'ic_gamma_diff_sserr';
gammaDiffSserr.sserrvar = 2;
gammaDiffSserr.expectedAnalysis = 'ic_diff_sserrvar';
gammaDiffSserr.gammascale = 1; %movable design; pins log-nu
scenarioMap.ic_gamma_diff_sserr = gammaDiffSserr;

% Gamma family, two-event difference score under the LOCATION-SCALE scale
% parameterization (analysis 7, gammascale = 2, NON-concurrent). Same design and
% same data as ic_gamma_diff_event; only the scale submodel differs, so the
% emitted mean submodel is byte-identical and the diff between the two goldens is
% the estimand-invariance argument made visible.
%
% CONTRAST WITH ic_gamma_ls_sserr ABOVE, which is deliberately excluded from the
% delegation pair list: there the group-level twin (analysis 1) has no
% location-scale variant, so the subject-level golden has nothing to match. Here
% BOTH ends of the delegation have one - psyrat_build_stan_diff_sserr_chisq
% forwards gammascale to psyrat_build_stan_diff_chisq - so the case-8 <- case-7
% byte-identity survives under location-scale and IS pinned in
% TestStanModelSnapshots.
gammaDiffLs = gammaDiff;
gammaDiffLs.name = 'ic_gamma_ls_diff_event';
gammaDiffLs.gammascale = 2;
scenarioMap.ic_gamma_ls_diff_event = gammaDiffLs;

% Gamma family, SUBJECT-LEVEL two-event difference score under the LOCATION-SCALE
% parameterization (analysis 8, gammascale = 2, NON-concurrent). Byte-identical to
% ic_gamma_ls_diff_event by delegation; the subject-level distinction is the
% downstream extraction (psyrat_gamma_extract_diff_sserr_ls reads each person's
% own log residual-SD effect off r_id[,3:4] instead of deriving it from their
% log-mean and log-nu).
gammaDiffSserrLs = gammaDiffSserr;
gammaDiffSserrLs.name = 'ic_gamma_ls_diff_sserr';
gammaDiffSserrLs.gammascale = 2;
scenarioMap.ic_gamma_ls_diff_sserr = gammaDiffSserrLs;

% Gamma family, two-event CONCURRENT difference score (analysis 7, diffrescor=2):
% the Gaussian-copula scaled-chi-square difference builder
% (psyrat_build_stan_diff_copula_chisq). diffwpcov=2 -> diffrescor=2 pairs the two
% events by (person, trial-position) and couples their scaled-chi-square residuals
% with a Gaussian copula correlation rho_e. Same isGamma positive balanced
% two-event table as the non-concurrent gamma difference (ic_gamma_diff_event).
gammaDiffCopula = gammaDiff;
gammaDiffCopula.name = 'ic_gamma_diff_copula';
gammaDiffCopula.diffwpcov = 2;
gammaDiffCopula.diffrescor = 2;
% Pin dispersion=1: the group-level default is now GLOBAL (2), so this per-person-nu
% scenario must select dispersion=1 explicitly to keep exercising the per-person
% builder (psyrat_build_stan_diff_copula_chisq) and its golden.
gammaDiffCopula.dispersion = 1;
scenarioMap.ic_gamma_diff_copula = gammaDiffCopula;

% Gamma family, GROUP-LEVEL GLOBAL fixed-dispersion concurrent one-facet difference
% (analysis 7, diffrescor=2, dispersion=2): the global-nu Gaussian-copula builder
% (psyrat_build_stan_diff_copula_fixnu_chisq). Same paired concurrent routing as
% ic_gamma_diff_copula, but dispersion=2 emits a 2-D person mean block and a single
% event-global log-nu (no per-person log-nu 4-D block), matching the owner reference.
gammaDiffCopulaFixnu = gammaDiffCopula;
gammaDiffCopulaFixnu.name = 'ic_gamma_diff_copula_fixnu';
gammaDiffCopulaFixnu.dispersion = 2;
scenarioMap.ic_gamma_diff_copula_fixnu = gammaDiffCopulaFixnu;

% Gamma family, single-occasion SUBJECT-LEVEL CONCURRENT two-event difference score
% (analysis 8, ic_diff_sserrvar, diffrescor=2): sserrvar=2 + diffest=2 + diffwpcov=2
% routes to case 8 with the concurrent path, which fits the case-7 Gaussian-copula
% model (psyrat_build_stan_diff_copula_chisq) on paired meas1/meas2. The subject-level
% residual (per-person copula residual cross-covariance) is a downstream extraction,
% so the golden is BYTE-IDENTICAL to ic_gamma_diff_copula. The distinct scenario pins
% the case-8 concurrent routing + builder selection.
gammaDiffSserrCopula = gammaDiff;
gammaDiffSserrCopula.name = 'ic_gamma_diff_sserr_copula';
gammaDiffSserrCopula.sserrvar = 2;
gammaDiffSserrCopula.diffwpcov = 2;
gammaDiffSserrCopula.diffrescor = 2;
gammaDiffSserrCopula.expectedAnalysis = 'ic_diff_sserrvar';
scenarioMap.ic_gamma_diff_sserr_copula = gammaDiffSserrCopula;

% Gamma family, crossed test-retest (analysis 5): the gamma test-retest builder
% (psyrat_build_stan_trt_chisq). Same isGamma positive table with a time/occasion
% column driving the two-facet routing.
gammaTrt = gamma1;
gammaTrt.name = 'ic_gamma_trt';
gammaTrt.includeTime = true;
gammaTrt.expectedAnalysis = 'trt';
scenarioMap.ic_gamma_trt = gammaTrt;

% Gamma family, crossed test-retest SUBJECT-LEVEL design (analysis 25,
% trt_sserrvar): sserrvar=2 + a time column routes to case 25. The gamma
% subject-level test-retest builder (psyrat_build_stan_trt_sserr_chisq) emits the
% SAME crossed model as the group-level gamma test-retest (case 5) - the
% person-varying log-nu already gives each participant their own dispersion, so no
% location-scale residual submodel is needed - and the per-subject residual is a
% downstream extraction difference (psyrat_gamma_extract_trt_sserr). The golden is
% therefore BYTE-IDENTICAL to ic_gamma_trt; the distinct scenario pins the case-25
% builder + golden and the trt_sserrvar routing. Same relationship as
% ic_gamma_sserr to ic_gamma_base and ic_gamma_diff_sserr to ic_gamma_diff_event.
gammaTrtSserr = gammaTrt;
gammaTrtSserr.name = 'ic_gamma_trt_sserr';
gammaTrtSserr.sserrvar = 2;
gammaTrtSserr.expectedAnalysis = 'trt_sserrvar';
gammaTrtSserr.gammascale = 1; %movable design; pins log-nu
scenarioMap.ic_gamma_trt_sserr = gammaTrtSserr;

% Gamma family, crossed test-retest SUBJECT-LEVEL under the LOCATION-SCALE
% parameterization (analysis 25, gammascale = 2). Two-facet twin of
% ic_gamma_ls_sserr: the mean submodel (occasion, trial, three interactions) is
% untouched and only the scale submodel changes, so the six observed-scale group
% components are identical to the log-nu form and only the per-participant
% residual estimand differs. Not byte-identical to ic_gamma_trt (analysis 5 has no
% location-scale variant) and so, like its one-facet sibling, deliberately absent
% from the delegation pair list.
gammaTrtSserrLs = gammaTrtSserr;
gammaTrtSserrLs.name = 'ic_gamma_ls_trt_sserr';
gammaTrtSserrLs.gammascale = 2;
scenarioMap.ic_gamma_ls_trt_sserr = gammaTrtSserrLs;

% Gamma family, two-facet (test-retest) two-event difference score (analysis 10):
% the gamma joint two-facet difference builder (psyrat_build_stan_diff_trt_chisq).
% Same isGamma positive table with event + time columns and diffest=2,
% non-concurrent (diffwpcov=1).
gammaDiffTrt = gamma1;
gammaDiffTrt.name = 'ic_gamma_diff_trt';
gammaDiffTrt.includeEvent = true;
gammaDiffTrt.includeTime = true;
gammaDiffTrt.diffest = 2;
gammaDiffTrt.expectedAnalysis = 'trt_diff';
gammaDiffTrt.gammascale = 1; %movable design; pins log-nu
scenarioMap.ic_gamma_diff_trt = gammaDiffTrt;

% Gamma family, two-facet (test-retest) two-event difference score under the
% LOCATION-SCALE parameterization (analysis 10, gammascale = 2, NON-concurrent).
% Same design and same data as ic_gamma_diff_trt; only the scale submodel differs,
% so the emitted mean submodel - the per-event intercepts, the person mean effects
% and ALL FIVE crossed factor blocks - is byte-identical, and the diff between the
% two goldens is the estimand-invariance argument made visible.
%
% NO DELEGATION PAIR, unlike ic_gamma_ls_diff_sserr <- ic_gamma_ls_diff_event.
% Analysis 10 has no static subject-level twin: psyrat_computevarcomp routes
% diffest==2 && sserrvar==2 to analysis 8 regardless of whether a time column is
% present, so there is no two-facet subject-level difference builder for a
% delegation golden to match.
gammaDiffTrtLs = gammaDiffTrt;
gammaDiffTrtLs.name = 'ic_gamma_ls_diff_trt';
gammaDiffTrtLs.gammascale = 2;
scenarioMap.ic_gamma_ls_diff_trt = gammaDiffTrtLs;

% Gamma family, two-facet CONCURRENT difference score (analysis 10, diffrescor=2):
% the Gaussian-copula scaled-chi-square two-facet difference builder
% (psyrat_build_stan_diff_trt_copula_chisq). diffwpcov=2 -> diffrescor=2 pairs the
% two events by (person, occasion, trial-position) and couples their scaled-chi-
% square residuals with a Gaussian copula correlation rho_e. Same isGamma positive
% event+time table as the non-concurrent gamma two-facet difference.
gammaDiffTrtCopula = gammaDiffTrt;
gammaDiffTrtCopula.name = 'ic_gamma_diff_trt_copula';
gammaDiffTrtCopula.diffwpcov = 2;
gammaDiffTrtCopula.diffrescor = 2;
% Pin dispersion=1: the group-level default is now GLOBAL (2), so this per-person-nu
% scenario must select dispersion=1 explicitly to keep exercising the per-person
% builder (psyrat_build_stan_diff_trt_copula_chisq) and its golden.
gammaDiffTrtCopula.dispersion = 1;
scenarioMap.ic_gamma_diff_trt_copula = gammaDiffTrtCopula;

% Gamma family, GROUP-LEVEL GLOBAL fixed-dispersion concurrent two-facet difference
% (analysis 10, diffrescor=2, dispersion=2): the global-nu two-facet Gaussian-copula
% builder (psyrat_build_stan_diff_trt_copula_fixnu_chisq). Same paired concurrent
% event+time routing as ic_gamma_diff_trt_copula, but dispersion=2 emits a 2-D
% person mean block and a single event-global log-nu, matching the owner reference.
gammaDiffTrtCopulaFixnu = gammaDiffTrtCopula;
gammaDiffTrtCopulaFixnu.name = 'ic_gamma_diff_trt_copula_fixnu';
gammaDiffTrtCopulaFixnu.dispersion = 2;
scenarioMap.ic_gamma_diff_trt_copula_fixnu = gammaDiffTrtCopulaFixnu;

% Gamma family, four-event difference-of-differences (analysis 9): the gamma
% joint four-cell DoD builder (psyrat_build_stan_dodiff_chisq). Reuses the DoD
% four-event table (already strictly positive, satisfying the family's support)
% via the isDoD branch and the matching dodmap; family='gamma' selects the
% scaled-chi-square builder over the Gaussian psyrat_build_stan_dodiff.
gammaDod = dod;
gammaDod.name = 'ic_gamma_dodiff';
gammaDod.family = 'gamma';
gammaDod.gammascale = 1; %movable design; pins log-nu
scenarioMap.ic_gamma_dodiff = gammaDod;

% Gamma family, four-event difference-of-differences under the LOCATION-SCALE
% parameterization (analysis 9, gammascale = 2, NON-concurrent). Same design and
% routing as ic_gamma_dodiff; only the scale submodel and its two priors differ,
% which is the point of snapshotting it separately. Diffing this golden against
% ic_gamma_dodiff.stan is the structural form of the estimand-invariance argument:
% the per-cell mean intercepts, the mean half of the 8-D person block, the whole
% 4-D trial block, the mu accumulation, both LKJ priors and the entire generated
% quantities block must be byte-identical, which is what licenses reusing
% psyrat_gamma_crosscov unchanged across all six cell pairs.
%
% Analysis 9 has NO subject-level twin (diffest==2 && sserrvar==2 routes to
% analysis 8) and no concurrent twin, so this scenario adds two goldens rather
% than four and contributes no delegation pair to the builder-twin check.
gammaDodLs = gammaDod;
gammaDodLs.name = 'ic_gamma_ls_dodiff';
gammaDodLs.gammascale = 2;
scenarioMap.ic_gamma_ls_dodiff = gammaDodLs;

% Person-specific dynamic NONCONCURRENT difference-of-differences, one-facet
% (analysis 28, dynrel=2 + diffest=3 + sserrvar=2, family='gamma',
% gammascale=2 - its only supported family/parameterization). Four gamma
% location-scale margins fitted jointly (8-D person block, 4-D crossed-trial
% block, per-cell dimension slopes on both submodels), long format, offsets as
% data, NO dependence module anywhere (the six residual covariances are fixed
% to zero by the estimand; the boundary test asserts the absent tokens). The
% emitted text is KDIM-independent (KDIM is data), so there is no two-dimension
% golden - the same convention as the case-19 margins scenario. The table
% (isDoDDynrel) carries deliberately UNEQUAL per-cell trial counts, which the
% long format must accept (the paired layout errors there) - that property is
% the reason this design does not reuse psyrat_build_diffdynrel_pairs.
dodDynrel = struct();
dodDynrel.name = 'dynrel_gamma_ls_dodiff_sserr_onefacet';
dodDynrel.includeGroup = false;
dodDynrel.includeEvent = true;
dodDynrel.includeTime = false;
dodDynrel.sserrvar = 2;
dodDynrel.diffest = 3;
dodDynrel.diffwpcov = 1;
dodDynrel.diffrescor = 1;
dodDynrel.family = 'gamma';
dodDynrel.gammascale = 2;
dodDynrel.expectedAnalysis = 'ic_dodiff_dynrel_sserrvar';
dodDynrel.expectedErrorId = '';
dodDynrel.isDoD = true;      %adds the dodmap argument
dodDynrel.isDynrel = true;   %adds the dynrel=2 argument
dodDynrel.isDoDDynrel = true; %selects the four-event + dimension table
dodDynrel.dodmap = {'E1', 'E2', 'E3', 'E4'};
scenarioMap.dynrel_gamma_ls_dodiff_sserr_onefacet = dodDynrel;

% The TWO-FACET (trials x occasions) twin, analysis 29: the one-facet margins
% plus four additional 4-D correlated location blocks (occ/tid/oid/to). The
% table adds a two-occasion time column (isDoDDynrelTrt); everything else -
% including the KDIM-independence of the emitted text - carries over.
dodDynrelTrt = dodDynrel;
dodDynrelTrt.name = 'dynrel_gamma_ls_dodiff_sserr_trt_onefacet';
dodDynrelTrt.includeTime = true;
dodDynrelTrt.expectedAnalysis = 'ic_dodiff_dynrel_sserrvar_trt';
dodDynrelTrt.isDoDDynrelTrt = true;
scenarioMap.dynrel_gamma_ls_dodiff_sserr_trt_onefacet = dodDynrelTrt;

% The GAUSSIAN arms of the same two analyses (owner ruling 2026-08-20;
% SCIENTIFIC_FORMULA_AUDIT.md section 26): identical routing inputs minus the
% family/gammascale pair (Gaussian is the default family), reusing the
% four-event + dimension tables. Golden names follow the Gaussian dynrel
% convention (no family infix). Reaching the builder end-to-end is also the
% acceptance pin on the removed family guard - until 2026-08-20 this exact
% configuration errored varargin:familyunsupported.
dodDynrelGauss = dodDynrel;
dodDynrelGauss.name = 'dynrel_dodiff_sserr_onefacet';
dodDynrelGauss = rmfield(dodDynrelGauss, {'family','gammascale'});
scenarioMap.dynrel_dodiff_sserr_onefacet = dodDynrelGauss;

dodDynrelGaussTrt = dodDynrelTrt;
dodDynrelGaussTrt.name = 'dynrel_dodiff_sserr_trt_onefacet';
dodDynrelGaussTrt = rmfield(dodDynrelGaussTrt, {'family','gammascale'});
scenarioMap.dynrel_dodiff_sserr_trt_onefacet = dodDynrelGaussTrt;

% Difference-score reliability WITH residual covariance (analysis 7, rescor
% branch: sserrvar=1, diffwpcov=2 -> diffrescor=2). The factory's only
% diffrescor=2 scenario sets sserrvar=2 and routes to analysis 8, so the
% analysis-7 bivariate (multi_normal_cholesky) variants would otherwise be
% unsnapshotted. This branch was unreachable until the B10 fix. Event-only is
% sufficient: the rescor Stan text does not depend on group. The factory's
% balanced two-event table supplies the within-participant pairs that diffwpcov=2
% preflight requires, so it is reused via makeScenarioTable (no isDoD flag).
rescor = struct();
rescor.name = 'ic_diff_rescor_event';
rescor.includeGroup = false;
rescor.includeEvent = true;
rescor.includeTime = false;
rescor.sserrvar = 1;
rescor.diffest = 2;
rescor.diffwpcov = 2;
rescor.diffrescor = 2;
rescor.expectedErrorId = '';

scenarioMap.ic_diff_rescor_event = rescor;

% Two-facet (test-retest) two-event difference-score reliability (analysis 10:
% diffest=2, sserrvar=1, WITH a time/occasion facet). These are the only
% event+time+diffest=2 scenarios, so without them the two-facet difference
% builders (psyrat_build_stan_diff_trt / _rescor) would be unsnapshotted. The
% non-concurrent form (diffwpcov=1) stacks the two conditions long with 0/1
% indicators; the concurrent form (diffwpcov=2 -> diffrescor=2) pairs them by
% (person, occasion, trial-position) with an estimated residual correlation. The
% factory's balanced two-event table, with the occasion column added, supplies
% the shared (id, occasion, trial-position) keys the concurrent pairing needs.
diffTrt = struct();
diffTrt.name = 'trt_diff_event_time';
diffTrt.includeGroup = false;
diffTrt.includeEvent = true;
diffTrt.includeTime = true;
diffTrt.sserrvar = 1;
diffTrt.diffest = 2;
diffTrt.diffwpcov = 1;
diffTrt.diffrescor = 1;
diffTrt.expectedErrorId = '';
scenarioMap.trt_diff_event_time = diffTrt;

diffTrtRescor = diffTrt;
diffTrtRescor.name = 'trt_diff_rescor_event_time';
diffTrtRescor.diffwpcov = 2;
diffTrtRescor.diffrescor = 2;
scenarioMap.trt_diff_rescor_event_time = diffTrtRescor;

% Custom-prior coverage (W6 end-state follow-up). Every factory/DoD scenario
% above uses default priors, so the population-level prior substitution in the
% builders (the sprintf('%g',...) lines) is only ever pinned at default values.
% These scenarios pass a NON-default 'priors' struct through the real pipeline so
% the substitution is captured for each prior-consuming builder family: single
% (psyrat_build_stan_single), facet (psyrat_build_stan_facet), test-retest
% (psyrat_build_stan_trt), two-event difference (psyrat_build_stan_diff),
% subject-level difference (psyrat_build_stan_diff_sserr), and two-facet
% difference (psyrat_build_stan_diff_trt), and the four-event difference-of-
% differences (psyrat_build_stan_dodiff, priors.dod, exposed 2026-09-05). The
% sserr builder (analysis 6, psyrat_build_stan_sserr) exposes a single scale,
% priors.sserr.sig_trl, that no scenario here varies; its remaining constants
% stay fixed.
% Partial structs are sufficient: psyrat_computevarcomp merges them over the
% defaults (psyrat_merge_priors), so only the tested family needs values.
singlePriors = struct('single', struct('mu', 123.5, 'sig_u', 33, 'sig_e', 27));
trtPriors = struct('trt', struct('sig_id', 12, 'sig_occ', 6, 'sig_trl', 11, ...
    'sig_err', 22, 'sig_trlxid', 4, 'sig_occxid', 3, 'sig_trlxocc', 2));
diffPriors = struct('diff', struct('b', 7, 'b_sigma', 8, 'sd_id', 9, ...
    'sd_trl', 11, 'sd_id_sigma', 3, 'rescor_mu', 2, 'sd_rescor', 1.5));
% Two-facet difference builder (psyrat_build_stan_diff_trt) prior substitution.
% Non-default values for every diff_trt scale so each appears in the golden.
diffTrtPriors = struct('diff_trt', struct('b', 7, 'b_sigma', 8, 'sd_id', 9, ...
    'sd_trl', 11, 'sd_occ', 6, 'sd_tid', 4, 'sd_oid', 3, 'sd_to', 2));

scenarioMap.ic_base_priors = localPriorScenario('ic_base_priors', ...
    false, false, false, 1, 1, 1, singlePriors);
scenarioMap.ic_event_priors = localPriorScenario('ic_event_priors', ...
    false, true, false, 1, 1, 1, singlePriors);
scenarioMap.trt_time_priors = localPriorScenario('trt_time_priors', ...
    false, false, true, 1, 1, 1, trtPriors);
scenarioMap.ic_diff_event_priors = localPriorScenario('ic_diff_event_priors', ...
    false, true, false, 1, 2, 1, diffPriors);
scenarioMap.ic_diff_group_event_sserr_priors = localPriorScenario( ...
    'ic_diff_group_event_sserr_priors', true, true, false, 2, 2, 2, diffPriors);
% diffest=2 WITH time -> two-facet difference builder, non-concurrent
% (diffwpcov=1). Pins the priors.diff_trt substitution.
scenarioMap.trt_diff_event_time_priors = localPriorScenario( ...
    'trt_diff_event_time_priors', false, true, true, 1, 2, 1, diffTrtPriors);
% Four-event difference-of-differences (analysis 9) with a NON-default
% priors.dod block: every one of the five exposed scales differs from its
% default so each substitution appears in the golden (lkj deliberately
% non-integer, so the %g rendering of a fractional shape is pinned too). Built
% from the DoD scenario above rather than localPriorScenario, which has no
% isDoD/dodmap arms.
dodPriors = struct('dod', struct('b_cell', 7, 'b_sigma_cell', 3, 'sd_id', 4, ...
    'sd_trial', 5, 'lkj', 1.5));
dodp = dod;
dodp.name = 'dod_event_priors';
dodp.priors = dodPriors;
scenarioMap.dod_event_priors = dodp;

% Dynamic/conditional reliability (Rast & Clayson, analysis 11). The one-facet
% location-scale builder (psyrat_build_stan_dynrel) is KDIM-generic, so the
% emitted Stan text is identical for one and two dimensions; a single scenario
% pins the model source (and its priors.dynrel substitution). The harness builds
% a dimension-bearing table inline (isDynrel) and passes 'dynrel', 2.
dynrel = struct();
dynrel.name = 'dynrel_onefacet';
dynrel.includeGroup = false;
dynrel.includeEvent = false;
dynrel.includeTime = false;
dynrel.sserrvar = 1;
dynrel.diffest = 1;
dynrel.diffwpcov = 1;
dynrel.diffrescor = 1;
dynrel.expectedErrorId = '';
dynrel.isDynrel = true;
scenarioMap.dynrel_onefacet = dynrel;

% Gamma / scaled-chi-square ONE-facet dynamic/conditional reliability (analysis
% 11, family='gamma'). Same dynrel routing (dynrel=2, one facet, a dimension
% column) but the gamma builder (psyrat_build_stan_dynrel_chisq) swaps the
% observation likelihood: log-linked mean + person-varying log-nu, both with
% dimension slopes (b, b_nu), and a scaled chi-square. Reuses the isDynrel
% dimension table; the runner shifts meas strictly positive for the gamma family.
dynrelGamma = dynrel;
dynrelGamma.name = 'dynrel_gamma_onefacet';
dynrelGamma.family = 'gamma';
dynrelGamma.gammascale = 1; %movable design; pins log-nu
scenarioMap.dynrel_gamma_onefacet = dynrelGamma;

% Gamma / scaled-chi-square ONE-facet SUBJECT-LEVEL dynamic/conditional
% reliability (analysis 26, ic_dynrel_sserrvar, family='gamma'): dynrel=2 +
% sserrvar=2 routes to case 26. Case 26 shares its case block AND its builder with
% the group-level case 11 (psyrat_build_stan_dynrel_chisq) - the person-varying
% log-nu already gives each participant their own dispersion - so the emitted Stan
% is BYTE-IDENTICAL to dynrel_gamma_onefacet. The subject-level distinction is a
% downstream extraction (ind_bs/ind_nu are transformed parameters of that same
% model) plus a conditional calculator (psyrat_ssrel_dynrel_gamma). The distinct
% scenario pins the case-26 routing and its golden.
dynrelGammaSserr = dynrelGamma;
dynrelGammaSserr.name = 'dynrel_gamma_sserr_onefacet';
dynrelGammaSserr.sserrvar = 2;
dynrelGammaSserr.expectedAnalysis = 'ic_dynrel_sserrvar';
dynrelGammaSserr.gammascale = 1; %movable design; pins log-nu
scenarioMap.dynrel_gamma_sserr_onefacet = dynrelGammaSserr;

% Gamma LOCATION-SCALE one-facet dynamic/conditional reliability (analysis 11,
% family='gamma', gammascale=2). Same routing and the same standata as
% dynrel_gamma_onefacet; the builder swaps the scale submodel from log-nu to a
% direct log residual-SD model and derives nu = 2*(mu/sigma)^2, so b_sigma is a
% PURE residual-SD slope instead of the mean/dispersion blend b_nu carries. The
% two scenarios are pinned side by side deliberately: the pair of goldens is the
% clearest statement of what the option actually changes, and a diff between them
% should show ONLY the scale submodel and its priors - never the data block, the
% mean submodel, the trial facet or the likelihood.
dynrelGammaLs = dynrelGamma;
dynrelGammaLs.name = 'dynrel_gamma_ls_onefacet';
dynrelGammaLs.gammascale = 2;
scenarioMap.dynrel_gamma_ls_onefacet = dynrelGammaLs;

% Gamma LOCATION-SCALE one-facet SUBJECT-LEVEL dynamic reliability (analysis 26,
% gammascale=2). As with the log-nu pair, case 26 shares its case block AND its
% builder with case 11, so the emitted Stan is BYTE-IDENTICAL to
% dynrel_gamma_ls_onefacet; the subject-level distinction is the downstream
% extraction (ind_bs/ind_sd) plus a conditional calculator. The distinct scenario
% pins the case-26 routing under gammascale=2.
dynrelGammaLsSserr = dynrelGammaLs;
dynrelGammaLsSserr.name = 'dynrel_gamma_ls_sserr_onefacet';
dynrelGammaLsSserr.sserrvar = 2;
dynrelGammaLsSserr.expectedAnalysis = 'ic_dynrel_sserrvar';
scenarioMap.dynrel_gamma_ls_sserr_onefacet = dynrelGammaLsSserr;

% Dynamic difference-score reliability (Rast & Clayson, group-level; analysis
% 12). dynrel=2 + diffest=2 routes to psyrat_build_stan_diffdynrel (the two-event
% difference model with event-specific dimension slopes on the mean and per-event
% residual). The harness builds a two-event dimension-bearing table inline
% (isDynrelDiff) and passes dynrel=2 + diffest=2. KDIM-generic, so one snapshot
% pins both the 1- and 2-dimension Stan source.
dynreldiff = struct();
dynreldiff.name = 'dynrel_diff_onefacet';
dynreldiff.includeGroup = false;
dynreldiff.includeEvent = true;
dynreldiff.includeTime = false;
dynreldiff.sserrvar = 1;
dynreldiff.diffest = 2;
dynreldiff.diffwpcov = 1;
dynreldiff.diffrescor = 1;
dynreldiff.expectedErrorId = '';
dynreldiff.isDynrel = true;
dynreldiff.isDynrelDiff = true;
scenarioMap.dynrel_diff_onefacet = dynreldiff;

% Gamma / scaled-chi-square ONE-facet GROUP-LEVEL dynamic/conditional DIFFERENCE
% score (analysis 12, family='gamma'). Same case-12 routing (dynrel=2, diffest=2,
% sserrvar=1, non-concurrent) but the gamma builder
% (psyrat_build_stan_diffdynrel_chisq) swaps the observation likelihood: per-event
% log-linked mean + per-event person-varying log-nu, BOTH carrying the
% event-crossed dimension slopes (b_dim, b_nu_dim), and a scaled chi-square. This
% replaces the Gaussian per-event log-residual slopes (b_sigma_dim) entirely.
% Reuses the isDynrelDiff two-event dimension table; the runner shifts meas
% strictly positive for the gamma family. KDIM-generic, so one snapshot pins both
% the 1- and 2-dimension Stan source.
dynreldiffGamma = dynreldiff;
dynreldiffGamma.name = 'dynrel_gamma_diff_onefacet';
dynreldiffGamma.family = 'gamma';
dynreldiffGamma.gammascale = 1; %movable design; pins log-nu
scenarioMap.dynrel_gamma_diff_onefacet = dynreldiffGamma;

% Gamma LOCATION-SCALE one-facet GROUP-LEVEL dynamic/conditional DIFFERENCE score
% (analysis 12, gammascale=2). Increment 8b. Same case-12 routing as the log-nu
% scenario above; the builder swaps the log-nu submodel for a per-event log
% residual SD carrying its own event-crossed dimension slopes (b_sigma /
% b_sigma_dim), with nu DERIVED as 2*(mu/sigma)^2.
%
% The two scenarios are pinned side by side deliberately, exactly as the 14/27
% pair is: a diff between dynrel_gamma_diff_onefacet and this one should show ONLY
% the scale submodel and its priors - never the data block, the transformed data
% Xdim rebuild, the mean submodel, the cross-event trial block, or the likelihood.
%
% COUNT THE PRIORS CORRECTLY when using this as an acceptance criterion. FOUR
% prior lines differ (b_nu->b_sigma, b_nu_dim->b_sigma_dim, and sd_id[3] and
% sd_id[4] separately), which is THREE distinct priors. Two of those lines sit
% inside the 4-D person block, so "the person block is unchanged" is true of its
% STRUCTURE - dimension 4, the Cholesky factor, the r_id construction - but NOT of
% its two dispersion-SD priors. An earlier draft of this comment said "its two
% priors" and excluded the person block wholesale; that would flag a correct diff
% as a regression. Undercounting this diff is a repeat failure in this workstream
% (increment 8a did it too), so the numbers are spelled out here rather than left
% to inspection. That is the structural counterpart
% of the argument that psyrat_gamma_crosscov is reused unchanged: all eight of its
% arguments are mean-submodel quantities, so if the mean submodel is byte-identical
% the cross-event covariances cannot depend on the parameterization. See
% gamma_scale_submodel_decision.md section I.3.
dynreldiffGammaLs = dynreldiffGamma;
dynreldiffGammaLs.name = 'dynrel_gamma_ls_diff_onefacet';
dynreldiffGammaLs.gammascale = 2;
scenarioMap.dynrel_gamma_ls_diff_onefacet = dynreldiffGammaLs;

% Dynamic difference-score reliability (Rast & Clayson, SUBJECT-level; analysis
% 13). dynrel=2 + diffest=2 + sserrvar=2 routes to
% psyrat_build_stan_diffdynrel_sserr (the group-level dynamic difference model
% with the person random-effect block widened to the joint 4-coefficient
% (0 + event | p | ID) on both the location and the log-residual submodels). The
% harness reuses the two-event dimension-bearing table (isDynrelDiff) and passes
% sserrvar=2. KDIM-generic, so one snapshot pins both the 1- and 2-dimension Stan.
dynreldiffss = struct();
dynreldiffss.name = 'dynrel_diff_sserr_onefacet';
dynreldiffss.includeGroup = false;
dynreldiffss.includeEvent = true;
dynreldiffss.includeTime = false;
dynreldiffss.sserrvar = 2;
dynreldiffss.diffest = 2;
dynreldiffss.diffwpcov = 1;
dynreldiffss.diffrescor = 1;
dynreldiffss.expectedErrorId = '';
dynreldiffss.isDynrel = true;
dynreldiffss.isDynrelDiff = true;
scenarioMap.dynrel_diff_sserr_onefacet = dynreldiffss;

% Gamma LOCATION-SCALE one-facet SUBJECT-LEVEL dynamic/conditional DIFFERENCE score
% (analysis 13, family='gamma', gammascale=2). Same case-13 routing as the Gaussian
% scenario above (dynrel=2, diffest=2, sserrvar=2, non-concurrent).
%
% THE DELEGATION THIS PINS RUNS THE OPPOSITE WAY FROM THE GAUSSIAN PAIR, and that
% is the point of snapshotting it. In the Gaussian family the subject-level builder
% is a genuinely DIFFERENT program: dynrel_diff_onefacet has a 2-D person block and
% dynrel_diff_sserr_onefacet widens it to 4 and adds an er_var_ss generated
% quantity. Under gamma there is nothing to widen - the location-scale case-12
% model already declares sd_id = [mean_1, mean_2, logsigma_1, logsigma_2] with a
% 4x4 person Cholesky and already applies r_id[,3:4] to the log residual SD - so
% psyrat_build_stan_diffdynrel_sserr_chisq is a ONE-LINE DELEGATION and this
% golden must come out BYTE-IDENTICAL to dynrel_gamma_ls_diff_onefacet.
%
% ACCEPTANCE CRITERION, and it is stricter than the 12-vs-12 log-nu/LS comparison
% described in the block above (NOT stricter than the sibling subject-level
% scenarios, which go through the identical byte-identity check in
% TestStanModelSnapshots' pairs list): a diff against dynrel_gamma_ls_diff_onefacet
% must be EMPTY. Not "only the scale submodel" - empty. If any line differs,
% the delegation has been broken and the subject-level path is fitting a different
% model from the one its per-participant conversion assumes. TestStanModelSnapshots
% pins the identity pair explicitly for the same reason the 27<-14 and
% sserr_onefacet<-onefacet pairs are pinned.
%
% There is deliberately NO log-nu counterpart scenario. Analysis 13 is
% location-scale ONLY (guarded in psyrat_computevarcomp's option block), so a
% gammascale=1 request for this design errors rather than emitting Stan. That
% rejection is pinned by TestGammaScaleOption, not here.
dynreldiffssGammaLs = dynreldiffss;
dynreldiffssGammaLs.name = 'dynrel_gamma_ls_diff_sserr_onefacet';
dynreldiffssGammaLs.family = 'gamma';
dynreldiffssGammaLs.gammascale = 2;
scenarioMap.dynrel_gamma_ls_diff_sserr_onefacet = dynreldiffssGammaLs;

% Dynamic/conditional reliability, trial + occasion two-facet (Rast & Clayson;
% analysis 14). dynrel=2 + a time/occasion facet (ntime>0) + diffest=1 +
% sserrvar=1 routes to psyrat_build_stan_dynrel_trt (the test-retest crossed
% design extended with the joint location+scale person Cholesky, dimension slopes
% on the mean and log-residual, and a dimension-conditioned person-specific
% residual). The harness builds an occasion + dimension table (isDynrelTrt) and
% passes dynrel=2. KDIM-generic, so one snapshot pins both the 1- and 2-dimension
% Stan.
dynreltrt = struct();
dynreltrt.name = 'dynrel_trt_onefacet';
dynreltrt.includeGroup = false;
dynreltrt.includeEvent = false;
dynreltrt.includeTime = true;
dynreltrt.sserrvar = 1;
dynreltrt.diffest = 1;
dynreltrt.diffwpcov = 1;
dynreltrt.diffrescor = 1;
dynreltrt.expectedErrorId = '';
dynreltrt.isDynrel = true;
dynreltrt.isDynrelTrt = true;
scenarioMap.dynrel_trt_onefacet = dynreltrt;

% Gamma / scaled-chi-square trial + occasion TWO-facet dynamic/conditional
% reliability (analysis 14, family='gamma'). Same dynrel_trt routing (dynrel=2 +
% time facet) but the gamma builder (psyrat_build_stan_dynrel_trt_chisq) swaps the
% observation likelihood: log-linked crossed mean + person-varying log-nu, both
% with dimension slopes (b, b_nu), and a scaled chi-square. Reuses the isDynrelTrt
% occasion+dimension table; the runner shifts meas strictly positive for gamma.
dynreltrtGamma = dynreltrt;
dynreltrtGamma.name = 'dynrel_gamma_trt_onefacet';
dynreltrtGamma.family = 'gamma';
dynreltrtGamma.gammascale = 1; %movable design; pins log-nu
scenarioMap.dynrel_gamma_trt_onefacet = dynreltrtGamma;

% Gamma / scaled-chi-square trial + occasion two-facet SUBJECT-LEVEL
% dynamic/conditional reliability (analysis 27, ic_dynrel_sserrvar_trt,
% family='gamma'): dynrel=2 + a time facet + sserrvar=2 routes to case 27. Case 27
% shares its case block AND its builder with the group-level case 14
% (psyrat_build_stan_dynrel_trt_chisq) - the person-varying log-nu already gives
% each participant their own dispersion - so the emitted Stan is BYTE-IDENTICAL to
% dynrel_gamma_trt_onefacet. The subject-level distinction is a downstream
% extraction (ind_bs/ind_nu are transformed parameters of that same model) plus a
% conditional calculator (psyrat_ssrel_dynrel_trt_gamma). The distinct scenario
% pins the case-27 routing and its golden.
dynreltrtGammaSserr = dynreltrtGamma;
dynreltrtGammaSserr.name = 'dynrel_gamma_sserr_trt_onefacet';
dynreltrtGammaSserr.sserrvar = 2;
dynreltrtGammaSserr.expectedAnalysis = 'ic_dynrel_sserrvar_trt';
dynreltrtGammaSserr.gammascale = 1; %movable design; pins log-nu
scenarioMap.dynrel_gamma_sserr_trt_onefacet = dynreltrtGammaSserr;

% Gamma LOCATION-SCALE trial + occasion TWO-facet dynamic/conditional reliability
% (analysis 14, family='gamma', gammascale=2). Same routing and the same standata
% as dynrel_gamma_trt_onefacet; the builder swaps the scale submodel from log-nu
% to a direct log residual-SD model and derives nu = 2*(mu/sigma)^2, so b_sigma is
% a PURE residual-SD slope instead of the mean/dispersion blend b_nu carries. The
% two scenarios are pinned side by side deliberately: a diff between them should
% show ONLY the scale submodel and its two priors - never the data block, the mean
% submodel, ANY OF THE FIVE CROSSED FACETS, or the likelihood. That last exclusion
% is what this increment adds over the one-facet pair, and it is the structural
% counterpart of the estimand-invariance argument: the six observed-score signal
% components are built from the log-mean submodel alone.
dynreltrtGammaLs = dynreltrtGamma;
dynreltrtGammaLs.name = 'dynrel_gamma_ls_trt_onefacet';
dynreltrtGammaLs.gammascale = 2;
scenarioMap.dynrel_gamma_ls_trt_onefacet = dynreltrtGammaLs;

% Gamma LOCATION-SCALE trial + occasion two-facet SUBJECT-LEVEL dynamic
% reliability (analysis 27, gammascale=2). As with the log-nu pair, case 27 shares
% its case block AND its builder with case 14, so the emitted Stan is
% BYTE-IDENTICAL to dynrel_gamma_ls_trt_onefacet; the subject-level distinction is
% the downstream extraction (ind_sd is a transformed parameter of that same model)
% plus a conditional calculator (psyrat_ssrel_dynrel_trt_gamma_ls). The distinct
% scenario pins the case-27 routing under gammascale=2.
dynreltrtGammaLsSserr = dynreltrtGammaLs;
dynreltrtGammaLsSserr.name = 'dynrel_gamma_ls_sserr_trt_onefacet';
dynreltrtGammaLsSserr.sserrvar = 2;
dynreltrtGammaLsSserr.expectedAnalysis = 'ic_dynrel_sserrvar_trt';
scenarioMap.dynrel_gamma_ls_sserr_trt_onefacet = dynreltrtGammaLsSserr;

% Dynamic/conditional reliability, trial + occasion two-facet GROUP-LEVEL
% difference score (Rast & Clayson; analysis 15). dynrel=2 + a time/occasion
% facet (ntime>0) + two events + diffest=2 + sserrvar=1 routes to
% psyrat_build_stan_diffdynrel_trt (the static two-facet difference model, case
% 10, extended with event-specific dimension slopes on the mean and per-event
% log-residual). The harness builds a two-event + occasion + dimension table
% (isDynrelDiffTrt) and passes dynrel=2 + diffest=2. KDIM-generic, so one
% snapshot pins both the 1- and 2-dimension Stan.
dynreldifftrt = struct();
dynreldifftrt.name = 'dynrel_diff_trt_onefacet';
dynreldifftrt.includeGroup = false;
dynreldifftrt.includeEvent = true;
dynreldifftrt.includeTime = true;
dynreldifftrt.sserrvar = 1;
dynreldifftrt.diffest = 2;
dynreldifftrt.diffwpcov = 1;
dynreldifftrt.diffrescor = 1;
dynreldifftrt.expectedErrorId = '';
dynreldifftrt.isDynrel = true;
dynreldifftrt.isDynrelDiffTrt = true;
scenarioMap.dynrel_diff_trt_onefacet = dynreldifftrt;

% Dynamic/conditional reliability, trial + occasion two-facet SUBJECT-LEVEL
% difference score (Rast & Clayson; analysis 16). dynrel=2 + a time/occasion facet
% (ntime>0) + two events + diffest=2 + sserrvar=2 routes to
% psyrat_build_stan_diffdynrel_sserr_trt (the group-level case-15 two-facet
% difference model with the case-13 joint 4-coefficient person block on the
% location AND scale submodels: per-person event-specific scale random effects and
% per-subject er_var_ss). The harness reuses the two-event + occasion + dimension
% table (isDynrelDiffSserrTrt -> localBuildDynrelDiffTrtTable) and passes dynrel=2 +
% diffest=2 + sserrvar=2. KDIM-generic, so one snapshot pins both the 1- and
% 2-dimension Stan source.
dynreldiffsstrt = struct();
dynreldiffsstrt.name = 'dynrel_diff_sserr_trt_onefacet';
dynreldiffsstrt.includeGroup = false;
dynreldiffsstrt.includeEvent = true;
dynreldiffsstrt.includeTime = true;
dynreldiffsstrt.sserrvar = 2;
dynreldiffsstrt.diffest = 2;
dynreldiffsstrt.diffwpcov = 1;
dynreldiffsstrt.diffrescor = 1;
dynreldiffsstrt.expectedErrorId = '';
dynreldiffsstrt.isDynrel = true;
dynreldiffsstrt.isDynrelDiffSserrTrt = true;
scenarioMap.dynrel_diff_sserr_trt_onefacet = dynreldiffsstrt;

% Dynamic/conditional reliability, one-facet group-level CONCURRENT difference
% (Rast & Clayson; analysis 17). dynrel=2 + diffest=2 + sserrvar=1 + diffwpcov=2
% (=> diffrescor=2) routes to psyrat_build_stan_diffdynrel_rescor (the static
% concurrent one-facet model, case 7, extended with event-specific dimension slopes
% on the mean and per-event log-residual and a per-row residual Cholesky). The
% harness reuses the balanced two-event dimension table (isDynrelDiffRescor ->
% localBuildDynrelDiffTable, which has equal event counts so events pair) and passes
% dynrel=2 + diffwpcov=2. KDIM-generic, so one snapshot pins both 1- and 2-dimension
% Stan source.
dynreldiffrescor = struct();
dynreldiffrescor.name = 'dynrel_diff_rescor_onefacet';
dynreldiffrescor.includeGroup = false;
dynreldiffrescor.includeEvent = true;
dynreldiffrescor.includeTime = false;
dynreldiffrescor.sserrvar = 1;
dynreldiffrescor.diffest = 2;
dynreldiffrescor.diffwpcov = 2;
dynreldiffrescor.diffrescor = 2;
dynreldiffrescor.expectedErrorId = '';
dynreldiffrescor.isDynrel = true;
dynreldiffrescor.isDynrelDiffRescor = true;
scenarioMap.dynrel_diff_rescor_onefacet = dynreldiffrescor;

% Dynamic/conditional reliability, trial + occasion two-facet group-level
% CONCURRENT difference (Rast & Clayson; analysis 18). dynrel=2 + a time/occasion
% facet + diffest=2 + sserrvar=1 + diffwpcov=2 routes to
% psyrat_build_stan_diffdynrel_trt_rescor (the static concurrent two-facet model,
% case 10, extended with the dimension slopes + per-row residual Cholesky). The
% harness reuses the balanced two-event + occasion + dimension table
% (isDynrelDiffRescorTrt -> localBuildDynrelDiffTrtTable) and passes dynrel=2 +
% diffwpcov=2. KDIM-generic, so one snapshot pins both 1- and 2-dimension Stan.
dynreldiffrescortrt = struct();
dynreldiffrescortrt.name = 'dynrel_diff_rescor_trt_onefacet';
dynreldiffrescortrt.includeGroup = false;
dynreldiffrescortrt.includeEvent = true;
dynreldiffrescortrt.includeTime = true;
dynreldiffrescortrt.sserrvar = 1;
dynreldiffrescortrt.diffest = 2;
dynreldiffrescortrt.diffwpcov = 2;
dynreldiffrescortrt.diffrescor = 2;
dynreldiffrescortrt.expectedErrorId = '';
dynreldiffrescortrt.isDynrel = true;
dynreldiffrescortrt.isDynrelDiffRescorTrt = true;
scenarioMap.dynrel_diff_rescor_trt_onefacet = dynreldiffrescortrt;

% Dynamic/conditional reliability, one-facet SUBJECT-LEVEL CONCURRENT difference
% (Rast & Clayson; analysis 19). dynrel=2 + diffest=2 + sserrvar=2 + diffwpcov=2
% routes to psyrat_build_stan_diffdynrel_sserr_rescor (the case-17 one-facet
% concurrent model with the case-13 joint 4-coefficient person block: per-person
% scale random effect + population residual correlation). The harness reuses the
% balanced two-event dimension table (isDynrelDiffSserrRescor ->
% localBuildDynrelDiffTable) and passes sserrvar=2 + diffwpcov=2.
dynreldiffssrescor = struct();
dynreldiffssrescor.name = 'dynrel_diff_sserr_rescor_onefacet';
dynreldiffssrescor.includeGroup = false;
dynreldiffssrescor.includeEvent = true;
dynreldiffssrescor.includeTime = false;
dynreldiffssrescor.sserrvar = 2;
dynreldiffssrescor.diffest = 2;
dynreldiffssrescor.diffwpcov = 2;
dynreldiffssrescor.diffrescor = 2;
dynreldiffssrescor.expectedErrorId = '';
dynreldiffssrescor.isDynrel = true;
dynreldiffssrescor.isDynrelDiffSserrRescor = true;
scenarioMap.dynrel_diff_sserr_rescor_onefacet = dynreldiffssrescor;

% The GAMMA twin of the scenario above (analysis 19, family='gamma',
% gammascale=2): the MODULAR CUT-COPULA stage-1 margins model
% (psyrat_build_stan_diffdynrel_sserr_margins_chisq).
%
% THIS SNAPSHOT IS LOAD-BEARING BEYOND THE USUAL "did the text move" CHECK. The
% two families differ not only in likelihood but in INFERENCE FRAMEWORK: the
% Gaussian model estimates the residual correlation inside HMC, this one must
% not. The golden therefore pins the ABSENCE of the copula from stage 1, and
% TestModularCutCopulaStanBoundary asserts that absence directly by token, so a
% future change that "optimized" the two stages back together would fail loudly
% rather than silently changing the estimand.
%
% gammascale is pinned explicitly even though 2 is this design's default and its
% only legal value: pinning it keeps the golden's provenance readable and keeps
% the scenario honest if the default ever moves.
dynreldiffssrescorGamma = dynreldiffssrescor;
dynreldiffssrescorGamma.name = 'dynrel_gamma_ls_diff_sserr_rescor_onefacet';
dynreldiffssrescorGamma.family = 'gamma';
dynreldiffssrescorGamma.gammascale = 2;
scenarioMap.dynrel_gamma_ls_diff_sserr_rescor_onefacet = dynreldiffssrescorGamma;

% Dynamic/conditional reliability, trial + occasion two-facet SUBJECT-LEVEL
% CONCURRENT difference (Rast & Clayson; analysis 20). dynrel=2 + a time/occasion
% facet + diffest=2 + sserrvar=2 + diffwpcov=2 routes to
% psyrat_build_stan_diffdynrel_sserr_trt_rescor (the case-18 two-facet concurrent
% model with the case-16 joint 4-coefficient person block). The harness reuses the
% balanced two-event + occasion + dimension table (isDynrelDiffSserrRescorTrt ->
% localBuildDynrelDiffTrtTable) and passes sserrvar=2 + diffwpcov=2.
dynreldiffssrescortrt = struct();
dynreldiffssrescortrt.name = 'dynrel_diff_sserr_rescor_trt_onefacet';
dynreldiffssrescortrt.includeGroup = false;
dynreldiffssrescortrt.includeEvent = true;
dynreldiffssrescortrt.includeTime = true;
dynreldiffssrescortrt.sserrvar = 2;
dynreldiffssrescortrt.diffest = 2;
dynreldiffssrescortrt.diffwpcov = 2;
dynreldiffssrescortrt.diffrescor = 2;
dynreldiffssrescortrt.expectedErrorId = '';
dynreldiffssrescortrt.isDynrel = true;
dynreldiffssrescortrt.isDynrelDiffSserrRescorTrt = true;
scenarioMap.dynrel_diff_sserr_rescor_trt_onefacet = dynreldiffssrescortrt;

% The GAMMA twin of the scenario above (analysis 20): the two-facet MODULAR
% CUT-COPULA stage-1 margins model
% (psyrat_build_stan_diffdynrel_sserr_trt_margins_chisq). Same reasoning as the
% one-facet gamma scenario; this one additionally pins that all five crossed
% facets load on the MEAN and none on the residual SD, which is the property the
% person-specific residual depends on.
dynreldiffssrescortrtGamma = dynreldiffssrescortrt;
dynreldiffssrescortrtGamma.name = 'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet';
dynreldiffssrescortrtGamma.family = 'gamma';
dynreldiffssrescortrtGamma.gammascale = 2;
scenarioMap.dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet = dynreldiffssrescortrtGamma;

% Dynamic/conditional reliability, one-facet SUBJECT-LEVEL CONCURRENT difference
% with a PER-SUBJECT residual correlation (Rast & Clayson; analysis 21, Stage 4e).
% dynrel=2 + diffest=2 + sserrvar=2 + diffwpcov=2 + ssrescor=2 routes to
% psyrat_build_stan_diffdynrel_sserr_rho (the case-19 model with the single
% population Lrescor replaced by the case-8 per-person rho hierarchy
% rho[s]=tanh(rescor_mu+sd_rescor*z)). The harness reuses the balanced two-event
% dimension table (isDynrelDiffSsRho -> localBuildDynrelDiffTable) and forwards
% ssrescor=2 alongside sserrvar=2 + diffwpcov=2.
dynreldiffssrho = struct();
dynreldiffssrho.name = 'dynrel_diff_sserr_ssrescor_onefacet';
dynreldiffssrho.includeGroup = false;
dynreldiffssrho.includeEvent = true;
dynreldiffssrho.includeTime = false;
dynreldiffssrho.sserrvar = 2;
dynreldiffssrho.diffest = 2;
dynreldiffssrho.diffwpcov = 2;
dynreldiffssrho.diffrescor = 2;
dynreldiffssrho.ssrescor = 2;
dynreldiffssrho.expectedErrorId = '';
dynreldiffssrho.isDynrel = true;
dynreldiffssrho.isDynrelDiffSsRho = true;
scenarioMap.dynrel_diff_sserr_ssrescor_onefacet = dynreldiffssrho;

% Dynamic/conditional reliability, trial + occasion two-facet SUBJECT-LEVEL
% CONCURRENT difference with a PER-SUBJECT residual correlation (analysis 22,
% Stage 4e). Adds a time/occasion facet; ssrescor=2 routes to
% psyrat_build_stan_diffdynrel_sserr_trt_rho (the case-20 model with the same
% per-person rho hierarchy). The harness reuses the balanced two-event + occasion +
% dimension table (isDynrelDiffSsRhoTrt -> localBuildDynrelDiffTrtTable) and forwards
% ssrescor=2.
dynreldiffssrhotrt = struct();
dynreldiffssrhotrt.name = 'dynrel_diff_sserr_ssrescor_trt_onefacet';
dynreldiffssrhotrt.includeGroup = false;
dynreldiffssrhotrt.includeEvent = true;
dynreldiffssrhotrt.includeTime = true;
dynreldiffssrhotrt.sserrvar = 2;
dynreldiffssrhotrt.diffest = 2;
dynreldiffssrhotrt.diffwpcov = 2;
dynreldiffssrhotrt.diffrescor = 2;
dynreldiffssrhotrt.ssrescor = 2;
dynreldiffssrhotrt.expectedErrorId = '';
dynreldiffssrhotrt.isDynrel = true;
dynreldiffssrhotrt.isDynrelDiffSsRhoTrt = true;
scenarioMap.dynrel_diff_sserr_ssrescor_trt_onefacet = dynreldiffssrhotrt;

% Nonparallel data-splits, single occasion (Rocha Table 4; analysis case 23).
% splits=3 with a per-split items weight column and NO group/event/time routes to
% psyrat_build_stan_splits_iins (the weighted observed-design p x (i:s) model: an
% explicit person x split interaction plus an n_i-weighted per-row residual
% sig_err ./ sqrt(weight)). The harness builds an id/meas/weight table inline
% (isSplits -> localBuildSplitsTable) with real n_i variation and passes splits=3.
% The n_i variation is required: the splits=3 preflight errors on all-equal n_i
% (localValidateSplits), and that preflight runs before the emit hook, so an
% all-equal weight would abort capture before any Stan is generated.
splits = struct();
splits.name = 'splits_iins_single';
splits.includeGroup = false;
splits.includeEvent = false;
splits.includeTime = false;
splits.sserrvar = 1;
splits.diffest = 1;
splits.diffwpcov = 1;
splits.diffrescor = 1;
splits.expectedErrorId = '';
splits.isSplits = true;
splits.splits = 3;
scenarioMap.splits_iins_single = splits;

% Nonparallel data-splits, test-retest (Rocha Table 5; analysis case 24). Adds a
% time/occasion facet (ntime>0) so splits=3 routes to
% psyrat_build_stan_splits_iins_trt (the weighted observed-design p x o x (i:s)
% model: an explicit person x occasion x split interaction plus the same
% n_i-weighted per-row residual). The harness reuses the id/meas/weight table with
% a time column added (isSplitsTrt -> localBuildSplitsTrtTable) and passes splits=3.
splitstrt = struct();
splitstrt.name = 'splits_iins_trt';
splitstrt.includeGroup = false;
splitstrt.includeEvent = false;
splitstrt.includeTime = true;
splitstrt.sserrvar = 1;
splitstrt.diffest = 1;
splitstrt.diffwpcov = 1;
splitstrt.diffrescor = 1;
splitstrt.expectedErrorId = '';
splitstrt.isSplits = true;
splitstrt.isSplitsTrt = true;
splitstrt.splits = 3;
scenarioMap.splits_iins_trt = splitstrt;
end

function cfg = localPriorScenario(name, includeGroup, includeEvent, ...
    includeTime, sserrvar, diffest, diffwpcov, priors)
%Build a custom-prior snapshot scenario struct. The data table is built by the
%integration factory's makeScenarioTable (no isDoD flag). diffrescor follows
%diffwpcov, matching the factory and the psyrat_computevarcomp flag coupling.
cfg = struct();
cfg.name = name;
cfg.includeGroup = includeGroup;
cfg.includeEvent = includeEvent;
cfg.includeTime = includeTime;
cfg.sserrvar = sserrvar;
cfg.diffest = diffest;
cfg.diffwpcov = diffwpcov;
cfg.diffrescor = diffwpcov;
cfg.expectedErrorId = '';
cfg.priors = priors;
end
