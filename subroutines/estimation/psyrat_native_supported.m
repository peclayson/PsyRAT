function engines = psyrat_native_supported(key)
%Return which native estimation engines implement a given design.
%
%engines = psyrat_native_supported(key)
%
%Input:
% key - a design identifier set at the psyrat_run_stan call site (see
%  psyrat_nativespec in psyrat_computevarcomp). Examples: 'sing' (one-facet
%  internal consistency), 'trt' (test-retest), 'splits_iins'.
%
%Called with NO argument, returns the full cellstr of design keys the registry
%recognizes (the canonical key list). TestNativeEngine iterates this list to
%assert the registry agrees with the psyrat_native_hmc / psyrat_native_fitlme
%dispatch switches, guarding against the three places a key is spelled (the
%psyrat_run_stan case site, this registry, and the engine switches) drifting apart.
%
%Output:
% engines - cellstr subset of {'hmc','fitlme'} naming the native engines that
%  can fit this design. Empty {} means no native engine implements it yet, so
%  the design is CmdStan-only.
%
%This registry is the SINGLE place that records native-engine coverage. It is
%consulted by psyrat_resolve_engine to (a) build the auto cascade and (b)
%validate a forced-engine request. Add a key here as each native fitter lands;
%until then a design stays CmdStan-only and the cascade/forced paths report it
%with a clear message.

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

if nargin == 0
    %Canonical list of every design key the switch below recognizes. Keep in sync
    %with the cases; the TestNativeEngine consistency test verifies each key here
    %resolves to a non-empty engine set and is dispatchable by that engine.
    engines = {'sing','sing_facet','splits','splits_trt','trt','sserr', ...
        'dynrel','diff','diff_rescor','diff_dynrel','diff_dynrel_sserr', ...
        'dynrel_trt','dod','diff_dynrel_trt','diff_dynrel_sserr_trt', ...
        'diff_dynrel_rescor','diff_dynrel_trt_rescor','diff_trt', ...
        'diff_trt_rescor','diff_dynrel_sserr_rescor', ...
        'diff_dynrel_sserr_trt_rescor','diff_dynrel_sserr_rho', ...
        'diff_dynrel_sserr_trt_rho','diff_sserr'};
    return;
end

switch char(key)
    case {'sing','sing_facet','splits','splits_trt'}
        % Single-residual (plain-Gaussian) designs fitlme can fit: one-facet
        % internal consistency (case 1 sing; cases 2-4 group/event sing_facet) and
        % the nonparallel-splits observed designs (case 23 splits = single
        % occasion; case 24 splits_trt = test-retest, weighted). These have NO
        % native HMC engine, so fitlme is their only non-CmdStan option. NOTE:
        % case 6 (sserr) is per-subject location-scale and is NOT fitlme-doable
        % (Phase 2 / HMC); plain test-retest (case 5 trt) has BOTH engines (below).
        engines = {'fitlme'};
    case 'trt'
        % Plain test-retest crossed variance components (case 5): a plain-Gaussian
        % six-component crossed random-effects ANOVA. This is the ONE design with
        % TWO native engines. The auto-cascade prefers CmdStan -> native HMC
        % (psyrat_build_stan_trt, the faithful Bayesian fit) -> native fitlme (REML
        % + parametric bootstrap; weak when sig_occ is poorly identified at small
        % occasion counts). HMC is listed first to document the cascade preference;
        % both remain reachable by explicit force (engine 3 = hmc, engine 2 =
        % fitlme). psyrat_resolve_engine checks hmc before fitlme regardless of
        % list order.
        engines = {'hmc','fitlme'};
    case 'sserr'
        % Subject-level error variance (case 6): the base one-facet location-scale
        % model (person location + per-subject log-residual via a 2-D Cholesky,
        % crossed trial main effect; no dimension predictors). A per-subject
        % heteroscedastic residual is not fitlme-doable, so this is HMC-only
        % (psyrat_native_hmc, matching psyrat_build_stan_sserr).
        engines = {'hmc'};
    case 'dynrel'
        % One-facet dynamic/conditional reliability (Rast & Clayson, case 11):
        % the sserr location-scale model EXTENDED with fixed dimension slopes on
        % the mean and log-residual. Also HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_dynrel); shares one log-posterior with sserr.
        engines = {'hmc'};
    case 'diff'
        % Two-event difference-score location-scale model, no residual covariance
        % (case 7). A univariate Gaussian over the stacked event data with an
        % event-indicator design and two correlated 2-D random effects (person,
        % trial). Per-event residual SDs and correlated REs are not fitlme-doable,
        % so HMC-only (psyrat_native_hmc, matching psyrat_build_stan_diff).
        engines = {'hmc'};
    case 'diff_rescor'
        % Two-event difference-score model WITH residual covariance (case 7,
        % diffrescor=2). A BIVARIATE (multi_normal_cholesky) outcome over the
        % paired event responses with constant per-event residual SDs and a
        % residual correlation (Lrescor), plus the same two correlated 2-D
        % person/trial random effects. A multivariate residual is not fitlme-
        % doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diff_rescor).
        engines = {'hmc'};
    case 'diff_dynrel'
        % Group-level dynamic/conditional reliability of a two-event difference
        % score (case 12). The no-rescor difference model EXTENDED with
        % event-specific fixed dimension slopes on both the mean and the per-event
        % log-residual (exactly as dynrel extends sserr). The location-scale
        % per-event residual is not fitlme-doable, so HMC-only (psyrat_native_hmc,
        % matching psyrat_build_stan_diffdynrel).
        engines = {'hmc'};
    case 'diff_dynrel_sserr'
        % Subject-level dynamic difference score (case 13). The group-level
        % dynamic difference with the person random-effect block widened to a 4-D
        % joint correlated block (location + log-residual per event) so each
        % participant has their own per-event residual (a 4x4 LKJ Cholesky). Not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_sserr).
        engines = {'hmc'};
    case 'dynrel_trt'
        % Two-facet (trial + occasion) dynamic reliability, non-difference (case
        % 14). The one-facet dynrel locscale model with five scalar crossed-facet
        % mean random effects (occ/trl/trlxid/occxid/trlxocc). The location-scale
        % per-person residual is not fitlme-doable, so HMC-only (psyrat_native_hmc,
        % matching psyrat_build_stan_dynrel_trt).
        engines = {'hmc'};
    case 'dod'
        % Task-free four-event difference-of-differences (case 9). A 4-cell
        % location-scale model with person AND trial random effects on both the
        % location and the scale (an 8-D joint block per grouping factor via an
        % 8x8 LKJ Cholesky). The per-cell heteroscedastic residual is not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_dodiff).
        engines = {'hmc'};
    case 'diff_dynrel_trt'
        % Group-level two-facet (trial + occasion) dynamic difference (case 15).
        % The one-facet group-level dynamic difference with six crossed 2x2
        % cross-condition blocks (id/trl/occ/tid/oid/to) + event-by-dimension
        % slopes. The per-event location-scale residual is not fitlme-doable, so
        % HMC-only (psyrat_native_hmc, matching psyrat_build_stan_diffdynrel_trt).
        engines = {'hmc'};
    case 'diff_dynrel_sserr_trt'
        % Subject-level two-facet dynamic difference (case 16). The group-level
        % two-facet dynamic difference with the person factor widened to a 4-D
        % joint location+scale block (a 4x4 LKJ Cholesky) so each participant has
        % their own per-event residual; the five crossed factors stay 2-D. Not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_sserr_trt).
        engines = {'hmc'};
    case 'diff_dynrel_rescor'
        % Concurrent one-facet group-level dynamic difference (case 17). The
        % bivariate-residual difference model (case 7*) with per-event dimension
        % slopes on the mean and the per-event log-residual SD (so the bivariate
        % residual Cholesky is per-row). A multivariate residual is not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_rescor).
        engines = {'hmc'};
    case 'diff_dynrel_trt_rescor'
        % Concurrent two-facet group-level dynamic difference (case 18). The
        % concurrent one-facet model (case 17) with the person+trial 2-D blocks
        % generalized to the six crossed cross-condition 2-D factors (id/trl/occ/
        % tid/oid/to); same per-row bivariate residual + dimension slopes. Not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_trt_rescor).
        engines = {'hmc'};
    case {'diff_trt','diff_trt_rescor'}
        % Static two-facet (test-retest) two-event difference (case 10). The
        % KDIM = 0 reduction of the dynamic two-facet difference (cases 15/18):
        % the same six crossed 2x2 cross-condition blocks (id/trl/occ/tid/oid/to)
        % with NO dimension slopes. 'diff_trt' is the non-concurrent sub-model
        % (stacked-long, group-level per-event residual); 'diff_trt_rescor' is the
        % concurrent sub-model (paired bivariate residual + rescor). Both reuse the
        % case-15/18 engines (local_hmc_diffdynrel_trt[_rescor]) at KDIM = 0. The
        % per-event location-scale / multivariate residual is not fitlme-doable, so
        % HMC-only (psyrat_native_hmc, matching psyrat_build_stan_diff_trt[_rescor]).
        engines = {'hmc'};
    case 'diff_dynrel_sserr_rescor'
        % Concurrent one-facet subject-level dynamic difference (case 19). The
        % concurrent one-facet group model (case 17) with the person block widened
        % to a 4-D joint location+scale block (a 4x4 LKJ Cholesky), so each
        % participant has their own per-event residual SDs; the residual
        % correlation stays population. Not fitlme-doable, so HMC-only
        % (psyrat_native_hmc, matching psyrat_build_stan_diffdynrel_sserr_rescor).
        engines = {'hmc'};
    case 'diff_dynrel_sserr_trt_rescor'
        % Concurrent two-facet subject-level dynamic difference (case 20). The
        % concurrent one-facet subject-level model (case 19) with the single trial
        % 2-D block generalized to the five crossed cross-condition 2-D factors
        % (trl/occ/tid/oid/to); the person factor stays the 4-D joint block. Not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_sserr_trt_rescor).
        engines = {'hmc'};
    case 'diff_dynrel_sserr_rho'
        % Concurrent one-facet subject-level dynamic difference with a PER-SUBJECT
        % residual correlation (case 21). The case-19 model with the population
        % residual correlation replaced by a hierarchical per-participant
        % rho[s] = tanh(rescor_mu + sd_rescor*z_rescor[s]). Not fitlme-doable, so
        % HMC-only (psyrat_native_hmc, matching psyrat_build_stan_diffdynrel_sserr_rho).
        engines = {'hmc'};
    case 'diff_dynrel_sserr_trt_rho'
        % Concurrent two-facet subject-level dynamic difference with a PER-SUBJECT
        % residual correlation (case 22). The case-20 model with the population
        % residual correlation replaced by the case-21 per-participant rho[s]
        % hierarchy. Not fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diffdynrel_sserr_trt_rho).
        engines = {'hmc'};
    case 'diff_sserr'
        % Static subject-level two-event difference (case 8). Three separate 2-D
        % blocks (person location, trial, per-subject scale) + a per-subject
        % residual correlation, over paired (bivariate) and unpaired (univariate)
        % observations. The per-subject heteroscedastic residual is not
        % fitlme-doable, so HMC-only (psyrat_native_hmc, matching
        % psyrat_build_stan_diff_sserr).
        engines = {'hmc'};
    otherwise
        engines = {};
end
end
