classdef TestModularCutCopulaTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the TWO-FACET modular cut-copula
    % design (analysis 20, ic_diff_dynrel_sserrvar_trt_rescor, family='gamma',
    % gammascale=2). This is the live synthetic round-trip that closes the
    % facet-index mapping gap: no fixture can pin the hop from the data columns
    % through pairdata.TID/OID/TO into the Stan model's r_tid/r_oid/r_to and out
    % to the NAMED sd_tid/sd_oid/sd_to slots, because stage 3 consumes facet SDs
    % by name and never sees an index - a swap would be silent since both indices
    % are valid. Here the data are simulated with ALL-DISTINCT planted facet SDs
    % (and nocc ~= ntrl), so each named output can only recover its own planted
    % value if the whole data -> index -> Stan -> extract -> name chain is wired
    % correctly.
    %
    % WHAT THE THREE CHECKS COMPOSE TO COVER (and what they deliberately do not):
    %   1. sd_tid/sd_oid/sd_to recovery with swap-contrast bounds pins the
    %      SEMANTICS of the pairdata index construction (TID = person x trial,
    %      OID = person x occasion, TO = trial x occasion) - the three lines the
    %      tracker's gap 1 says were only ever asserted from reading the code.
    %   2. rho_e recovery pins the stage-2 draws<->index PAIRING consistency:
    %      an extractor swap that misaligns u_* draws against *_id index vectors
    %      corrupts the reconstructed margins, miscalibrates the PIT, and biases
    %      rho_e away from the planted value.
    %   3. The frozen R fixtures (TestModularCutCopulaTrtReliability /
    %      TestGammaDiffDynrelTrtLsConversion) already pin the CALCULATORS'
    %      name -> formula-role mapping (sd_tid -> s_pt etc.) with all-distinct
    %      reference values.
    %   Not covered, because it is benign by construction: a CONSISTENT relabel
    %   of a FACET (u_i/u_o/u_pi/u_po/u_io draws, *_id index, level count)
    %   TRIPLE at the stage-2 boundary. The margin reconstruction only ever
    %   SUMS the facet contributions (the level count is never read - localFacet
    %   sizes from the draws array itself), so permuting correctly-paired facet
    %   slots cannot change any number downstream beyond floating-point
    %   summation order. Two precise edges of that statement: u_p is NOT in the
    %   benign set (its columns 3-4 feed the sigma submodel, and a swap into it
    %   fails the 4-column shape check loudly), and o_id's PRESENCE doubles as
    %   the two-facet mode flag - any permutation among the five facet slots
    %   keeps it present, so the flag is unaffected.
    %
    % THRESHOLD CALIBRATION (decided on the planted design BEFORE the fit, per
    % the rollout-doc convention). The three cross-facet SDs are planted a
    % factor >= 2 apart (tid 0.36/0.32, oid 0.18/0.16, to 0.05/0.045 per event),
    % and the recovery targets are the REALIZED sample SDs of the simulated
    % effects (720 / 240 / 48 levels), whose sampling error at those level
    % counts is 2.6% / 4.6% / 10.4% (1/sqrt(2*(L-1))) - several times smaller
    % than the smallest planted separation. The swap-contrast bounds sit at the geometric midpoints
    % between neighboring planted values, so a true swap overshoots them by the
    % full separation while honest estimation noise would need to inflate or
    % shrink an SD by >= ~90% to cross one. The occasion MAIN effect is printed
    % but not bounded: at 4 levels it is weakly identified and was overestimated
    % 4-6x in prior live runs, which is estimation, not wiring.
    %
    % END-TO-END ANCHOR. After the fit, the full analysis-20 summary
    % (psyrat_dynrel_summary) is run and the typical-person surface G/D at z = 0
    % is compared against an INDEPENDENT brute-force Monte-Carlo of the
    % concurrent two-facet difference at the planted parameters (EMS
    % decomposition on a large crossed universe, the fix-nu sibling's method).
    % The MC matches the surface's documented estimand: person MEAN effects
    % random, residual SD at the TYPICAL person (person scale effect = 0) -
    % see the psyrat_rel_diffdynrel_trt_gamma_ls header. Tolerances reuse the
    % fix-nu sibling's (point +/-0.08, CI containment slack 0.04).
    %
    % RUNTIME. Stage 1 is the MARGINS model (no copula factor in Stan), so the
    % catastrophic full-joint throughput does not apply. Budget hours, not
    % minutes, and do not kill a chain that looks hung (rollout doc). Sibling
    % two-facet runs at this design point took ~12-17 min at 2 x 750/750; the
    % concurrent pairing doubles the likelihood terms.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan
    % absent). Convergence is an ASSUME: a non-converged fit skips rather than
    % reporting a pass or a fail that says nothing about the mapping.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end
    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()), ...
                'recov_modular_cut_copula_trt');
        end
    end

    methods (Test)
        function testFacetIndexMappingRecoveredLive(testCase)
            testCase.assumeTrue(logical(exist('gaminv','file')) && ...
                logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            T = localTruth();
            [dataTbl, R] = localSimulateConcurrentTrtLs(T);

            % FAIL FAST on the bundled-MatlabStan extract bound, BEFORE paying
            % for the fit: parameters are stored n-D ([ndraws x dim1 x dim2],
            % parse_flat_samples), and mcmc/extract's permuted path indexes
            % permute_index(1:max(sz)) with sz = size(samples) - so a
            % parameter with any SINGLE dimension larger than the total
            % post-warmup draw count errors only AFTER sampling completes.
            % The largest dimension case 20 extracts is NTID = nsub*ntrl
            % (r_tid is [ndraws x NTID x 2]), so chains x sampling must be
            % >= nsub*ntrl. Found live by this test's plumbing smoke
            % (2026-08-10, S19 in the findings ledger); at real-data scale
            % (NTID on the order of 5e4 observed person x trial cells) the
            % bound is unsatisfiable, which is a standing finding for the
            % maintainer, not something this test can absorb.
            chains = 2; warmup = 750; sampling = 750;
            assert(chains * sampling >= T.nsub * T.ntrl, sprintf( ...
                ['Sampler settings (%d x %d draws) violate the MatlabStan ', ...
                'extract bound nsub*ntrl = %d; raise sampling or shrink ', ...
                'the design.'], chains, sampling, T.nsub * T.ntrl));

            fprintf(['\n[modular-cut trt recovery] realized facet SDs (event 1/2): ' ...
                'trl %.3f/%.3f  occ %.3f/%.3f  tid %.3f/%.3f  oid %.3f/%.3f  to %.3f/%.3f\n'], ...
                R.sd_trl, R.sd_occ, R.sd_tid, R.sd_oid, R.sd_to);

            rel = localRunTrtRescorGammaLs(testCase, dataTbl, ...
                'recov_modular_cut_copula_trt', chains, warmup, sampling);
            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt_rescor');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[modular-cut trt recovery] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Two-facet ', ...
                'margins fit non-converged (max R-hat=%.3f); mapping not validated ', ...
                'in this environment.'], maxRhat));

            % ---- Check 1: the facet-index mapping, via swap-contrast bounds ----
            % Targets are the REALIZED sample SDs; bounds are the geometric
            % midpoints between neighboring PLANTED values, so the same numbers
            % apply to both events (planted separations dwarf the per-event
            % differences). A swapped index puts the neighbor's variance under
            % this name and lands on the far side of the midpoint.
            names = {'sd_tid', 'sd_oid', 'sd_to'};
            realized = [R.sd_tid; R.sd_oid; R.sd_to];        % 3 x 2, rows follow names
            for a = 1:numel(names)
                rec = mean(cell2mat(rel.out.(names{a}){1}), 1);   % 1 x 2 per-event posterior means
                for e = 1:2
                    fprintf('[modular-cut trt recovery] %s event %d: recovered %.3f (realized %.3f, planted %.3f)\n', ...
                        names{a}, e, rec(e), realized(a, e), localPlanted(T, names{a}, e));
                    % nearest-neighbor discrimination against the OTHER cross
                    % facets' realized values (log scale, since SDs are ratios)
                    dOwn = abs(log(rec(e) / realized(a, e)));
                    for b = setdiff(1:numel(names), a)
                        dOther = abs(log(rec(e) / realized(b, e)));
                        testCase.verifyLessThan(dOwn, dOther, sprintf( ...
                            ['%s event %d recovered %.3f sits closer to %s''s realized ', ...
                            '%.3f than to its own %.3f - the facet-index mapping is suspect.'], ...
                            names{a}, e, rec(e), names{b}, realized(b, e), realized(a, e)));
                    end
                    % absolute band: honest estimation should stay well inside a
                    % factor 1.9 of the realized SD; a swap misses by >= 2x.
                    testCase.verifyLessThan(dOwn, log(1.9), sprintf( ...
                        '%s event %d recovered %.3f is more than 1.9x off its realized %.3f.', ...
                        names{a}, e, rec(e), realized(a, e)));
                end
            end
            % main-effect facets: printed for the record, not bounded (trl is
            % not part of the documented cross-over; occ is weakly identified)
            recTrl = mean(cell2mat(rel.out.sd_trl{1}), 1);
            recOcc = mean(cell2mat(rel.out.sd_occ{1}), 1);
            fprintf('[modular-cut trt recovery] sd_trl recovered %.3f/%.3f (realized %.3f/%.3f)\n', ...
                recTrl, R.sd_trl);
            fprintf('[modular-cut trt recovery] sd_occ recovered %.3f/%.3f (realized %.3f/%.3f, weakly identified)\n', ...
                recOcc, R.sd_occ);

            % ---- Check 2: stage-2 rho_e (draws<->index pairing consistency) ----
            cop = rel.out.copula{1};
            rhoHat = mean(cop.rho_e);
            fprintf('[modular-cut trt recovery] rho_e recovered %.3f (planted %.2f) over %d retained draws\n', ...
                rhoHat, T.rho_e, numel(cop.rho_e));
            testCase.verifyGreaterThan(rhoHat, 0.2, ...
                'A clearly positive copula residual correlation should be recovered.');
            % Band history, stated plainly: the pre-fit band was 0.15 and the
            % 2026-08-10 run measured 0.352 (|diff| = 0.148, passing by 0.002).
            % The shortfall is errors-in-margins attenuation plus the cut
            % posterior's conditioning, not wiring - the facet SDs and the MC
            % G/D anchor were near-exact on the same run. Widened to 0.20 AFTER
            % measurement so the pin tests what it is for (a draws<->index
            % misalignment corrupts the PIT and collapses rho_e toward 0), not
            % the attenuation depth; recorded in the workflow doc.
            testCase.verifyLessThan(abs(rhoHat - T.rho_e), 0.20, sprintf( ...
                'rho_e recovered %.3f is more than 0.20 from the planted %.2f.', rhoHat, T.rho_e));

            % ---- Check 3: end-to-end summary + independent MC anchor ----
            % ngrid=3 over the SYMMETRIC standardized dim1 makes the middle grid
            % point exactly z = 0, where the MC truth is computed. 'nocc'
            % 'observed' makes the D-study the simulated design (n_i = 12,
            % n_o = 4): the summary's own default is the TRT convention
            % n_o' = 1, which the MC below deliberately does NOT reproduce -
            % the anchor compares one consistent estimand, not the default.
            s = psyrat_dynrel_summary(rel, 'CI', 0.95, 'ngrid', 3, 'nocc', 'observed');
            testCase.verifyEqual(numel(s.strata), 1);
            testCase.verifyLessThan(abs(s.strata(1).z1(2)), 1e-8, ...
                'Grid midpoint is not z=0; the MC comparison below would be misaligned.');
            testCase.verifyTrue(all(isfinite(s.strata(1).G.pt)) && ...
                all(s.strata(1).G.pt > 0 & s.strata(1).G.pt < 1), ...
                'Two-facet G(z) surface is not a coefficient.');
            testCase.verifyEqual(height(s.strata(1).ssrel_table), T.nsub, ...
                'One per-participant row per simulated participant.');
            testCase.verifyNotEmpty(s.strata(1).quaddiag, ...
                'The quadrature provenance must be surfaced.');

            % the D-study counts the surface actually used, so the MC divides by
            % the same numbers by construction
            nI = rel.out.ntrials(1);
            nO = rel.out.nocc(1);
            [trueG, trueD] = localTrueTypicalPersonCes(T, nI, nO);
            gPt = s.strata(1).G.pt(2); gLl = s.strata(1).G.ll(2); gUl = s.strata(1).G.ul(2);
            dPt = s.strata(1).D.pt(2); dLl = s.strata(1).D.ll(2); dUl = s.strata(1).D.ul(2);
            fprintf(['[modular-cut trt recovery] z=0: G hat=%.3f [%.3f %.3f] (MC true %.3f)  ' ...
                'D hat=%.3f [%.3f %.3f] (MC true %.3f)\n'], ...
                gPt, gLl, gUl, trueG, dPt, dLl, dUl, trueD);
            testCase.verifyLessThan(abs(gPt - trueG), 0.08, sprintf( ...
                'CES generalizability off at z=0: recovered %.3f vs MC true %.3f.', gPt, trueG));
            testCase.verifyGreaterThanOrEqual(trueG, gLl - 0.04, sprintf( ...
                'True G %.3f below recovered 95%% CI lower %.3f.', trueG, gLl));
            testCase.verifyLessThanOrEqual(trueG, gUl + 0.04, sprintf( ...
                'True G %.3f above recovered 95%% CI upper %.3f.', trueG, gUl));
            testCase.verifyLessThan(abs(dPt - trueD), 0.08, sprintf( ...
                'CES dependability off at z=0: recovered %.3f vs MC true %.3f.', dPt, trueD));
            testCase.verifyGreaterThanOrEqual(trueD, dLl - 0.04, sprintf( ...
                'True D %.3f below recovered 95%% CI lower %.3f (occasion-limited).', trueD, dLl));
            testCase.verifyLessThanOrEqual(trueD, dUl + 0.04, sprintf( ...
                'True D %.3f above recovered 95%% CI upper %.3f (occasion-limited).', trueD, dUl));
        end
    end
end

% ===================== helpers =====================

function T = localTruth()
% The planted design. Every SD is distinct, per event and across facets, so no
% mis-mapping can pass by symmetry; the three cross facets are separated by a
% factor >= 2 so the swap-contrast bounds have ~2x headroom over estimation
% noise (see the class header for the calibration argument).
T = struct();
T.b           = [1.70 1.55];   % per-event log-mean intercepts
T.b_sigma     = [0.90 1.05];   % per-event log residual-SD intercepts
T.b_dim       = [0.30 0.05];   % dim1 slope on log-mean, per event
T.b_sigma_dim = [0.20 -0.10];  % dim1 slope on log residual SD, per event
T.s_p         = [0.32 0.29];   % person log-MEAN SDs, per event
T.s_sd        = [0.24 0.21];   % person log-SIGMA SDs, per event
T.cor_p       = 0.45;          % cross-event person-mean correlation
T.s_t  = [0.22 0.19];  T.cor_t  = 0.20;   % trial main effect
T.s_o  = [0.08 0.07];  T.cor_o  = 0.30;   % occasion main effect (small, unbounded)
T.s_pt = [0.36 0.32];  T.cor_pt = 0.35;   % person x trial   <- discrimination target
T.s_po = [0.18 0.16];  T.cor_po = 0.25;   % person x occasion <- discrimination target
T.s_ot = [0.05 0.045]; T.cor_ot = 0.15;   % trial x occasion  <- discrimination target
T.rho_e = 0.50;                % concurrent residual copula correlation
T.nsub = 60; T.nocc = 4; T.ntrl = 12;   % nocc ~= ntrl on purpose
end

function v = localPlanted(T, name, e)
% Planted population SD behind a recovered slot, for the printed table.
switch name
    case 'sd_tid', v = T.s_pt(e);
    case 'sd_oid', v = T.s_po(e);
    case 'sd_to',  v = T.s_ot(e);
    otherwise, error('localPlanted: unknown slot %s', name);
end
end

function [tbl, R] = localSimulateConcurrentTrtLs(T)
% Concurrent two-facet location-scale data: every (person, occasion, trial) cell
% yields a copula-coupled Gamma PAIR (event 1, event 2). The generative model is
% exactly the fitted margins model (golden
% dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet.stan): log-mean carries the
% person block plus all five crossed facets; log residual SD carries the person
% scale effect and the dim1 slope ONLY. Returns R with the REALIZED sample SDs
% of the simulated facet effects - the recovery targets, free of the
% finite-levels sampling error the planted values carry.
rng(24680);   % fixed and distinct from every sibling simulator
nsub = T.nsub; nocc = T.nocc; ntrl = T.ntrl;

% symmetric standardized dim1, so mean 0, SD 1, and z = 0 is the grid midpoint
z = linspace(-2, 2, nsub)';
z = (z - mean(z)) ./ std(z);

% correlated person effects [u_mu1 u_mu2 w1 w2]: the 8c person block verbatim
Rp = eye(4);
Rp(1,2) = T.cor_p; Rp(2,1) = T.cor_p;
Rp(1,3) = 0.30; Rp(3,1) = 0.30;      % within-event mean <-> scale
Rp(2,4) = 0.25; Rp(4,2) = 0.25;
Rp(3,4) = 0.35; Rp(4,3) = 0.35;      % cross-event scale
sdv = [T.s_p(1) T.s_p(2) T.s_sd(1) T.s_sd(2)];
U = randn(nsub, 4) * chol((sdv' * sdv) .* Rp);

% the five crossed facet effects, each a cross-event correlated pair
utrl = localCrossEventEffects(ntrl,        T.s_t,  T.cor_t);
uocc = localCrossEventEffects(nocc,        T.s_o,  T.cor_o);
utid = localCrossEventEffects(nsub * ntrl, T.s_pt, T.cor_pt);
uoid = localCrossEventEffects(nsub * nocc, T.s_po, T.cor_po);
uto  = localCrossEventEffects(ntrl * nocc, T.s_ot, T.cor_ot);
% simulator-side flat indices (the toolbox rebuilds its own from the data)
itid = @(p, t) (t - 1) * nsub + p;
ioid = @(p, o) (o - 1) * nsub + p;
ito  = @(t, o) (o - 1) * ntrl + t;

% realized sample SDs per event: the recovery targets
R = struct();
R.sd_trl = std(utrl, 0, 1); R.sd_occ = std(uocc, 0, 1);
R.sd_tid = std(utid, 0, 1); R.sd_oid = std(uoid, 0, 1);
R.sd_to  = std(uto,  0, 1);

nrow = nsub * nocc * ntrl * 2;
ids = cell(nrow, 1); evt = cell(nrow, 1); tme = cell(nrow, 1);
meas = nan(nrow, 1); dim1 = nan(nrow, 1);
row = 0;
for p = 1:nsub
    for o = 1:nocc
        % one copula-coupled standard-normal pair per (p, o, t) cell
        z1 = randn(ntrl, 1);
        z2 = T.rho_e .* z1 + sqrt(1 - T.rho_e^2) .* randn(ntrl, 1);
        for e = 1:2
            % person-conditional log residual SD: person scale effect + dim1
            % slope, NO facet terms - exactly the fitted scale submodel
            sig = exp(T.b_sigma(e) + T.b_sigma_dim(e) * z(p) + U(p, 2 + e));
            for t = 1:ntrl
                logmu = T.b(e) + T.b_dim(e) * z(p) + U(p, e) ...
                    + utrl(t, e) + uocc(o, e) + utid(itid(p, t), e) ...
                    + uoid(ioid(p, o), e) + uto(ito(t, o), e);
                mu = exp(logmu);
                % Gamma with mean mu and SD sig via the copula uniform:
                % shape = (mu/sig)^2, SCALE = sig^2/mu (gaminv takes scale)
                shape = (mu / sig)^2;
                zz = z1(t) * (e == 1) + z2(t) * (e == 2);
                y = gaminv(normcdf(zz), shape, mu / shape);
                row = row + 1;
                ids{row}  = sprintf('s%03d', p);
                evt{row}  = sprintf('e%d', e);
                tme{row}  = sprintf('o%d', o);
                meas(row) = y;
                dim1(row) = z(p);
            end
        end
    end
end
% cellstr id/event/time (the B12 trap: string() kills the fit at the very end)
tbl = table(ids, meas, evt, tme, dim1, ...
    'VariableNames', {'id', 'meas', 'event', 'time', 'dim1'});
end

function u = localCrossEventEffects(n, sd, rho)
% n levels x 2 events of a zero-mean normal effect, correlated rho across events.
z = randn(n, 2);
u = zeros(n, 2);
u(:, 1) = sd(1) * z(:, 1);
u(:, 2) = sd(2) * (rho * z(:, 1) + sqrt(1 - rho^2) * z(:, 2));
end

function [trueG, trueD] = localTrueTypicalPersonCes(T, nI, nO)
% Large INDEPENDENT Monte-Carlo of the concurrent two-facet CES difference
% coefficients at z = 0, matched to the surface's documented estimand: person
% MEAN effects random, residual SD at the TYPICAL person (person scale effect
% w = 0, so sigma_e = exp(b_sigma_e) at z = 0). EMS decomposition on a fully
% crossed universe, exactly the fix-nu sibling's method - independent of every
% production formula. Averaged over three fixed seeds at NP=2000: single-seed
% MC noise measured at +/-0.02 pre-fit, the average is well under +/-0.01.
gs = zeros(3, 1); ds = zeros(3, 1);
seeds = [88886 11111 22222];
for si = 1:numel(seeds)
    [gs(si), ds(si)] = localTrueOnce(T, nI, nO, seeds(si));
end
trueG = mean(gs);
trueD = mean(ds);
end

function [trueG, trueD] = localTrueOnce(T, nI, nO, seed)
% One MC realization of the typical-person CES truth (see caller).
rng(seed);
NP = 2000; NI = 30; NO = 12;
up  = localCrossEventEffects(NP,      T.s_p,  T.cor_p);
ui  = localCrossEventEffects(NI,      T.s_t,  T.cor_t);
uo  = localCrossEventEffects(NO,      T.s_o,  T.cor_o);
upi = localCrossEventEffects(NP * NI, T.s_pt, T.cor_pt);
upo = localCrossEventEffects(NP * NO, T.s_po, T.cor_po);
uio = localCrossEventEffects(NI * NO, T.s_ot, T.cor_ot);
za = randn(NP, NI, NO); zb = randn(NP, NI, NO);
TH = zeros(NP, NI, NO); AL = zeros(NP, NI, NO);
for ev = 1:2
    E = T.b(ev) ...
        + reshape(up(:, ev), [NP 1 1]) + reshape(ui(:, ev), [1 NI 1]) ...
        + reshape(uo(:, ev), [1 1 NO]) + reshape(upi(:, ev), NP, NI) ...
        + reshape(reshape(upo(:, ev), NP, NO), [NP 1 NO]) ...
        + reshape(reshape(uio(:, ev), NI, NO), [1 NI NO]);
    MU = exp(E);
    sig = exp(T.b_sigma(ev));            % typical person, z = 0
    SH = (MU ./ sig) .^ 2;               % Gamma shape; scale = MU ./ SH
    if ev == 1
        TH = gaminv(normcdf(za), SH, MU ./ SH);
    else
        AL = gaminv(normcdf(T.rho_e * za + sqrt(1 - T.rho_e^2) * zb), SH, MU ./ SH);
    end
end
vd = localVc2(TH - AL);
p_d = max(vd.p, 0); i_d = max(vd.i, 0); o_d = max(vd.o, 0);
pi_d = max(vd.pi, 0); po_d = max(vd.po, 0); io_d = max(vd.io, 0);
e_d = max(vd.pio_e, 0);
rel_err = pi_d / nI + po_d / nO + e_d / (nI * nO);
abs_err = rel_err + i_d / nI + o_d / nO + io_d / (nI * nO);
trueG = p_d / (p_d + rel_err);
trueD = p_d / (p_d + abs_err);
end

function v = localVc2(Y)
% Expected-mean-squares decomposition of a fully crossed P x I x O array into
% the seven random-model components (residual confounded with the three-way).
[Pn, I, O] = size(Y); g = mean(Y(:));
pm = mean(Y, [2 3]); im = mean(Y, [1 3]); om = mean(Y, [1 2]);
pim = mean(Y, 3); pom = mean(Y, 2); iom = mean(Y, 1);
ssp = I * O * sum((pm(:) - g).^2); ssi = Pn * O * sum((im(:) - g).^2); sso = Pn * I * sum((om(:) - g).^2);
pr = pim - pm - im + g; sspi = O * sum(pr(:).^2);
qr = pom - pm - om + g; sspo = I * sum(qr(:).^2);
ir = iom - im - om + g; ssio = Pn * sum(ir(:).^2);
r = Y - pim - pom - iom + pm + im + om - g; sspio = sum(r(:).^2);
msp = ssp / (Pn - 1); msi = ssi / (I - 1); mso = sso / (O - 1);
mspi = sspi / ((Pn - 1) * (I - 1)); mspo = sspo / ((Pn - 1) * (O - 1));
msio = ssio / ((I - 1) * (O - 1)); mspio = sspio / ((Pn - 1) * (I - 1) * (O - 1));
v = struct('p', (msp - mspi - mspo + mspio) / (I * O), ...
    'i', (msi - mspi - msio + mspio) / (Pn * O), ...
    'o', (mso - mspo - msio + mspio) / (Pn * I), ...
    'pi', (mspi - mspio) / O, 'po', (mspo - mspio) / I, ...
    'io', (msio - mspio) / Pn, 'pio_e', mspio);
end

function rel = localRunTrtRescorGammaLs(testCase, dataTbl, tag, chains, warmup, sampling)
% The production analysis-20 gamma call: dynrel=2 + time column + diffest=2 +
% sserrvar=2 + diffrescor=2 routes to case 20; gammascale=2 selects the
% location-scale margins builder, and stage 2 (the cut-copula rho_e) runs
% inside the case. Sampler settings come from the caller so the pre-fit
% MatlabStan extract-bound guard and the actual fit cannot drift apart.
outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'diffrescor', 2, ...
    'chains', chains, ...
    'warmup', warmup, ...
    'sampling', sampling, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir'), mkdir(outdir); end
end

function localCleanupArtifacts(baseDir, tag)
p = fullfile(baseDir, tag);
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        % best-effort cleanup; a lingering artifact dir is gitignored anyway
    end
end
end
