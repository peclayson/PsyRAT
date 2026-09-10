classdef TestGammaDiffDynrelSserrLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE Gamma /
    % scaled-chi-square SUBJECT-LEVEL dynamic DIFFERENCE score (analysis 13,
    % dynrel = 2 + diffest = 2 + sserrvar = 2 + family = 'gamma' +
    % gammascale = 2). Increment 8c.
    %
    % WHAT THIS CHECKS THAT THE OFFLINE SUITE CANNOT. The accuracy suite is
    % regression-only - its oracle duplicates the production formulas - and
    % TestGammaDiffDynrelSserrLsConversion pins the conversion arithmetic against
    % a closed form written from section J. Neither shows that the FITTED model
    % recovers the generative truth. This does: data are simulated from the
    % location-scale generative model with known parameters, fit with CmdStan, and
    % the per-participant coefficients compared against their own planted truth.
    %
    % THE ESTIMAND UNDER TEST IS THE PER-PARTICIPANT TABLE, not the surface. Case
    % 12 already covers the surface, and this design emits the identical Stan
    % model (the builder is a one-line delegation), so re-testing the surface would
    % re-test case 12. The whole difference between 12 and 13 is
    % psyrat_ssrel_diffdynrel_gamma_ls plus the psyrat_dynrel_summary branch that
    % calls it, so the run is routed THROUGH psyrat_dynrel_summary rather than
    % calling the calculator directly - that is the new code.
    %
    % THE DESIGN POINT IS SIZED FOR THE PERSON BLOCK. Increment 8b lost a full
    % 80-minute fit to this: all four dimension slopes recovered essentially
    % perfectly while the surface sat above truth, which localized the fault to the
    % person block rather than the estimand (LKJ(1) shrinks cor_p; s_p absorbs the
    % slack when an event carries a large residual against few trials). Here slots
    % 3-4 of that same 4-D person block are no longer nuisance parameters - they
    % ARE the estimand, because they are what makes participants differ. So the
    % 8b-corrected sizing (250 x 40) is carried over rather than reduced, and
    % s_sd is planted well away from zero: at s_sd -> 0 every participant collapses
    % to the same coefficient and the per-participant table would reproduce the
    % group value while looking perfectly healthy.
    %
    % EXPECT PER-SUBJECT POINTS ON BOTH SIDES OF THE GROUP CURVE. Since S17-FIX
    % the per-participant error POOLS the mean-surface person x trial term
    % exp(2*alpha(z))*K with the participant's own conditional variance,
    % matching the group surface's construction. AT MATCHED TRIAL COUNTS -
    % which this simulation plants deliberately, every participant sharing the
    % curve's n' - the pooled term is then common to both sides and drops out
    % of the error-term comparison exactly. What remains, and remains
    % two-signed: the curve marginalizes the conditional term over the person
    % residual-SD population (exp(2*log_sigma + 2*s_sd^2)) while each
    % participant uses their own realized w_s. (At UNMATCHED counts the pooled
    % term does not fully cancel - it is divided by different n on the two
    % sides - but that remainder is sign-aligned with the counts mechanism, so
    % the two-signed conclusion survives in production plots too.) The test
    % asserts the gap in both directions and its monotone gradient in the
    % participant's own scale effect, so a future change to the subject-level
    % estimand fails here loudly.
    %
    % Budget ~80 minutes. CmdStan-gated (skips Incomplete, not Failed, when
    % CmdStan/MatlabStan are absent).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            base = localCmdStanOutputDir(testCase.projectRoot());
            if isfolder(base)
                rmdir(base, 's');
            end
        end
    end

    methods (Test)

        function testPerParticipantCoefficientsRecovered(testCase)
            T = localTruth();
            rng(20260806, 'twister');   % data generation only; estimation seed is fixed at 12345
            [tbl, W, z] = localSimulateGammaDiffDynrelSserrLs(T);

            rel = localRunGammaDiffDynrelSserrLs(testCase, tbl, ...
                'recov_gamma_diffdynrel_sserr_ls');

            % Route through the summary: that is the new branch, and it is where a
            % misrouted analysis-13 run would silently produce the GROUP result
            % with no per-participant table at all.
            summ = psyrat_dynrel_summary(rel, 'CI', 0.95, 'ngrid', 9);
            testCase.assertTrue(isfield(summ.strata(1), 'ssrel_table'), ...
                ['No per-participant table was produced. A case-13 gamma run ' ...
                'that falls into the group-level branch fails exactly this way, ' ...
                'and nothing else about it looks wrong.']);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), T.nsub);

            % ---- per-participant truth: the POOLED estimand (S17-FIX) ----
            trueG = localTruePerSubjectGPooled(T, z, W);

            covered = (trueG >= tab.gen_ll) & (trueG <= tab.gen_ul);
            frac = mean(covered);
            fprintf('\n--- case-13 LS recovery: per-participant G ---\n');
            fprintf('  95%% CI coverage across %d participants: %.3f\n', T.nsub, frac);
            fprintf('  true  G: min %.4f  median %.4f  max %.4f\n', ...
                min(trueG), median(trueG), max(trueG));
            fprintf('  est   G: min %.4f  median %.4f  max %.4f\n', ...
                min(tab.gen_pt), median(tab.gen_pt), max(tab.gen_pt));

            % Coverage, not exactness: these are 250 separate posterior intervals,
            % so a handful of misses is expected even under a correct model. The
            % floor is deliberately below 0.95 for that reason, but well above what
            % a systematically wrong conversion would produce.
            testCase.verifyGreaterThan(frac, 0.85, ...
                sprintf(['Only %.1f%% of participants'' true G fell inside their ' ...
                'credible interval; a conversion error shows up here as a ' ...
                'systematic miss, not a few tail misses.'], 100*frac));

            % ---- the SPREAD across participants is the point of the design ----
            % A converter that ignored the per-person scale effects would return
            % nearly identical coefficients for everyone and still pass a coverage
            % test built on wide intervals. Pin the dispersion and the ordering.
            testCase.verifyGreaterThan(std(tab.gen_pt), 0.02, ...
                ['Per-participant coefficients are nearly constant. Either the ' ...
                'person scale effects are not reaching the coefficient, or ' ...
                's_sd was planted too small.']);
            r = localSpearman(tab.gen_pt, trueG);
            fprintf('  Spearman(est, true) across participants: %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.80, ...
                ['The per-participant ORDERING must be recovered, not just the ' ...
                'level - that ordering is what a subject-level readout is for.']);

            % The participants with the largest planted residuals must come out
            % least reliable. Checks the SIGN of the mapping, which a transposed or
            % mis-indexed er_ss matrix would break while preserving spread.
            %
            % MEASURE IT CONTROLLING FOR z, AND CALIBRATE AGAINST TRUTH. The first
            % version of this check used the raw Spearman against the planted scale
            % effects with a hard -0.5 floor, and it FAILED at -0.440 on a run whose
            % coverage (0.968) and ordering (0.929) were both healthy. The threshold
            % was the defect, not the fit: the marginal statistic is diluted by z,
            % which drives most of G's variation here (true G vs z = +0.74), so even
            % the PLANTED TRUTH only reaches -0.646 - a -0.5 floor left almost no
            % room for estimation noise. Within z bands the truth is -0.88 to -0.97,
            % which is the relationship actually being claimed.
            %
            % So: residualize both against z, and require the estimate to retain
            % most of the truth's own statistic on the same data. Self-calibrating,
            % so it cannot pass vacuously and does not encode a guessed constant.
            ownResid = W(:,1) + W(:,2);
            spEst  = localBandedSpearman(tab.gen_pt, ownResid, z);
            spTrue = localBandedSpearman(trueG,      ownResid, z);
            fprintf('  within-z-band Spearman(G, own scale): est %.3f, true %.3f\n', ...
                spEst, spTrue);
            testCase.assertLessThan(spTrue, -0.5, ...
                ['the planted design must itself show the relationship, or the ' ...
                'check below is vacuous']);
            testCase.verifyLessThan(spEst, 0.5*spTrue, ...
                ['Larger own residual must mean lower reliability: the estimated ' ...
                'within-z relationship must retain at least half the strength of ' ...
                'the planted one.']);

            % ---- dependability <= generalizability, per participant ----
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9), ...
                'Absolute error cannot be smaller than relative error.');

            % ---- the gap against the group curve is NOT one-signed ----
            % (Two-signed before S17-FIX and still two-signed after it: at the
            % matched trial counts this simulation plants, the pooled person x
            % trial term is common to the curve and the table and drops out of
            % the error-term comparison, so the mechanism is unchanged.)
            % An earlier version of this test picked the participant nearest z = 0
            % and asserted their coefficient exceeded the curve. That assertion is
            % FALSE and would have failed here for a reason that is not a bug: the
            % curve marginalizes the conditional term over the person residual-SD
            % population (exp(2*log_sigma + 2*s_sd^2)) while a participant uses
            % their own exp(2*(log_sigma + w_s)), so any participant whose own
            % scale effect exceeds s_sd^2 sits BELOW. With s_sd = 0.24/0.21 planted
            % here, a large share of the sample does. Losing an 80-minute fit to a
            % wrong assertion is the trap increment 8b already fell into once.
            %
            % What IS testable is that the sign of the gap tracks each
            % participant's own scale effect. Every participant has the same
            % trial counts in this simulation, so the replication mechanism is
            % neutralized and the residual mechanism is isolated.
            grp = psyrat_rel_diffdynrel_gamma_ls( ...
                'b',            localDrawMat(rel,'b'), ...
                'b_sigma',      localDrawMat(rel,'b_sigma'), ...
                'b_dim',        localDrawMat(rel,'b_dim'), ...
                'b_sigma_dim',  localDrawMat(rel,'b_sigma_dim'), ...
                'sd_id',        localDrawMat(rel,'sd_id'), ...
                'sd_trl',       localDrawMat(rel,'sd_trl'), ...
                'cor_id',       localDrawMat(rel,'cor_id'), ...
                'cor_i',        localDrawMat(rel,'cor_i'), ...
                'ndim',1, 'z1',z(:)', 'obs',[T.ntrl T.ntrl], 'CI',0.95);
            gap = tab.gen_pt - grp.G.pt(:);
            ownScale = W(:,1) + W(:,2);     % planted per-participant scale effects
            fprintf('  gap: %.1f%% of participants above the curve, %.1f%% below\n', ...
                100*mean(gap > 0), 100*mean(gap < 0));
            fprintf('  gap: Spearman(gap, own scale effect) = %.3f\n', ...
                localSpearman(gap, ownScale));

            testCase.verifyTrue(any(gap > 0) && any(gap < 0), ...
                ['The per-subject-minus-curve gap must occur in BOTH ' ...
                'directions at this design point. If every participant landed ' ...
                'on one side, either s_sd collapsed or the curve and the table ' ...
                'stopped being different estimands.']);
            % Calibrated against the same statistic on the planted truth rather
            % than a guessed constant, for the reason recorded on the own-residual
            % check above: a hardcoded floor here came in at -0.821 against -0.80,
            % which is too thin to survive a different seed.
            spGapEst  = localSpearman(gap, ownScale);
            spGapTrue = localSpearman(trueG - grp.G.pt(:), ownScale);
            fprintf('  gap: Spearman(gap, own scale) true %.3f\n', spGapTrue);
            testCase.assertLessThan(spGapTrue, -0.5, ...
                'the planted design must itself show the own-scale gradient');
            testCase.verifyLessThan(spGapEst, 0.5*spGapTrue, ...
                ['The gap must DECREASE in the participant''s own scale ' ...
                'effect - that monotone relationship is what a change to the ' ...
                'subject-level estimand would break.']);

            % The signed S17-FIX comparison: against the PRE-FIX
            % (conditional-only) estimand, pooling each participant's share of
            % the nonnegative mean-surface term can only LOWER their
            % coefficient. Keeping the un-pooled oracle alive pins the
            % direction of what S17-FIX changed.
            trueGunpooled = localTruePerSubjectG(T, z, W);
            testCase.verifyTrue(all(trueG <= trueGunpooled + 1e-12), ...
                ['S17-FIX: pooling the mean-surface term must never RAISE a ' ...
                'participant''s coefficient - the direction is signed, as a ' ...
                'weak inequality.']);

            % ---- the person scale SDs, which are what make participants differ ----
            sdid = localDrawMat(rel, 'sd_id');
            localReportParam(testCase, sdid(:,3), T.s_sd(1), 's_sd event 1');
            localReportParam(testCase, sdid(:,4), T.s_sd(2), 's_sd event 2');
        end

    end
end

% =========================================================================
function T = localTruth()
% The audited design point, carried over from increment 8b's CORRECTED sizing.
% Changing any of these invalidates the header's argument; re-run the offline
% sweep before touching them.
T = struct();
T.b           = [1.75 1.60];
T.b_sigma     = [0.90 1.05];
T.b_dim       = [0.30 0.05];
T.b_sigma_dim = [0.20 -0.10];
T.s_p         = [0.32 0.29];   % person log-MEAN SDs, per event
T.s_sd        = [0.24 0.21];   % person log-SIGMA SDs, per event - the estimand
T.s_i         = [0.11 0.13];   % trial log-mean SDs, per event
T.cor_p       = 0.45;          % cross-event person-mean correlation
T.cor_i       = 0.20;          % cross-event trial correlation
T.nsub        = 250;
T.ntrl        = 40;            % per event; n' = [40 40]
end

% =========================================================================
function [tbl, W, z] = localSimulateGammaDiffDynrelSserrLs(T)
% Two-event non-concurrent data from the location-scale generative model.
% Returns the per-participant scale effects W (nsub x 2) and dimension values z
% as well, because the per-participant truth is conditional on them - that is
% what distinguishes this from the case-12 simulation.
%
% The events are NON-CONCURRENT, so each event gets its OWN trial index and the
% conditional draws are independent given (mu, sigma). That is the assumption the
% estimand rests on, so it is built into the simulation rather than assumed.
nsub = T.nsub; ntrl = T.ntrl;

z = linspace(-2, 2, nsub)';
z = (z - mean(z)) ./ std(z);

% correlated person effects: [u_mu1 u_mu2 w1 w2]
R = eye(4);
R(1,2) = T.cor_p; R(2,1) = T.cor_p;
R(1,3) = 0.30; R(3,1) = 0.30;      % within-event mean <-> scale
R(2,4) = 0.25; R(4,2) = 0.25;
R(3,4) = 0.35; R(4,3) = 0.35;      % cross-event scale
sdv = [T.s_p(1) T.s_p(2) T.s_sd(1) T.s_sd(2)];
Sig = (sdv' * sdv) .* R;
U = randn(nsub, 4) * chol(Sig);
W = U(:, 3:4);                      % the per-participant scale effects

% correlated trial effects across events
Ri = [1 T.cor_i; T.cor_i 1];
sdi = [T.s_i(1) T.s_i(2)];
Sigi = (sdi' * sdi) .* Ri;
V = randn(ntrl, 2) * chol(Sigi);

ids = cell(nsub*ntrl*2, 1);
evt = cell(nsub*ntrl*2, 1);
meas = nan(nsub*ntrl*2, 1);
dim1 = nan(nsub*ntrl*2, 1);
row = 0;
for s = 1:nsub
    for e = 1:2
        logmu_p  = T.b(e)       + T.b_dim(e)       * z(s) + U(s, e);
        logsig_p = T.b_sigma(e) + T.b_sigma_dim(e) * z(s) + U(s, 2 + e);
        sig = exp(logsig_p);
        for t = 1:ntrl
            mu = exp(logmu_p + V(t, e));
            % scaled chi-square with mean mu and Var = sigma^2, i.e.
            % Gamma(shape = nu/2, rate = nu/(2*mu)) with nu = 2*(mu/sigma)^2.
            % An EXACT Gamma draw, not a Wilson-Hilferty chi-square: under
            % location-scale nu is derived and non-integer, so the small-nu tail
            % is genuinely visited (increment 7 rejected the approximation).
            nu = 2 * (mu / sig)^2;
            y = localGammaRand(nu/2) * (2*mu/nu);
            row = row + 1;
            ids{row}  = sprintf('s%03d', s);
            evt{row}  = sprintf('e%d', e);
            meas(row) = y;
            dim1(row) = z(s);
        end
    end
end
tbl = table(ids, meas, evt, dim1, ...
    'VariableNames', {'id','meas','event','dim1'});
end

% =========================================================================
function G = localTruePerSubjectG(T, z, W)
% Closed-form PRE-S17-FIX (conditional-only) per-participant difference-score G.
%
% NOT the production estimand any more: since S17-FIX the shipped residual
% pools the mean-surface person x trial term (see localTruePerSubjectGPooled,
% which is now the primary truth). This oracle is kept alive for exactly one
% assertion - the signed direction of what S17-FIX changed: pooling a
% nonnegative term can only lower each coefficient. Written from the section-J
% formulas directly rather than by calling the production converter.
n = numel(z);
G = nan(n,1);
for s = 1:n
    zz = z(s);
    aA = T.b(1) + T.b_dim(1)*zz;   aB = T.b(2) + T.b_dim(2)*zz;

    % population between-person components at this participant's z
    spA = T.s_p(1); spB = T.s_p(2); siA = T.s_i(1); siB = T.s_i(2);
    p2A = exp(2*aA + siA^2) * exp(spA^2) * (exp(spA^2) - 1);
    p2B = exp(2*aB + siB^2) * exp(spB^2) * (exp(spB^2) - 1);
    EYA = exp(aA + 0.5*(spA^2 + siA^2));
    EYB = exp(aB + 0.5*(spB^2 + siB^2));
    covP = EYA * EYB * (exp(T.cor_p * spA * spB) - 1);
    sigU = p2A + p2B - 2*covP;

    % this participant's own conditional residual per event
    sA = exp(2*(T.b_sigma(1) + T.b_sigma_dim(1)*zz + W(s,1)));
    sB = exp(2*(T.b_sigma(2) + T.b_sigma_dim(2)*zz + W(s,2)));
    relerr = sA/T.ntrl + sB/T.ntrl;

    G(s) = sigU / (sigU + relerr);
end
end

% =========================================================================
function G = localTruePerSubjectGPooled(T, z, W)
% Closed-form true PER-PARTICIPANT difference-score G under the POOLED estimand
% - the PRODUCTION estimand since S17-FIX, matching the owner's reference
% implementation: each event's residual is the population mean-surface
% person x trial term exp(2*alpha_k(z))*K_k plus the participant's own
% conditional variance. This is the primary truth the recovery is scored
% against. Written from the section-J formulas directly rather than by calling
% the production converter, so it is an independent check of the estimand and
% not a restatement of the implementation.
n = numel(z);
G = nan(n,1);
for s = 1:n
    zz = z(s);
    aA = T.b(1) + T.b_dim(1)*zz;   aB = T.b(2) + T.b_dim(2)*zz;
    spA = T.s_p(1); spB = T.s_p(2); siA = T.s_i(1); siB = T.s_i(2);
    p2A = exp(2*aA + siA^2) * exp(spA^2) * (exp(spA^2) - 1);
    p2B = exp(2*aB + siB^2) * exp(spB^2) * (exp(spB^2) - 1);
    EYA = exp(aA + 0.5*(spA^2 + siA^2));
    EYB = exp(aB + 0.5*(spB^2 + siB^2));
    covP = EYA * EYB * (exp(T.cor_p * spA * spB) - 1);
    sigU = p2A + p2B - 2*covP;

    sA = exp(2*(T.b_sigma(1) + T.b_sigma_dim(1)*zz + W(s,1)));
    sB = exp(2*(T.b_sigma(2) + T.b_sigma_dim(2)*zz + W(s,2)));

    % the POPULATION-AVERAGE mean-surface person x trial term per event,
    % exp(2*alpha_k)*K_k - the term the production pooling adds (S17-FIX).
    % This is the estimand, not an approximation: the per-person share V_s
    % would be a different, computable choice (J.6 derives it from u_s), and
    % neither the reference implementation nor PsyRAT uses it.
    KA = exp(spA^2 + siA^2)*(exp(spA^2) - 1)*(exp(siA^2) - 1);
    KB = exp(spB^2 + siB^2)*(exp(spB^2) - 1)*(exp(siB^2) - 1);
    vA = exp(2*aA)*KA;
    vB = exp(2*aB)*KB;

    relerr = (sA + vA)/T.ntrl + (sB + vB)/T.ntrl;
    G(s) = sigU / (sigU + relerr);
end
end

% =========================================================================
function M = localDrawMat(rel, field)
% Stratum 1 of a single-stratum (ungrouped, single-event-pair) design. cor_id
% comes back draws x 4 x 4, the rest draws x N; cell2mat preserves both.
M = cell2mat(rel.out.(field){1});
end

% =========================================================================
function localReportParam(testCase, draws, truth, label)
q = quantile(draws, [0.025 0.975]);
fprintf('  %-18s true %+.3f  est %+.3f [%+.3f, %+.3f]\n', ...
    label, truth, mean(draws), q(1), q(2));
testCase.verifyGreaterThanOrEqual(truth, q(1), ...
    sprintf('True %s below the credible interval.', label));
testCase.verifyLessThanOrEqual(truth, q(2), ...
    sprintf('True %s above the credible interval.', label));
end

% =========================================================================
function rel = localRunGammaDiffDynrelSserrLs(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~isfolder(outDir)
    mkdir(outDir);
end
rel = psyrat_computevarcomp('data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'chains', 4, ...
    'warmup', 1000, ...
    'sampling', 1000, ...
    'seed', 12345, ...
    'showgui', 0, ...
    'verbose', 1, ...
    'cmdstanoutdir', outDir);
end

% =========================================================================
function d = localCmdStanOutputDir(projectRoot)
d = fullfile(projectRoot, 'tests', 'cmdstan_out');
end

% =========================================================================
function r = localBandedSpearman(a, b, z, nband)
% Spearman(a,b) computed WITHIN quantile bands of z and averaged, so a covariate
% that drives most of a's marginal variation cannot dilute the relationship being
% tested. Bands rather than a linear partial correlation because the dependence on
% z is nonlinear here (everything scales as exp(2*alpha(z))).
if nargin < 4, nband = 5; end
edges = quantile(z(:), linspace(0,1,nband+1));
vals = nan(nband,1);
for k = 1:nband
    m = z(:) >= edges(k) & z(:) <= edges(k+1);
    if sum(m) >= 10
        vals(k) = localSpearman(a(m), b(m));
    end
end
r = mean(vals(~isnan(vals)));
end

function r = localSpearman(a, b)
% Spearman rank correlation in BASE MATLAB. corr() and tiedrank() are both
% Statistics Toolbox, and CI installs --products=MATLAB only, so they are not
% available; corrcoef is base. Ties are averaged, which matters here only in
% principle (posterior point estimates are continuous) but is cheap to do right.
r = corrcoef(localRank(a(:)), localRank(b(:)));
r = r(1,2);
end

function rk = localRank(x)
n = numel(x);
[~, ord] = sort(x);
rk = zeros(n,1);
rk(ord) = 1:n;
% average ranks within tied groups
[xs, ~] = sort(x);
k = 1;
while k <= n
    j = k;
    while j < n && xs(j+1) == xs(k)
        j = j + 1;
    end
    if j > k
        rk(ord(k:j)) = mean(k:j);
    end
    k = j + 1;
end
end

% =========================================================================
function x = localGammaRand(shape)
% Marsaglia-Tsang exact Gamma(shape,1) sampler, so the test does not depend on
% the Statistics Toolbox (CI installs base MATLAB only). Scale is applied by the
% caller. Boosted for shape < 1 via the standard u^(1/shape) trick.
if shape < 1
    x = localGammaRand(shape + 1) * rand()^(1/shape);
    return;
end
d = shape - 1/3;
c = 1/sqrt(9*d);
while true
    v = -1;
    while v <= 0
        zz = randn();
        v = (1 + c*zz)^3;
    end
    u = rand();
    if u < 1 - 0.0331*zz^4
        x = d*v; return;
    end
    if log(u) < 0.5*zz^2 + d*(1 - v + log(v))
        x = d*v; return;
    end
end
end
