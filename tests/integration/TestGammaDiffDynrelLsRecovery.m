classdef TestGammaDiffDynrelLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE Gamma /
    % scaled-chi-square one-facet GROUP-LEVEL DYNAMIC DIFFERENCE design
    % (analysis 12, family='gamma', gammascale=2). Increment 8b of the
    % gammascale rollout; the first location-scale design that has cross-event
    % covariance AND dimension slopes on both submodels.
    %
    % Sibling of TestGammaDiffDynrelRecovery, which covers the log-nu twin.
    % Two-event non-concurrent data are simulated from a log-linked per-event
    % mean and a per-event, person-varying log RESIDUAL SD, both conditioned on
    % a standardized subject-level dimension z:
    %
    %   log mu_ep  = b_e       + b_dim_e       * z_p + u_mu,ep + u_tau,ei
    %   log sig_ep = b_sigma_e + b_sigma_dim_e * z_p + w_ep
    %   nu         = 2*(mu/sigma)^2,   Var(Y | mu, sigma) = sigma^2 exactly
    %
    % with the person block [u_mu,1 u_mu,2 w_1 w_2] correlated (4-D) and the
    % trial block [u_tau,1 u_tau,2] correlated across events.
    %
    % WHAT IS ASSERTED. The per-event mean slopes b_dim, the per-event
    % residual-SD slopes b_sigma_dim, and - decisively - the z-conditioned
    % DIFFERENCE coefficient surface must cover the closed-form true surface.
    %
    % b_sigma_dim IS ASSERTED, following increment 8a's reasoning: it is a
    % directly SAMPLED parameter of the fitted model, not a derived contrast of
    % separately estimated per-cell residuals. Increments 6 and 7 left their
    % static residual-SD CONTRAST as a printed diagnostic because that quantity
    % compounds per-cell estimation errors. The distinguishing variable is
    % sampled-parameter vs derived-difference, not the design family.
    %
    % Base MATLAB only (Marsaglia-Tsang gamma draw, no Statistics Toolbox), so it
    % survives CI's --products=MATLAB image. Under location-scale
    % nu = 2*(mu/sigma)^2 is DERIVED and non-integer, so the small-nu tail is
    % genuinely visited and a Wilson-Hilferty or sum-of-squared-normals
    % approximation would misstate the data-generating process.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).
    %
    % ---- THE DESIGN POINT, AND THE FOUR TRAPS IT AVOIDS ----
    % Planted (see the offline sweep recorded in gamma_locationscale_rollout.md,
    % increment 8b):
    %
    %   b           = [1.75  1.60]     per-event log-mean intercepts
    %   b_sigma     = [0.90  1.05]     per-event log residual-SD intercepts
    %   b_dim       = [0.30  0.05]     mean slopes, ASYMMETRIC ACROSS events
    %   b_sigma_dim = [0.20 -0.10]     scale slopes, ASYMMETRIC WITHIN the pair
    %   n' per event = [20 20]
    %
    % Two events give MORE cancellation routes than increment 8a's single-event
    % design, so each was checked offline before any fit:
    %
    %  (1) G ON THE CEILING. The first sweep sat at G ~ 0.99 at every z, where the
    %      coefficient is insensitive to everything and a recovery test asserts
    %      nothing. Raising the residual intercepts to ~0.9 brings the surface to
    %      G = 0.688 -> 0.910 across z in [-2, +2] (range 0.222).
    %  (2) SYMMETRIC SCALE SLOPES. With b_sigma_dim equal across events the two
    %      residual SDs move together, so the CONTRAST - the quantity this
    %      parameterization exists to separate - is nearly blind: the sweep gave a
    %      contrast range of 0.093 against 1.07 for the asymmetric variant, while G
    %      still looked healthy. The planted slopes have OPPOSITE SIGNS, and the
    %      resulting per-event SDs move in opposite directions (+123% vs -33%),
    %      with the contrast crossing zero (-1.84 -> +1.33).
    %  (3) NO MEAN ASYMMETRY. Ruling 15 says asymmetric event mean slopes
    %      legitimately move the coefficient. b_dim = [0.30 0.05] plants that
    %      asymmetry so the test would catch an implementation that suppressed it.
    %  (4) A DEGENERATE LOG-NU TWIN. The implied log-nu slopes are
    %      b_nu_dim = 2*(b_dim - b_sigma_dim) = [0.200  0.300]. Both are well away
    %      from zero, so the paired cross-parameterization fit has a real
    %      dispersion slope to recover on BOTH events. Increment 8a's trap was a
    %      point where that implied slope vanished.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            base = localCmdStanOutputDir(testCase.projectRoot());
            for tag = {'recov_gamma_diffdynrel_ls'}
                d = fullfile(base, tag{1});
                if exist(d, 'dir'), rmdir(d, 's'); end
            end
        end
    end

    methods (Test)

        function testDifferenceSurfaceAndSlopesRecovered(testCase)
            T = localTruth();
            rng(20260805, 'twister');   % data generation only; estimation seed is fixed at 12345
            tbl = localSimulateGammaDiffDynrelLs(T);

            rel = localRunDiffDynrelGammaLs(testCase, tbl, 'recov_gamma_diffdynrel_ls');

            % Draws straight out of REL.out and into the calculator, matching the
            % idiom of TestGammaDiffDynrelRecovery (the log-nu twin). Going
            % through psyrat_dynrel_summary would work too but adds the stratum
            % indirection for no gain in a single-stratum design.
            bb    = localDrawMat(rel, 'b');
            bsig  = localDrawMat(rel, 'b_sigma');
            bdim  = localDrawMat(rel, 'b_dim');
            bsdim = localDrawMat(rel, 'b_sigma_dim');
            sdid  = localDrawMat(rel, 'sd_id');
            sdtrl = localDrawMat(rel, 'sd_trl');
            corid = localDrawMat(rel, 'cor_id');
            cori  = localDrawMat(rel, 'cor_i');

            % ---- the coefficient surface against closed-form truth ----
            zq = [-1 0 1];
            trueG = localTrueG(T, zq);
            rs = psyrat_rel_diffdynrel_gamma_ls('b',bb,'b_sigma',bsig, ...
                'b_dim',bdim,'b_sigma_dim',bsdim,'sd_id',sdid,'sd_trl',sdtrl, ...
                'cor_id',corid,'cor_i',cori,'ndim',1,'z1',zq, ...
                'obs',[T.ntrl T.ntrl],'CI',0.95);

            fprintf('\n--- case-12 LS recovery: difference G surface ---\n');
            for k = 1:numel(zq)
                fprintf('  z=%+0.1f  true %.4f   est %.4f [%.4f, %.4f]\n', ...
                    zq(k), trueG(k), rs.G.pt(k), rs.G.ll(k), rs.G.ul(k));
                testCase.verifyGreaterThanOrEqual(trueG(k), rs.G.ll(k), ...
                    sprintf('True G(%+0.1f) below the credible interval.', zq(k)));
                testCase.verifyLessThanOrEqual(trueG(k), rs.G.ul(k), ...
                    sprintf('True G(%+0.1f) above the credible interval.', zq(k)));
            end

            % The surface must RISE across z: event 1's residual SD grows with z
            % but event 2's shrinks faster in the contrast, and the mean slopes
            % are both positive. Pins the direction, not just the coverage.
            testCase.verifyGreaterThan(rs.G.pt(end), rs.G.pt(1), ...
                'Recovered G(z) should rise across z at this design point.');

            % ---- the slopes, both submodels, per event ----
            fprintf('--- slopes (odd cols event 1, even cols event 2) ---\n');
            localReportSlope(testCase, bdim(:,1), T.b_dim(1),       'b_dim event 1');
            localReportSlope(testCase, bdim(:,2), T.b_dim(2),       'b_dim event 2');
            localReportSlope(testCase, bsdim(:,1), T.b_sigma_dim(1), 'b_sigma_dim event 1');
            localReportSlope(testCase, bsdim(:,2), T.b_sigma_dim(2), 'b_sigma_dim event 2');

            % The ASYMMETRY of the scale slopes is the interpretive payload: the
            % two events' residual SDs must be recovered as moving in OPPOSITE
            % directions, which is what a log-nu fit cannot report.
            contrast = bsdim(:,1) - bsdim(:,2);
            trueContrast = T.b_sigma_dim(1) - T.b_sigma_dim(2);
            q = quantile(contrast, [0.025 0.975]);
            fprintf('  b_sigma_dim CONTRAST: true %+.3f  est %+.3f [%+.3f, %+.3f]\n', ...
                trueContrast, mean(contrast), q(1), q(2));
            testCase.verifyGreaterThanOrEqual(trueContrast, q(1), ...
                'True scale-slope contrast below the credible interval.');
            testCase.verifyLessThanOrEqual(trueContrast, q(2), ...
                'True scale-slope contrast above the credible interval.');
            testCase.verifyGreaterThan(mean(contrast > 0), 0.95, ...
                ['The recovered scale slopes must be ordered event1 > event2 - ' ...
                'that opposite-direction movement is what this parameterization ' ...
                'exists to expose.']);
        end

    end
end

% =========================================================================
function T = localTruth()
% The audited design point. Changing any of these invalidates the four-trap
% argument in the class header; re-run the offline sweep if you do.
T = struct();
T.b           = [1.75 1.60];
T.b_sigma     = [0.90 1.05];
T.b_dim       = [0.30 0.05];
T.b_sigma_dim = [0.20 -0.10];
T.s_p         = [0.32 0.29];   % person log-MEAN SDs, per event
T.s_sd        = [0.24 0.21];   % person log-SIGMA SDs, per event
T.s_i         = [0.11 0.13];   % trial log-mean SDs, per event
T.cor_p       = 0.45;          % cross-event person-mean correlation
T.cor_i       = 0.20;          % cross-event trial correlation

% SIZED FOR THE PERSON BLOCK, NOT THE SLOPES. The first run of this test used
% nsub = 160, ntrl = 20 and FAILED the surface assertion while recovering all
% four dimension slopes essentially perfectly. The diagnosis (recorded in
% gamma_locationscale_rollout.md, increment 8b) was that the arithmetic is right
% and the PERSON BLOCK was the weak link:
%   - cor_p recovered 0.40 vs 0.45. LKJ(1) on a 4-D block puts a marginal prior
%     on each correlation that shrinks toward 0, and 160 people did not overcome
%     it. A smaller cor_p inflates sigma^2_p,Delta = A + B - 2*cov, raising G.
%   - s_p(2) recovered 0.33 vs 0.29. Event 2 carries the larger residual
%     (sigma ~ 2.9) against only 20 trials, so its person mean is estimated
%     loosely and the person SD absorbs the slack.
% Together those two moved G(-1) by +0.049 of a +0.072 total gap, verified by
% perturbing each in isolation through the closed form.
%
% More trials per event attacks the second directly and more people the first, so
% both were raised. This is the principled fix; loosening the assertion would
% have hidden a real identification property behind a green test.
T.nsub        = 250;
T.ntrl        = 40;            % per event; n' = [40 40]
end

% =========================================================================
function tbl = localSimulateGammaDiffDynrelLs(T)
% Two-event non-concurrent data from the location-scale generative model. The
% events are NON-CONCURRENT, so each event gets its OWN trial index and the
% conditional draws are independent given (mu, sigma) - that is the assumption
% the estimand rests on, so it is built into the simulation rather than assumed.
nsub = T.nsub; ntrl = T.ntrl;

% subject-level dimension, standardized across the whole sample
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

% correlated trial effects across events, one set of levels per event index
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
            % Gamma(shape = nu/2, rate = nu/(2*mu)) with nu = 2*(mu/sigma)^2
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
function G = localTrueG(T, zq)
% Closed-form true difference-score G at each z, from the SAME converters the
% production path uses. This is a regression-style truth in the sense that it
% reuses psyrat_gamma_varcomps_ls / psyrat_gamma_crosscov, but the PARAMETERS are
% the planted ones rather than anything estimated, so it does check that the
% fitted model recovers the generative truth.
G = nan(size(zq));
for k = 1:numel(zq)
    zz = zq(k);
    aA = T.b(1) + T.b_dim(1)*zz;   aB = T.b(2) + T.b_dim(2)*zz;
    lA = T.b_sigma(1) + T.b_sigma_dim(1)*zz;
    lB = T.b_sigma(2) + T.b_sigma_dim(2)*zz;
    vA = psyrat_gamma_varcomps_ls(aA, T.s_p(1), T.s_i(1), lA, T.s_sd(1));
    vB = psyrat_gamma_varcomps_ls(aB, T.s_p(2), T.s_i(2), lB, T.s_sd(2));
    cc = psyrat_gamma_crosscov(aA, T.s_p(1), T.s_i(1), ...
                               aB, T.s_p(2), T.s_i(2), T.cor_p, T.cor_i);
    vp  = vA.sigma_p2 + vB.sigma_p2 - 2*cc.cov_p_obs;
    err = vA.sigma_pi_e2/T.ntrl + vB.sigma_pi_e2/T.ntrl;
    G(k) = vp / (vp + err);
end
end

% =========================================================================
function M = localDrawMat(rel, field)
% Stratum 1 of a single-stratum (ungrouped, single-event-pair) design. cor_id
% comes back draws x 4 x 4, the rest draws x N; cell2mat preserves both.
M = cell2mat(rel.out.(field){1});
end

% =========================================================================
function localReportSlope(testCase, draws, truth, label)
q = quantile(draws, [0.025 0.975]);
fprintf('  %-22s true %+.3f  est %+.3f [%+.3f, %+.3f]\n', ...
    label, truth, mean(draws), q(1), q(2));
testCase.verifyGreaterThanOrEqual(truth, q(1), ...
    sprintf('True %s below the credible interval.', label));
testCase.verifyLessThanOrEqual(truth, q(2), ...
    sprintf('True %s above the credible interval.', label));
end

% =========================================================================
function rel = localRunDiffDynrelGammaLs(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

% =========================================================================
function d = localCmdStanOutputDir(projectRoot)
d = fullfile(projectRoot, 'tests', 'cmdstan_out');
if ~exist(d, 'dir'), mkdir(d); end
end

% =========================================================================
function x = localGammaRand(shape)
% Exact Gamma(shape, scale = 1) draw by the Marsaglia-Tsang (2000) method, using
% only base-MATLAB randn/rand so the test runs on a CI image installed with
% --products=MATLAB. Statistics Toolbox gamrnd would be simpler and is not
% available there.
%
% Under location-scale nu = 2*(mu/sigma)^2 is DERIVED and non-integer, and this
% design visits small nu in the tail, so an approximation here would misstate the
% data-generating process in a known-truth test.
if shape < 1
    x = localGammaRand(shape + 1) * rand^(1/shape);
    return;
end
d = shape - 1/3;
c = 1 / sqrt(9*d);
while true
    v = -1;
    while v <= 0
        z = randn;
        v = (1 + c*z)^3;
    end
    u = rand;
    if log(u) < 0.5*z^2 + d - d*v + d*log(v)
        x = d * v;
        return;
    end
end
end
