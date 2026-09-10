classdef TestDiffRescorNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the two-event difference-
    % score model WITH residual covariance (case 7; diffest = 2, diffrescor = 2).
    % This is the bivariate (multi_normal_cholesky) sibling of the no-rescor
    % difference model (TestDiffNativeHmcRecovery): paired event responses with
    % constant per-event residual SDs exp(b_sigma), a residual correlation
    % (Lrescor), and the same two correlated 2-D person/trial random effects on
    % the two event means. It exercises the native HMC engine through the full
    % psyrat_computevarcomp pipeline, including the diff_rescor native_extra
    % priors threading and the rescor/err_varcov extraction.
    %
    % The CmdStan path for this exact design is validated separately by
    % TestRescorDiffRecovery (the B10 guard); this test confirms the native HMC
    % engine recovers the SAME known truth without CmdStan. The data generator is
    % the concurrent-pair simulator shared with that test (each (subject,trial)
    % yields paired event-A/event-B measurements with subject random intercepts
    % and bivariate residuals correlated at truth.rho).
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients). It needs
    % the Statistics and Machine Learning Toolbox (hmcSampler) and the Deep Learning
    % Toolbox (dlarray); when either is absent it skips cleanly (Incomplete). Slow
    % (it samples a ~90-parameter model), so it lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (two events, concurrent pairs) ----
            % Subject SD = 3 per event, residual SD = 2 per event, residual
            % correlation 0.4, no trial-position effect (the trial facet should
            % shrink toward 0, not absorb the residual). Mirrors the CmdStan
            % recovery regime in TestRescorDiffRecovery.
            truth = struct('rho', 0.4, 'sd_id', [3 3], 'sigma', [2 2], ...
                'sd_trl', [0 0], 'b', [5 3]);
            tbl = localSimulateConcurrentPairs(truth, 28, 16, 20260627);

            % ---- run native HMC through the full pipeline (diffest=2, diffwpcov=2) ----
            rel = psyrat_computevarcomp( ...
                'data', tbl, ...
                'sserrvar', 1, ...
                'diffest', 2, ...
                'diffwpcov', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws (rescor layout) ----
            b = cell2mat(rel.out.b{1});            % draws x 2 (per-event mean)
            bsig = cell2mat(rel.out.b_sigma{1});   % draws x 2 (per-event log-resid)
            sdid = cell2mat(rel.out.sd_id{1});     % draws x 2
            sdtrl = cell2mat(rel.out.sd_trl{1});   % draws x 2
            idvc = cell2mat(rel.out.id_varcov{1}); % draws x 2 x 2
            rescor = cell2mat(rel.out.rescor{1});  % draws x 1
            errvc = cell2mat(rel.out.err_varcov{1}); % draws x 2 x 2
            nd = size(b, 1);

            % the rescor branch must populate rescor + err_varcov (the no-rescor
            % branch does NOT, so this also guards the branch routing)
            testCase.verifyEqual(size(b), [nd 2]);
            testCase.verifyEqual(size(rescor), [nd 1]);
            testCase.verifyEqual(size(errvc), [nd 2 2]);
            testCase.verifyEqual(size(idvc), [nd 2 2]);

            sigmaHat = mean(exp(bsig), 1);

            % ---- recovery: posterior mean within a generous absolute band of the
            % population truth (same philosophy as the other native recovery tests
            % and TestRescorDiffRecovery's bounds). ----
            % residual correlation (the headline rescor parameter), true 0.4
            testCase.verifyGreaterThan(mean(rescor), 0.15, ...
                sprintf('Residual correlation under-recovered (rho_hat=%.3f, true=0.40).', mean(rescor)));
            testCase.verifyLessThan(mean(rescor), 0.65, ...
                sprintf('Residual correlation over-recovered (rho_hat=%.3f, true=0.40).', mean(rescor)));
            % per-event means, true [5 3]
            testCase.verifyEqual(mean(b(:, 1)), truth.b(1), 'AbsTol', 1.2);
            testCase.verifyEqual(mean(b(:, 2)), truth.b(2), 'AbsTol', 1.2);
            % residual SDs, true 2 per event
            testCase.verifyGreaterThan(min(sigmaHat), 1.4, ...
                sprintf('Residual SD under-recovered (sigma=[%.2f %.2f], true=2).', sigmaHat(1), sigmaHat(2)));
            testCase.verifyLessThan(max(sigmaHat), 2.6, ...
                sprintf('Residual SD over-recovered (sigma=[%.2f %.2f], true=2).', sigmaHat(1), sigmaHat(2)));
            % subject SDs, true 3 per event
            testCase.verifyGreaterThan(min(mean(sdid, 1)), 2.0, ...
                'Subject SD under-recovered (true 3).');
            testCase.verifyLessThan(max(mean(sdid, 1)), 4.0, ...
                'Subject SD over-recovered (true 3).');
            % trial facet should be small (true 0), i.e. not absorbing residual
            testCase.verifyLessThan(max(mean(sdtrl, 1)), 1.2, ...
                'Trial SD inflated (true 0); trial facet may be confounding with residual.');

            % ---- internal consistency of the rescor mapping ----
            % err_varcov diagonal must equal sigma^2 = exp(b_sigma)^2 per draw
            testCase.verifyEqual(mean(errvc(:, 1, 1)), mean(exp(bsig(:, 1)).^2), 'RelTol', 1e-6);
            testCase.verifyEqual(mean(errvc(:, 2, 2)), mean(exp(bsig(:, 2)).^2), 'RelTol', 1e-6);
            % err_varcov off-diagonal must equal rescor * sigma1 * sigma2 per draw
            offExpected = mean(rescor .* exp(bsig(:, 1)) .* exp(bsig(:, 2)));
            testCase.verifyEqual(mean(errvc(:, 1, 2)), offExpected, 'RelTol', 1e-6);
            testCase.verifyEqual(mean(errvc(:, 2, 1)), offExpected, 'RelTol', 1e-6);
            % id_varcov diagonal must equal sd_id^2 per draw
            testCase.verifyEqual(mean(idvc(:, 1, 1)), mean(sdid(:, 1).^2), 'RelTol', 1e-6);
            testCase.verifyEqual(mean(idvc(:, 2, 2)), mean(sdid(:, 2).^2), 'RelTol', 1e-6);
        end
    end
end

function tbl = localSimulateConcurrentPairs(truth, nSub, nTrial, seed)
% Deterministic concurrent two-event data: each (subject,trial) yields paired
% event-A/event-B measurements with subject random intercepts and bivariate
% residuals correlated at truth.rho. No trial-position effect (truth.sd_trl=0).
% Identical generator to TestRescorDiffRecovery so the native HMC engine is
% checked against the same known truth as the CmdStan path.
rng(seed);
sigma = truth.sigma;
covE = [sigma(1)^2, truth.rho*sigma(1)*sigma(2); ...
        truth.rho*sigma(1)*sigma(2), sigma(2)^2];
cholE = chol(covE, 'lower');

ids = cell(nSub*nTrial*2, 1);
evs = cell(nSub*nTrial*2, 1);
meas = zeros(nSub*nTrial*2, 1);
k = 0;
for s = 1:nSub
    a = truth.sd_id(:) .* randn(2, 1);
    for t = 1:nTrial
        e = cholE * randn(2, 1);
        k = k + 1;
        ids{k} = sprintf('s%02d', s); evs{k} = 'A';
        meas(k) = truth.b(1) + a(1) + e(1);
        k = k + 1;
        ids{k} = sprintf('s%02d', s); evs{k} = 'B';
        meas(k) = truth.b(2) + a(2) + e(2);
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});
end
