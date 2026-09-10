classdef TestSserrNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the subject-level error-
    % variance model (case 6; sserrvar = 2). This is the base one-facet location-
    % scale model (person location + per-subject log-residual via a 2-D Cholesky,
    % plus a crossed trial main effect) with NO dimension predictors (KDIM = 0) -
    % the same native engine that fits dynrel (case 11), exercised on its KDIM = 0
    % path through the full psyrat_computevarcomp pipeline.
    %
    % Requires NO CmdStan: the native HMC engine is pure MATLAB (hmcSampler + AD
    % gradients via dlarray). It needs the Statistics and Machine Learning Toolbox
    % (hmcSampler) and the Deep Learning Toolbox (dlarray); when either is absent
    % the test skips cleanly (Incomplete). It is slow and lives in tests/integration
    % so the unit lane excludes it.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (no dimension predictors) ----
            rng(20260626, 'twister');
            nsub = 40; ntrl = 12;
            grand = 4.0; sigLog0 = log(1.2);
            sigP = 1.5; sigDP = 0.40; rho = 0.5; sigI = 0.7;

            S = [sigP^2, rho*sigP*sigDP; rho*sigP*sigDP, sigDP^2];
            L = chol(S, 'lower');
            e = L * randn(2, nsub);
            nu_p = e(1, :)';
            del_p = e(2, :)';
            beta_t = sigI * randn(ntrl, 1);

            ids = cell(nsub * ntrl, 1);
            meas = zeros(nsub * ntrl, 1);
            r = 0;
            for pp = 1:nsub
                for tt = 1:ntrl
                    r = r + 1;
                    mu = grand + nu_p(pp) + beta_t(tt);
                    sg = exp(sigLog0 + del_p(pp));
                    ids{r} = sprintf('S%02d', pp);
                    meas(r) = mu + sg * randn;
                end
            end
            dataTbl = table(ids, meas, 'VariableNames', {'id', 'meas'});

            % ---- run native HMC through the full pipeline (sserrvar = 2) ----
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'sserrvar', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_sserrvar');
            testCase.verifyEqual(rel.engine, 'hmc');

            gro = cell2mat(rel.out.gro_sds{1});      % draws x 2
            log0 = rel.out.pop_sdlog(:, 1);          % Intercept_sigma
            sigi = rel.out.sig_trl(:, 1);
            nd = size(gro, 1);
            testCase.verifyEqual(size(gro), [nd 2]);
            testCase.verifyEqual(size(cell2mat(rel.out.ind_bs{1})), [nd nsub]);
            testCase.verifyEqual(size(cell2mat(rel.out.ind_sdlog{1})), [nd nsub]);

            % ---- recovery: posterior mean within a generous absolute band of
            % the population truth (same philosophy as the CmdStan recovery tests;
            % the bands absorb the small-N realized-vs-population sampling gap and
            % MCMC variation while catching a gross engine regression).
            testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 0.6);   % sigma_p
            testCase.verifyEqual(mean(gro(:, 2)), sigDP, 'AbsTol', 0.3);  % sigma_dp
            testCase.verifyEqual(mean(log0), sigLog0, 'AbsTol', 0.25);    % Intercept_sigma
            testCase.verifyEqual(mean(sigi), sigI, 'AbsTol', 0.7);        % sigma_i (trial)
        end
    end
end
