classdef TestDiffNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the two-event difference-
    % score model without residual covariance (case 7; diffest = 2, diffrescor = 1).
    % This is a univariate Gaussian over the stacked event data with an event-
    % indicator design (per-event mean b and log-residual b_sigma) and two
    % correlated 2-D random effects (person, trial), each a non-centered 2x2
    % Cholesky. It exercises the native HMC engine on a multivariate-style design
    % (distinct from the location-scale dynrel/sserr family) through the full
    % psyrat_computevarcomp pipeline, including the native_extra priors threading.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients). It needs
    % the Statistics and Machine Learning Toolbox (hmcSampler) and the Deep Learning
    % Toolbox (dlarray); when either is absent it skips cleanly (Incomplete). Slow
    % (it samples a larger model), so it lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (two events) ----
            rng(20260626, 'twister');
            nsub = 40; ntrl = 12;
            b1 = 5.0; b2 = 3.0; bs1 = log(1.3); bs2 = log(1.0);
            sdId = [2.0 1.4]; rhoId = 0.4;
            sdTrl = [1.2 0.9]; rhoTrl = 0.3;

            Sid = [sdId(1)^2, rhoId*sdId(1)*sdId(2); rhoId*sdId(1)*sdId(2), sdId(2)^2];
            Lid = chol(Sid, 'lower'); r1 = Lid * randn(2, nsub);
            Strl = [sdTrl(1)^2, rhoTrl*sdTrl(1)*sdTrl(2); rhoTrl*sdTrl(1)*sdTrl(2), sdTrl(2)^2];
            Ltrl = chol(Strl, 'lower'); r2 = Ltrl * randn(2, ntrl);

            ids = {}; meas = []; ev = {};
            for s = 1:nsub
                for t = 1:ntrl
                    ids{end+1,1} = sprintf('S%02d', s); %#ok<AGROW>
                    ev{end+1,1} = 'E1'; %#ok<AGROW>
                    meas(end+1,1) = b1 + r1(1,s) + r2(1,t) + exp(bs1)*randn; %#ok<AGROW>
                    ids{end+1,1} = sprintf('S%02d', s); %#ok<AGROW>
                    ev{end+1,1} = 'E2'; %#ok<AGROW>
                    meas(end+1,1) = b2 + r1(2,s) + r2(2,t) + exp(bs2)*randn; %#ok<AGROW>
                end
            end
            dataTbl = table(ids, meas, ev, 'VariableNames', {'id', 'meas', 'event'});

            % ---- run native HMC through the full pipeline (diffest = 2) ----
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'diffest', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws (no-rescor layout) ----
            b = cell2mat(rel.out.b{1});            % draws x 2 (per-event mean)
            bsig = cell2mat(rel.out.b_sigma{1});   % draws x 2 (per-event log-resid)
            sdid = cell2mat(rel.out.sd_id{1});     % draws x 2
            sdtrl = cell2mat(rel.out.sd_trl{1});   % draws x 2
            idvc = cell2mat(rel.out.id_varcov{1}); % draws x 2 x 2
            nd = size(b, 1);
            testCase.verifyEqual(size(b), [nd 2]);
            testCase.verifyEqual(size(bsig), [nd 2]);
            testCase.verifyEqual(size(idvc), [nd 2 2]);

            % ---- recovery: posterior mean within a generous absolute band of
            % the population truth (same philosophy as the other recovery tests).
            % Reference recovery at this fixed config/seed: b_sigma ~ [0.26 0.01],
            % sd_id ~ [2.0 1.4], sd_trl ~ [1.5 1.2]; the per-event means and the
            % weakly-identified trial SDs/correlations use wider bands.
            testCase.verifyEqual(mean(b(:, 1)), b1, 'AbsTol', 1.2);
            testCase.verifyEqual(mean(b(:, 2)), b2, 'AbsTol', 1.0);
            testCase.verifyEqual(mean(bsig(:, 1)), bs1, 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsig(:, 2)), bs2, 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:, 1)), sdId(1), 'AbsTol', 0.7);
            testCase.verifyEqual(mean(sdid(:, 2)), sdId(2), 'AbsTol', 0.6);
            testCase.verifyEqual(mean(sdtrl(:, 1)), sdTrl(1), 'AbsTol', 0.9);
            testCase.verifyEqual(mean(sdtrl(:, 2)), sdTrl(2), 'AbsTol', 0.9);

            % id_varcov diagonal must equal sd_id^2 (consistency of the mapping)
            testCase.verifyEqual(mean(idvc(:, 1, 1)), mean(sdid(:, 1).^2), 'RelTol', 1e-6);
            testCase.verifyEqual(mean(idvc(:, 2, 2)), mean(sdid(:, 2).^2), 'RelTol', 1e-6);
        end
    end
end
