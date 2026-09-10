classdef TestDiffDynrelNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the group-level dynamic
    % difference reliability model (Rast & Clayson, analysis 12; dynrel=2 +
    % diffest=2). This is the no-rescor two-event difference model EXTENDED with
    % event-specific fixed dimension slopes on BOTH the mean and the per-event
    % log-residual (exactly as dynrel extends sserr), matching
    % psyrat_build_stan_diffdynrel.
    %
    % The CmdStan path for this exact design is validated separately by
    % TestDiffDynrelRecovery; this test confirms the native HMC engine recovers
    % the SAME known truth without CmdStan, through the full psyrat_computevarcomp
    % pipeline (including the native_extra dual-prior threading: diff for the base
    % block, dynrel for the slope scales). The data generator is the shared
    % gain/loss simulator with correlated person/trial event effects and a
    % per-event residual whose log-SD depends on a subject-level standardized
    % dimension.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients). It needs
    % the Statistics and Machine Learning Toolbox (hmcSampler) and the Deep Learning
    % Toolbox (dlarray); when either is absent it skips cleanly (Incomplete). Slow
    % (samples a model with two crossed 2-D random effects), so it lives in
    % tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (gain/loss, one standardized dimension) ----
            rng(424242, 'twister');
            nsub = 50; ntrl = 20;
            Sp = [16 6; 6 9];                   % person event var-cov (gain,loss)
            mu = [5 2]; bmean = [1.0 0.5];      % cell means + mean dim slopes
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10]; % resid log-SD intercept + slope
            st = [1.5 1.2];                     % trial SDs (independent)
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % sigma^2_p_delta = 13

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            beta_t = [randn(ntrl,1)*st(1), randn(ntrl,1)*st(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain', 'loss'};
            for p = 1:nsub
                for e = 1:2
                    for t = 1:ntrl
                        mn = mu(e) + bmean(e)*zsub(p) + nu_p(p,e) + beta_t(t,e);
                        sg = exp(a(e) + bsig(e)*zsub(p));
                        ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                        ev{end+1,1} = events{e}; %#ok<AGROW>
                        meas(end+1,1) = mn + randn*sg; %#ok<AGROW>
                        dim1(end+1,1) = zsub(p); %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, dim1, ...
                'VariableNames', {'id', 'meas', 'event', 'dim1'});

            % ---- run native HMC through the full pipeline (dynrel=2, diffest=2) ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'dynrel', 2, ...
                'diffest', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});       % draws x 2 x 2
            trlvc = cell2mat(rel.out.trl_varcov{1});     % draws x 2 x 2
            bsg = cell2mat(rel.out.b_sigma{1});          % draws x 2 (resid intercepts)
            bsd = cell2mat(rel.out.b_sigma_dim{1});      % draws x KDIM (scale slopes)
            bmd = cell2mat(rel.out.b_dim{1});            % draws x KDIM (mean slopes)
            nd = size(bsg, 1);
            testCase.verifyEqual(size(bsg), [nd 2]);
            testCase.verifyEqual(size(bsd), [nd 2]);   % KDIM = 2 for one dimension
            testCase.verifyEqual(size(bmd), [nd 2]);
            testCase.verifyEqual(size(idvc), [nd 2 2]);

            % ---- recovery: generous absolute bands on the population truth
            % (mirrors TestDiffDynrelRecovery's bounds). ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 6.0);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.20);

            % ---- dynamic difference reliability G_delta(z) vs the closed form ----
            zg = [-1.5 0 1.5]; nprime = 15;
            out = psyrat_rel_diffdynrel('id_varcov', idvc, ...
                'trl_varcov', trlvc, 'b_sigma', bsg, ...
                'b_sigma_dim', bsd, 'ndim', 1, 'z1', zg, 'obs', nprime, 'CI', .95);
            for k = 1:numel(zg)
                wp1 = exp(a(1) + bsig(1)*zg(k))^2;
                wp2 = exp(a(2) + bsig(2)*zg(k))^2;
                Gtr = uniTrue / (uniTrue + (wp1 + wp2)/nprime);
                testCase.verifyEqual(out.G.pt(k), Gtr, 'AbsTol', 0.12);
            end
        end
    end
end
