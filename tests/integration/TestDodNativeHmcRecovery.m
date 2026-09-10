classdef TestDodNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the task-free four-event
    % difference-of-differences model (analysis 9; diffest = 3). This is a 4-cell
    % location-scale model with person AND trial random effects on both the
    % location (cells 1-4) and the scale (cells 1-4) - an 8-D joint block per
    % grouping factor via an 8x8 LKJ Cholesky (psyrat_corr_chol_lkj). Matches
    % psyrat_build_stan_dodiff.
    %
    % There is no CmdStan-path recovery test for this design to mirror, so the
    % generator is written here: a balanced person x cell x trial design with
    % INDEPENDENT person-cell and trial-cell location effects (so the location
    % var-covariances are diagonal), a per-cell residual, and small per-cell
    % scale random effects. The model must recover the person-cell location
    % variances (id_varcov diagonal), the population per-cell residual variances
    % (err_var_Cell = exp(2*b_sigma_cell)), and -- the headline DoD quantity --
    % the difference-of-differences universe-score variance c' Sigma_id c with the
    % contrast c = [1 -1 -1 1] (which, for independent cells, equals the sum of the
    % four person-cell location variances).
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. The 8-D person/trial
    % blocks (two 8x8 LKJ Choleskys) make this a high-dimensional model, so it is
    % slow and lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (four cells) ----
            rng(20260628, 'twister');
            nsub = 30; ntrl = 8; events = {'e1','e2','e3','e4'};
            cellMu = [5 3 4 2];
            cellLogResid = [log(2) log(2) log(1.8) log(1.8)];
            sdLoc = [3 2.5 2 2.5];   % person-cell LOCATION SDs (independent cells)
            sdSc  = [0.3 0.3 0.3 0.3]; % person-cell scale SDs (small)
            sdTl  = [1 1 0.9 0.9];   % trial-cell LOCATION SDs
            sdTs  = [0.2 0.2 0.2 0.2]; % trial-cell scale SDs

            ids = {}; meas = []; ev = {};
            for p = 1:nsub
                aP = [sdLoc sdSc]' .* randn(8, 1);   % person 8-D effects (4 loc, 4 scale)
                aT = [sdTl sdTs]' .* randn(8, ntrl); % trial 8-D effects (by position)
                for e = 1:4
                    for t = 1:ntrl
                        mu = cellMu(e) + aP(e) + aT(e, t);
                        lsg = cellLogResid(e) + aP(4+e) + aT(4+e, t);
                        ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                        ev{end+1,1} = events{e}; %#ok<AGROW>
                        meas(end+1,1) = mu + exp(lsg) * randn; %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, 'VariableNames', {'id', 'meas', 'event'});

            % ---- run native HMC through the full pipeline (diffest = 3) ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'diffest', 3, ...
                'dodmap', events, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_dodiff');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            bsig = rel.out.b_sigma{1};      % draws x 4 (b_sigma_cell)
            idvc = rel.out.id_varcov{1};    % draws x 4 x 4 (person-cell location var-cov)
            ervar = rel.out.err_var{1};     % draws x 4 (exp(2*b_sigma_cell))
            corid = rel.out.cor_id{1};      % draws x 4 x 4
            nd = size(bsig, 1);
            testCase.verifyEqual(size(bsig), [nd 4]);
            testCase.verifyEqual(size(idvc), [nd 4 4]);

            % cor_id is a correlation block (unit diagonal) on a sample draw
            testCase.verifyEqual([corid(1,1,1) corid(1,2,2) corid(1,3,3) corid(1,4,4)], ...
                [1 1 1 1], 'AbsTol', 1e-6);

            % ---- recovery: person-cell LOCATION variances (id_varcov diagonal) ----
            for c = 1:4
                testCase.verifyEqual(mean(idvc(:, c, c)), sdLoc(c)^2, ...
                    'AbsTol', max(3.0, 0.5*sdLoc(c)^2), ...
                    sprintf('id_varcov(%d,%d) off (true %.2f).', c, c, sdLoc(c)^2));
            end

            % ---- recovery: population per-cell residual variances ----
            for c = 1:4
                testCase.verifyEqual(mean(ervar(:, c)), exp(2*cellLogResid(c)), ...
                    'AbsTol', max(1.5, 0.4*exp(2*cellLogResid(c))), ...
                    sprintf('err_var_Cell(%d) off (true %.2f).', c, exp(2*cellLogResid(c))));
            end

            % ---- headline DoD universe-score variance: c' Sigma_id c, c=[1 -1 -1 1].
            % For independent cells Sigma_id is diagonal, so the truth is the sum of
            % the four person-cell location variances. ----
            cvec = [1 -1 -1 1].';
            uniDraws = zeros(nd, 1);
            for d = 1:nd
                Sig = squeeze(idvc(d, :, :));
                uniDraws(d) = cvec.' * Sig * cvec;
            end
            uniTrue = sum(sdLoc.^2);
            testCase.verifyEqual(mean(uniDraws), uniTrue, 'AbsTol', max(6.0, 0.4*uniTrue));
            testCase.verifyGreaterThan(mean(uniDraws), 0);
        end
    end
end
