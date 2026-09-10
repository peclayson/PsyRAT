classdef TestSubjDiffDynrelNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the SUBJECT-LEVEL dynamic
    % difference reliability model (Rast & Clayson, analysis 13; dynrel=2 +
    % diffest=2 + sserrvar=2). This is the group-level dynamic difference
    % (TestDiffDynrelNativeHmcRecovery) with the person random-effect block widened
    % to a 4-D joint correlated block over {location event1, location event2, scale
    % event1, scale event2}, so each participant has their own per-event residual
    % via a 4x4 LKJ Cholesky (psyrat_corr_chol_lkj). Matches
    % psyrat_build_stan_diffdynrel_sserr.
    %
    % The CmdStan path for this exact design is validated separately by
    % TestSubjDiffDynrelRecovery; this test confirms the native HMC engine recovers
    % the SAME known truth without CmdStan, through the full psyrat_computevarcomp
    % pipeline, including the 4x4 LKJ person block and the per-participant
    % phi_delta,s(z) table (psyrat_dynrel_summary). The data generator is the
    % shared gain/loss simulator with correlated person/trial LOCATION effects and
    % a per-person event-specific residual whose log-SD depends on a subject-level
    % standardized dimension AND a person-specific scale random effect.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients). It needs
    % the Statistics and Machine Learning Toolbox (hmcSampler) and the Deep Learning
    % Toolbox (dlarray); when either is absent it skips cleanly (Incomplete). It is
    % the largest native HMC model (a 4-D person block + a 2-D trial block), so it
    % is slow and lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (gain/loss, one standardized dimension, per-person
            % residual-scale random effect) ----
            rng(515151, 'twister');
            nsub = 45; ntrl = 20;
            Sp = [16 6; 6 9];                   % person LOCATION event var-cov
            mu = [5 2]; bmean = [1.0 0.5];      % cell means + mean dim slopes
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10]; % resid log-SD intercept + slope
            sds = [0.5 0.4];                    % between-person residual-scale SDs
            st = [1.5 1.2];                     % trial SDs (independent)
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % sigma^2_p_delta = 13

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)]; % person scale RE
            beta_t = [randn(ntrl,1)*st(1), randn(ntrl,1)*st(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain', 'loss'};
            for p = 1:nsub
                for e = 1:2
                    for t = 1:ntrl
                        mn = mu(e) + bmean(e)*zsub(p) + nu_p(p,e) + beta_t(t,e);
                        sg = exp(a(e) + bsig(e)*zsub(p) + delta(p,e));
                        ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                        ev{end+1,1} = events{e}; %#ok<AGROW>
                        meas(end+1,1) = mn + randn*sg; %#ok<AGROW>
                        dim1(end+1,1) = zsub(p); %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, dim1, ...
                'VariableNames', {'id', 'meas', 'event', 'dim1'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'dynrel', 2, ...
                'diffest', 2, ...
                'sserrvar', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});   % draws x 2 x 2 (location block)
            bsg = cell2mat(rel.out.b_sigma{1});      % draws x 2 (resid intercepts)
            bsd = cell2mat(rel.out.b_sigma_dim{1});  % draws x KDIM (scale slopes)
            sdid = cell2mat(rel.out.sd_id{1});       % draws x 4 (1-2 location, 3-4 scale)
            Lid = cell2mat(rel.out.L_id{1});         % draws x 4 x 4
            nd = size(bsg, 1);
            testCase.verifyEqual(size(sdid), [nd 4]);
            testCase.verifyEqual(size(Lid), [nd 4 4]);
            testCase.verifyEqual(size(idvc), [nd 2 2]);

            % L_id is a valid correlation Cholesky (unit row norms) on a sample draw
            Ld = squeeze(Lid(1,:,:));
            R = Ld*Ld.';
            testCase.verifyEqual(diag(R), ones(4,1), 'AbsTol', 1e-6);

            % ---- recovery: generous absolute bands on the population truth ----
            % difference universe-score variance (from the location block)
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 6.0);
            % per-event residual intercepts
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.35);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.35);
            % event-specific scale slopes (weakly identified; generous band)
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.30);
            % between-person residual-scale SDs (sd_id columns 3-4) -- the defining
            % SUBJECT-LEVEL parameter of this model
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.30);

            % ---- per-participant phi_delta,s(z) table: mean estimated phi vs the
            % closed-form truth, and every coefficient a proper reliability ----
            summ = psyrat_dynrel_summary(rel, 'CI', .95);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);

            genTrue = zeros(nsub,1);
            for p = 1:nsub
                wp1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1))^2;
                wp2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2))^2;
                relerr = wp1/ntrl + wp2/ntrl;
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.12);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end
    end
end
