classdef TestSubjDiffDynrelRescorNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the CONCURRENT one-facet
    % SUBJECT-LEVEL dynamic difference reliability model (Rast & Clayson, analysis
    % 19; dynrel=2 + diffest=2 + sserrvar=2 + diffwpcov=2). This is the concurrent
    % one-facet group model (case 17) with the person block widened to a 4-D joint
    % location+scale block (a 4x4 LKJ Cholesky, as case 13), so each participant
    % has their own per-event residual SDs (the per-subject scale REs enter the
    % per-row sigma); the residual correlation (Lrescor) stays population. Matches
    % psyrat_build_stan_diffdynrel_sserr_rescor.
    %
    % Generator uses INDEPENDENT person event LOCATION effects (so the difference
    % universe variance recovers at modest N) plus the per-person residual-scale
    % random effect (the defining subject-level quantity) and the bivariate
    % residual correlation. One-facet (no crossed factors), so moderate-dimensional.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth ----
            rng(40404041, 'twister');
            nsub = 35; ntrl = 15;
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rho = 0.4;
            sdP = [3 2.5];                    % person event LOCATION SDs (independent)
            sds = [0.5 0.4];                  % between-person residual-SCALE SDs
            trlSD = [1.0 0.8];
            uniTrue = sdP(1)^2 + sdP(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)]; % person scale RE
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                    s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                    Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2];
                    e = Lr * randn(2,1);
                    m1 = mu(1) + bmean(1)*zsub(p) + nu_p(p,1) + trlE(t,1) + e(1);
                    m2 = mu(2) + bmean(2)*zsub(p) + nu_p(p,2) + trlE(t,2) + e(2);
                    ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                    ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                    meas = [meas; m1; m2];                                   %#ok<AGROW>
                    dim1 = [dim1; zsub(p); zsub(p)];                         %#ok<AGROW>
                end
            end
            T = table(ids, meas, ev, dim1, ...
                'VariableNames', {'id','meas','event','dim1'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'dynrel', 2, ...
                'diffest', 2, ...
                'sserrvar', 2, ...
                'diffwpcov', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_rescor');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            rc   = cell2mat(rel.out.rescor{1});
            sdid = cell2mat(rel.out.sd_id{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(sdid), [nd 4]);
            testCase.verifyEqual(size(rc), [nd 1]);

            % ---- recovery: difference universe variance + headline rescor +
            % the per-subject residual-scale SDs (sd_id columns 3-4) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.22);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.22);
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.30);

            % ---- per-participant reliability table is well-formed (in (0,1)) ----
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'obs', ntrl);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end
    end
end
