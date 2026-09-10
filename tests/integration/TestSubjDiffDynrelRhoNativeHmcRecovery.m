classdef TestSubjDiffDynrelRhoNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the concurrent one-facet
    % subject-level dynamic difference model with a PER-SUBJECT residual
    % correlation (Rast & Clayson, analysis 21; dynrel=2 + diffest=2 + sserrvar=2 +
    % diffwpcov=2 + ssrescor=2). This is the case-19 model with the population
    % residual correlation replaced by a hierarchical per-participant correlation
    % rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]) (the case-8 rho hierarchy).
    % Matches psyrat_build_stan_diffdynrel_sserr_rho.
    %
    % Generator uses INDEPENDENT person event LOCATION effects (so the difference
    % universe variance recovers at modest N), a per-person residual-SCALE random
    % effect, and a per-person residual correlation drawn from the rescor_mu /
    % sd_rescor hierarchy. The per-subject correlations are individually weakly
    % identified (few pairs per subject), so the POPULATION-average correlation
    % (mean rho) is asserted, not the per-subject values.
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
            rng(40404042, 'twister');
            nsub = 35; ntrl = 15;
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rescorMu = atanh(0.4);            % population mean of atanh(rho)
            sdRescor = 0.3;                   % between-person SD of atanh(rho)
            sdP = [3 2.5];                    % person event LOCATION SDs (independent)
            sds = [0.5 0.4];                  % between-person residual-SCALE SDs
            trlSD = [1.0 0.8];
            uniTrue = sdP(1)^2 + sdP(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)];
            rho_s = tanh(rescorMu + sdRescor*randn(nsub,1));  % per-person correlation
            popRhoTrue = mean(rho_s);
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                    s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                    rp = rho_s(p);
                    Lr = [s1 0; rp*s2 sqrt(1-rp^2)*s2];
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
                'ssrescor', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_rho');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            rc   = cell2mat(rel.out.rescor{1});   % per-draw mean(rho)
            sdid = cell2mat(rel.out.sd_id{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(sdid), [nd 4]);

            % ---- recovery (generous bands) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            % population-average residual correlation (the per-subject values are
            % weakly identified, so assert the mean only, with a wide band)
            testCase.verifyEqual(mean(rc), popRhoTrue, 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            % the per-subject residual-scale SDs (sd_id columns 3-4)
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.30);

            % ---- per-participant reliability table is well-formed ----
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'obs', ntrl);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end
    end
end
