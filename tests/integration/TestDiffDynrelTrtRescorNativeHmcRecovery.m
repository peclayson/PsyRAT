classdef TestDiffDynrelTrtRescorNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the CONCURRENT two-facet
    % group-level dynamic difference reliability model (Rast & Clayson, analysis
    % 18; dynrel=2 + diffest=2 + diffwpcov=2 + a time/occasion facet). This is the
    % concurrent one-facet model (case 17) with the person+trial 2-D blocks
    % generalized to the six crossed cross-condition 2-D factors (id/trl/occ/tid/
    % oid/to, as case 15); same per-row bivariate residual (Lrescor) + per-event
    % dimension slopes. Matches psyrat_build_stan_diffdynrel_trt_rescor.
    %
    % The CmdStan path is validated by TestDiffDynrelRescorRecovery
    % (testTwoFacetConcurrentRecovers, nsub=72). The generator here uses
    % INDEPENDENT person/facet event LOCATION effects (so the LKJ(2) prior agrees
    % with the truth and the difference universe variance recovers at a modest
    % subject count), while keeping the bivariate residual correlation. The
    % well-identified components and the relative G_delta(z) surface are asserted.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. High-dimensional (six
    % crossed factors), so slow.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (independent person/facet event location effects;
            % bivariate residual with correlation rho) ----
            rng(50505050, 'twister');
            nsub = 25; nocc = 3; ntrl = 4;
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rho = 0.4;
            sdP = [3 2.5];                    % person event LOCATION SDs (independent)
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = sdP(1)^2 + sdP(2)^2;
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for o = 1:nocc
                    for t = 1:ntrl
                        s1 = exp(a(1) + bsig(1)*zsub(p));
                        s2 = exp(a(2) + bsig(2)*zsub(p));
                        Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2];
                        e = Lr * randn(2,1);
                        cm1 = nu_p(p,1) + occE(o,1) + trlE(t,1) + ptE(p,t,1) + poE(p,o,1) + toE(t,o,1);
                        cm2 = nu_p(p,2) + occE(o,2) + trlE(t,2) + ptE(p,t,2) + poE(p,o,2) + toE(t,o,2);
                        m1 = mu(1) + bmean(1)*zsub(p) + cm1 + e(1);
                        m2 = mu(2) + bmean(2)*zsub(p) + cm2 + e(2);
                        ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                        ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                        tm  = [tm; {sprintf('t%d',o)}; {sprintf('t%d',o)}];       %#ok<AGROW>
                        meas = [meas; m1; m2];                                   %#ok<AGROW>
                        dim1 = [dim1; zsub(p); zsub(p)];                         %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, tm, dim1, ...
                'VariableNames', {'id','meas','event','time','dim1'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'dynrel', 2, ...
                'diffest', 2, ...
                'diffwpcov', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_trt_rescor');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc  = cell2mat(rel.out.id_varcov{1});
            tidvc = cell2mat(rel.out.tid_varcov{1});
            oidvc = cell2mat(rel.out.oid_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            bsd = cell2mat(rel.out.b_sigma_dim{1});
            rc  = cell2mat(rel.out.rescor{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(idvc), [nd 2 2]);
            testCase.verifyEqual(size(rc), [nd 1]);

            % ---- well-identified components (generous bands) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.22);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.22);

            bpiHat = mean(tidvc(:,1,1) + tidvc(:,2,2) - 2*tidvc(:,2,1));
            bpoHat = mean(oidvc(:,1,1) + oidvc(:,2,2) - 2*oidvc(:,2,1));
            testCase.verifyEqual(bpiHat, bpiTrue, 'AbsTol', max(1.5, 0.5*bpiTrue));
            testCase.verifyEqual(bpoHat, bpoTrue, 'AbsTol', max(1.5, 0.5*bpoTrue));

            % ---- two-facet difference G(z) recovered vs truth closed form (with
            % the residual covariance term) ----
            zg = [-1.5 0 1.5]; nt = ntrl; no = nocc;
            out = psyrat_rel_diffdynrel_trt('id_varcov', idvc, ...
                'trl_varcov', cell2mat(rel.out.trl_varcov{1}), ...
                'occ_varcov', cell2mat(rel.out.occ_varcov{1}), ...
                'tid_varcov', tidvc, 'oid_varcov', oidvc, ...
                'to_varcov', cell2mat(rel.out.to_varcov{1}), ...
                'b_sigma', bsg, 'b_sigma_dim', bsd, 'rescor', rc, ...
                'ndim', 1, 'z1', zg, 'obs', nt, 'nocc', no, 'reltype', 3, 'CI', .95);
            for k = 1:numel(zg)
                s1 = exp(a(1) + bsig(1)*zg(k)); s2 = exp(a(2) + bsig(2)*zg(k));
                res = s1^2 + s2^2 - 2*rho*s1*s2;        % difference residual var
                Gtr = uniTrue / (uniTrue + bpiTrue/nt + bpoTrue/no + res/(nt*no));
                testCase.verifyEqual(out.G.pt(k), Gtr, 'AbsTol', 0.12);
            end
        end
    end
end
