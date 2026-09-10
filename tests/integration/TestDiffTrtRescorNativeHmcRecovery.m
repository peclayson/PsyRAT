classdef TestDiffTrtRescorNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the STATIC CONCURRENT
    % two-facet (test-retest) two-event difference reliability model (analysis
    % 10, concurrent: diffest=2 + sserrvar=1 + diffwpcov=2 + a time/occasion
    % facet). This is the KDIM=0 reduction of the dynamic concurrent two-facet
    % difference (case 18): the six crossed cross-condition 2-D factors (id/trl/
    % occ/tid/oid/to) in concurrent (paired) form with NO dimension slopes, plus
    % a per-row bivariate residual (Lrescor) that collapses to a constant
    % bivariate Cholesky when KDIM=0. Matches psyrat_build_stan_diff_trt_rescor;
    % reuses local_hmc_diffdynrel_trt_rescor at KDIM=0.
    %
    % As in the case-18 recovery test (TestDiffDynrelTrtRescorNativeHmcRecovery),
    % the generator uses INDEPENDENT person/facet event LOCATION effects (so the
    % LKJ(2) prior agrees with the truth and the difference universe variance
    % recovers at a modest subject count), while keeping the bivariate residual
    % correlation. The well-identified components, the residual correlation, and
    % the static difference G coefficient are asserted.
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
            % constant bivariate residual with correlation rho, no dimension) ----
            rng(60606060, 'twister');
            nsub = 25; nocc = 3; ntrl = 4;
            mu = [5 2];
            a = [log(3) log(2.5)];           % per-event residual log-SDs (constant)
            rho = 0.4;
            sdP = [3 2.5];                    % person event LOCATION SDs (independent)
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = sdP(1)^2 + sdP(2)^2;
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;

            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            s1 = exp(a(1)); s2 = exp(a(2));
            Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2];   % constant bivariate Cholesky

            ids = {}; meas = []; ev = {}; tm = {}; events = {'gain','loss'};
            for p = 1:nsub
                for o = 1:nocc
                    for t = 1:ntrl
                        e = Lr * randn(2,1);
                        cm1 = nu_p(p,1) + occE(o,1) + trlE(t,1) + ptE(p,t,1) + poE(p,o,1) + toE(t,o,1);
                        cm2 = nu_p(p,2) + occE(o,2) + trlE(t,2) + ptE(p,t,2) + poE(p,o,2) + toE(t,o,2);
                        m1 = mu(1) + cm1 + e(1);
                        m2 = mu(2) + cm2 + e(2);
                        ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                        ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                        tm  = [tm; {sprintf('t%d',o)}; {sprintf('t%d',o)}];       %#ok<AGROW>
                        meas = [meas; m1; m2];                                   %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, tm, ...
                'VariableNames', {'id','meas','event','time'});

            % ---- run native HMC through the full pipeline (diffest=2 + time
            % facet -> analysis 10; diffwpcov=2 -> concurrent) ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'diffest', 2, ...
                'diffwpcov', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'trt_diff');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc  = cell2mat(rel.out.id_varcov{1});
            tidvc = cell2mat(rel.out.tid_varcov{1});
            oidvc = cell2mat(rel.out.oid_varcov{1});
            trlvc = cell2mat(rel.out.trl_varcov{1});
            occvc = cell2mat(rel.out.occ_varcov{1});
            tovc  = cell2mat(rel.out.to_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            rc  = cell2mat(rel.out.rescor{1});
            evc = cell2mat(rel.out.err_varcov{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(idvc), [nd 2 2]);
            testCase.verifyEqual(size(rc), [nd 1]);
            testCase.verifyEqual(size(evc), [nd 2 2]);

            % ---- well-identified components (generous bands) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);

            bpiHat = mean(tidvc(:,1,1) + tidvc(:,2,2) - 2*tidvc(:,2,1));
            bpoHat = mean(oidvc(:,1,1) + oidvc(:,2,2) - 2*oidvc(:,2,1));
            testCase.verifyEqual(bpiHat, bpiTrue, 'AbsTol', max(1.5, 0.5*bpiTrue));
            testCase.verifyEqual(bpoHat, bpoTrue, 'AbsTol', max(1.5, 0.5*bpoTrue));

            % err_varcov off-diagonal is rescor reconstructed at the constant
            % per-event residual SDs (the term psyrat_trtdiff_ercov uses).
            erCovHat = mean(evc(:,2,1));
            testCase.verifyEqual(erCovHat, rho*s1*s2, 'AbsTol', max(0.8, 0.5*abs(rho*s1*s2)));

            % ---- static two-facet difference G recovered vs truth closed form
            % (concurrent -> residual covariance subtracted) ----
            nt = ntrl; no = nocc;
            er_cov = rc(:) .* exp(bsg(:,1)) .* exp(bsg(:,2));   % as psyrat_trtdiff_ercov
            ds = psyrat_diffrel_trt('bp', idvc, 'bpi', tidvc, 'bpo', oidvc, ...
                'bt', trlvc, 'bo', occvc, 'boi', tovc, ...
                'er_var', bsg, 'er_cov', er_cov, ...
                'obs', [nt nt], 'nocc', [no no], 'reltype', 3, 'est', 'gen', 'CI', .95);
            res = s1^2 + s2^2 - 2*rho*s1*s2;     % Var(e1 - e2) with covariance
            Gtr = uniTrue / (uniTrue + bpiTrue/nt + bpoTrue/no + res/(nt*no));
            testCase.verifyEqual(ds.pt, Gtr, 'AbsTol', 0.12);
        end
    end
end
