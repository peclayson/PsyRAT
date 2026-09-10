classdef TestDiffTrtNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the STATIC two-facet
    % (test-retest) two-event difference reliability model (analysis 10,
    % NON-concurrent: diffest=2 + sserrvar=1 + a time/occasion facet, default
    % diffwpcov=1). This is the KDIM=0 reduction of the dynamic two-facet
    % difference (case 15): the same six crossed 2x2 cross-condition blocks
    % (id/trl/occ/tid/oid/to), each a 2x2 per-event Cholesky, with NO dimension
    % slopes. The two conditions are stacked long with 0/1 event indicators and
    % a group-level per-event residual SD (no residual covariance). Matches
    % psyrat_build_stan_diff_trt; reuses local_hmc_diffdynrel_trt at KDIM=0.
    %
    % As in the case-15 recovery test (TestDiffDynrelTrtNativeHmcRecovery), the
    % generator uses INDEPENDENT person/facet event effects (diagonal
    % covariances) so the LKJ(2) prior AGREES with the truth and the difference
    % universe-score variance sigma^2_p_delta = bp1 + bp2 recovers at a modest
    % subject count. The well-identified components and the static difference
    % G coefficient are asserted; the few-level facets only with broad bounds.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. Six crossed factors
    % make this a high-dimensional model, so it is slow and lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (INDEPENDENT person/facet event effects; constant
            % per-event residual SD, no dimension) ----
            rng(20260627, 'twister');
            nsub = 28; nocc = 3; ntrl = 4;
            mu = [5 2];
            a = [log(3) log(2.5)];           % per-event residual log-SDs (constant)
            sdP   = [3 2.5];                  % person event LOCATION SDs (independent)
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = sdP(1)^2 + sdP(2)^2;   % sigma^2_p_delta (independent events)
            bpiTrue = ptSD(1)^2 + ptSD(2)^2; % person x trial difference variance
            bpoTrue = poSD(1)^2 + poSD(2)^2; % person x occasion difference variance

            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; events = {'gain','loss'};
            for p = 1:nsub
                for e = 1:2
                    for o = 1:nocc
                        for t = 1:ntrl
                            mn = mu(e) + nu_p(p,e) + occE(o,e) + trlE(t,e) ...
                                + ptE(p,t,e) + poE(p,o,e) + toE(t,o,e);
                            sg = exp(a(e));
                            ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                            ev{end+1,1} = events{e};            %#ok<AGROW>
                            tm{end+1,1} = sprintf('t%d', o);    %#ok<AGROW>
                            meas(end+1,1) = mn + randn*sg;      %#ok<AGROW>
                        end
                    end
                end
            end
            T = table(ids, meas, ev, tm, ...
                'VariableNames', {'id', 'meas', 'event', 'time'});

            % ---- run native HMC through the full pipeline (diffest=2 + time
            % facet -> analysis 10; default diffwpcov=1 -> non-concurrent) ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'diffest', 2, ...
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
            nd = size(bsg, 1);
            testCase.verifyEqual(size(idvc), [nd 2 2]);
            testCase.verifyEqual(size(bsg), [nd 2]);

            % ---- well-identified components (generous bands) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);

            % person x trial / person x occasion difference variances
            bpiHat = mean(tidvc(:,1,1) + tidvc(:,2,2) - 2*tidvc(:,2,1));
            bpoHat = mean(oidvc(:,1,1) + oidvc(:,2,2) - 2*oidvc(:,2,1));
            testCase.verifyEqual(bpiHat, bpiTrue, 'AbsTol', max(1.5, 0.5*bpiTrue));
            testCase.verifyEqual(bpoHat, bpoTrue, 'AbsTol', max(1.5, 0.5*bpoTrue));

            % ---- static two-facet difference G recovered vs truth closed form
            % (non-concurrent -> residual covariance is 0) ----
            nt = ntrl; no = nocc;
            ds = psyrat_diffrel_trt('bp', idvc, 'bpi', tidvc, 'bpo', oidvc, ...
                'bt', trlvc, 'bo', occvc, 'boi', tovc, ...
                'er_var', bsg, 'er_cov', zeros(nd,1), ...
                'obs', [nt nt], 'nocc', [no no], 'reltype', 3, 'est', 'gen', 'CI', .95);
            res = exp(a(1))^2 + exp(a(2))^2;     % Var(e1 - e2), independent residuals
            Gtr = uniTrue / (uniTrue + bpiTrue/nt + bpoTrue/no + res/(nt*no));
            testCase.verifyEqual(ds.pt, Gtr, 'AbsTol', 0.12);
        end
    end
end
