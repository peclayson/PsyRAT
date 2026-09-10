classdef TestDiffDynrelRescorNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the CONCURRENT one-facet
    % group-level dynamic difference reliability model (Rast & Clayson, analysis
    % 17; dynrel=2 + diffest=2 + diffwpcov=2). This is the bivariate-residual
    % difference model (case 7*) with per-event fixed dimension slopes on the mean
    % and the per-event log-residual SD, so the bivariate residual Cholesky is
    % per-row (dimension-dependent). Matches psyrat_build_stan_diffdynrel_rescor.
    %
    % The CmdStan path is validated by TestDiffDynrelRescorRecovery (nsub=72). The
    % generator here uses INDEPENDENT person event LOCATION effects (so the LKJ(2)
    % prior agrees with the truth and the difference universe variance recovers at
    % a modest subject count), while keeping the bivariate residual correlation
    % (the headline rescor parameter) and the dimension-dependent residual SDs.
    % This is a one-facet model (no crossed factors), so it is the lowest-
    % dimensional concurrent variant.
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

            % ---- known truth (independent person event location effects; bivariate
            % residual with correlation rho) ----
            rng(40404040, 'twister');
            nsub = 35; ntrl = 15;
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rho = 0.4;                            % TRUE residual correlation
            sdP = [3 2.5];                        % person event LOCATION SDs (independent)
            trlSD = [1.0 0.8];
            uniTrue = sdP(1)^2 + sdP(2)^2;        % sigma^2_p_delta (independent events)

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p));
                    s2 = exp(a(2) + bsig(2)*zsub(p));
                    Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2]; % chol of residual cov
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
                'diffwpcov', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 300, ...
                'sampling', 300, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_rescor');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            rc   = cell2mat(rel.out.rescor{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(idvc), [nd 2 2]);
            testCase.verifyEqual(size(rc), [nd 1]);

            % ---- recovery: difference universe variance + the headline rescor ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            % interleaved b_sigma_dim: col 1 = event 1, col 2 = event 2
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.22);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.22);

            % ---- dynamic surface G_delta(z) vs the closed form WITH the residual
            % covariance term ----
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'obs', ntrl);
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                s1 = exp(a(1) + bsig(1)*zv); s2 = exp(a(2) + bsig(2)*zv);
                res = s1^2 + s2^2 - 2*rho*s1*s2;       % difference residual var
                Gtr = uniTrue / (uniTrue + res/ntrl);
                testCase.verifyEqual(st.G.pt(idx), Gtr, 'AbsTol', 0.12);
            end
            testCase.verifyTrue(all(st.G.pt(:) > 0 & st.G.pt(:) < 1));
        end
    end
end
