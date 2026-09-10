classdef TestDiffSserrNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the STATIC subject-level
    % two-event difference-score model (analysis 8; diffest=2 + sserrvar=2 +
    % diffwpcov=2). This is the last HMC design: three separate 2-D blocks (person
    % location, trial, per-subject scale) + a per-subject residual correlation
    % rho[s]=tanh(rescor_mu+sd_rescor*z_rescor[s]) over PAIRED (bivariate) and
    % UNPAIRED (univariate) observations. Matches psyrat_build_stan_diff_sserr.
    %
    % Generator: balanced concurrent gain/loss with INDEPENDENT person event
    % location effects (so the difference universe variance recovers at modest N),
    % a per-person residual-SCALE random effect, and a per-person residual
    % correlation from the rescor_mu / sd_rescor hierarchy. The per-subject
    % correlations are weakly identified, so the population-average correlation is
    % asserted.
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
            rng(80808080, 'twister');
            nsub = 35; ntrl = 12;
            b = [5 2];
            a = [log(3) log(2.5)];           % per-event residual log-SD intercept
            rescorMu = atanh(0.4); sdRescor = 0.3;
            sdLoc = [3 2.5];                 % person event LOCATION SDs (independent)
            sdSc  = [0.5 0.4];               % between-person residual-SCALE SDs
            sdTl  = [1.0 0.9];               % trial location SDs
            uniTrue = sdLoc(1)^2 + sdLoc(2)^2;   % difference universe variance

            nu_p   = [randn(nsub,1)*sdLoc(1), randn(nsub,1)*sdLoc(2)];
            delta  = [randn(nsub,1)*sdSc(1), randn(nsub,1)*sdSc(2)];  % person scale RE
            rho_s  = tanh(rescorMu + sdRescor*randn(nsub,1));
            popRhoTrue = mean(rho_s);
            trlL   = [randn(ntrl,1)*sdTl(1), randn(ntrl,1)*sdTl(2)];

            ids = {}; meas = []; ev = {}; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + delta(p,1));
                    s2 = exp(a(2) + delta(p,2));
                    rp = rho_s(p);
                    Lr = [s1 0; rp*s2 sqrt(1-rp^2)*s2];
                    e = Lr * randn(2,1);
                    m1 = b(1) + nu_p(p,1) + trlL(t,1) + e(1);
                    m2 = b(2) + nu_p(p,2) + trlL(t,2) + e(2);
                    ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                    ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                    meas = [meas; m1; m2];                                   %#ok<AGROW>
                end
            end
            T = table(ids, meas, ev, 'VariableNames', {'id','meas','event'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
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

            testCase.verifyEqual(rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc = cell2mat(rel.out.id_varcov{1});   % draws x 2 x 2 (location var-cov)
            sdid = cell2mat(rel.out.sd_id{1});       % draws x 2 (location SDs)
            rc   = cell2mat(rel.out.rescor{1});      % draws x 1 (mean rho)
            er1  = rel.out.er_var_ss1{1};            % draws x NSUB (per-subject log-resid e1)
            wp   = rel.out.wp_cov_ss{1};             % draws x NSUB
            nd = size(idvc, 1);
            testCase.verifyEqual(size(idvc), [nd 2 2]);
            testCase.verifyEqual(size(er1), [nd nsub]);
            testCase.verifyEqual(size(wp), [nd nsub]);

            % ---- recovery (generous bands) ----
            % person-cell location variances + the difference universe variance
            testCase.verifyEqual(mean(idvc(:,1,1)), sdLoc(1)^2, 'AbsTol', max(3.0, 0.5*sdLoc(1)^2));
            testCase.verifyEqual(mean(idvc(:,2,2)), sdLoc(2)^2, 'AbsTol', max(3.0, 0.5*sdLoc(2)^2));
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(6.0, 0.4*uniTrue));
            % location SDs
            testCase.verifyEqual(mean(sdid(:,1)), sdLoc(1), 'AbsTol', 0.8);
            testCase.verifyEqual(mean(sdid(:,2)), sdLoc(2), 'AbsTol', 0.8);
            % population-average residual correlation (weakly identified -> wide band)
            testCase.verifyEqual(mean(rc), popRhoTrue, 'AbsTol', 0.25);
            % per-event residual log-SD intercept (the per-subject scale REs average
            % to ~0, so mean over subjects/draws of er_var_ss ~ a)
            testCase.verifyEqual(mean(er1(:)), a(1), 'AbsTol', 0.30);
        end
    end
end
