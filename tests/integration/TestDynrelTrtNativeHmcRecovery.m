classdef TestDynrelTrtNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the trial + occasion
    % two-facet dynamic/conditional reliability model (Rast & Clayson, analysis
    % 14; dynrel = 2 with a time/occasion facet). This is the one-facet dynrel
    % locscale model with the single trial main effect generalized to five scalar
    % crossed-facet mean random effects (occ, trl, trial x person, occasion x
    % person, trial x occasion). Matches psyrat_build_stan_dynrel_trt.
    %
    % The CmdStan path for this exact design is validated separately by
    % TestDynrelTrtRecovery; this test confirms the native HMC engine recovers the
    % SAME known truth without CmdStan, through the full psyrat_computevarcomp
    % pipeline (the five-facet data prep, the 2-D person location-scale block, and
    % the two-facet G(z)/D(z) surface via psyrat_rel_dynrel_trt). The data
    % generator is the shared balanced person x occasion x trial simulator.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. The person x trial /
    % person x occasion raw effects make this a high-dimensional model, so it is
    % slow and lives in tests/integration. Data are kept modest (vs the CmdStan
    % test's nsub50/nocc4/ntrl6) to bound runtime; the well-identified person-level
    % components and the relative G(z) surface are asserted, the few-level facets
    % only with broad sanity bounds (as in the CmdStan test).

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth ----
            rng(20260622, 'twister');
            nsub = 35; nocc = 4; ntrl = 5;
            grand = 15;
            sigP = 4; sigDP = 0.25;
            sigOcc = 1.5; sigTrl = 1.5; sigTxP = 1.0; sigOxP = 1.0; sigTxO = 0.8;
            sigLog0 = log(3); bMean = 1.5; bSigma = 0.3;

            zraw = randn(nsub, 1);
            zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = sigP * randn(nsub, 1);
            del_p = sigDP * randn(nsub, 1);
            occ_main = sigOcc * randn(nocc, 1);
            trl_main = sigTrl * randn(ntrl, 1);
            pt = sigTxP * randn(nsub, ntrl);
            po = sigOxP * randn(nsub, nocc);
            to = sigTxO * randn(ntrl, nocc);

            n = nsub * nocc * ntrl;
            ids = cell(n, 1); tm = cell(n, 1);
            meas = zeros(n, 1); dim1 = zeros(n, 1);
            r = 0;
            for pp = 1:nsub
                for oo = 1:nocc
                    for tt = 1:ntrl
                        r = r + 1;
                        mu = grand + bMean * zsub(pp) + nu_p(pp) ...
                            + occ_main(oo) + trl_main(tt) ...
                            + pt(pp, tt) + po(pp, oo) + to(tt, oo);
                        sg = exp(sigLog0 + bSigma * zsub(pp) + del_p(pp));
                        ids{r} = sprintf('S%02d', pp);
                        tm{r} = sprintf('t%d', oo);
                        meas(r) = mu + sg * randn;
                        dim1(r) = zsub(pp);
                    end
                end
            end
            dataTbl = table(ids, meas, tm, dim1, ...
                'VariableNames', {'id', 'meas', 'time', 'dim1'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'dynrel', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_dynrel_trt');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- posterior means ----
            gro = cell2mat(rel.out.gro_sds{1});        % draws x 2 (sigma_p, sigma_dp)
            log0 = rel.out.pop_sdlog(:, 1);
            bsig = cell2mat(rel.out.b_sigma{1});       % draws x 1
            sOcc = rel.out.sig_occ(:, 1);
            sTrl = rel.out.sig_trl(:, 1);
            sTxP = rel.out.sig_trlxid(:, 1);
            sOxP = rel.out.sig_occxid(:, 1);
            sTxO = rel.out.sig_trlxocc(:, 1);

            % ---- recovery: well-identified components (generous bands) ----
            testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 1.5);
            testCase.verifyEqual(mean(log0), sigLog0, 'AbsTol', 0.35);
            testCase.verifyEqual(mean(bsig(:, 1)), bSigma, 'AbsTol', 0.20);
            testCase.verifyEqual(mean(gro(:, 2)), sigDP, 'AbsTol', 0.35);
            testCase.verifyEqual(mean(sTxP), sigTxP, 'AbsTol', 0.8);
            testCase.verifyEqual(mean(sOxP), sigOxP, 'AbsTol', 0.8);

            % ---- few-level facets: broad sanity bounds only ----
            testCase.verifyGreaterThan(mean(sOcc), 0);
            testCase.verifyLessThan(mean(sOcc), sigOcc + 3);
            testCase.verifyGreaterThan(mean(sTrl), 0);
            testCase.verifyLessThan(mean(sTrl), sigTrl + 3);
            testCase.verifyGreaterThan(mean(sTxO), 0);
            testCase.verifyLessThan(mean(sTxO), sigTxO + 3);

            % ---- two-facet G(z)/D(z) recovered vs truth closed form (reltype 3) ----
            zg = [-1.5 0 1.5]; nt = ntrl; no = nocc;
            out = psyrat_rel_dynrel_trt('sig_p', gro(:, 1), 'sig_log0', log0, ...
                'b_sigma', bsig(:, 1), 'sig_occ', sOcc, 'sig_trl', sTrl, ...
                'sig_trlxid', sTxP, 'sig_occxid', sOxP, 'sig_trlxocc', sTxO, ...
                'ndim', 1, 'z1', zg, 'obs', nt, 'nocc', no, 'reltype', 3, ...
                'CI', .95);
            for a = 1:numel(zg)
                sige = exp(sigLog0 + bSigma * zg(a));
                uni = sigP^2;
                gerr = sigTxP^2/nt + sigOxP^2/no + sige^2/(nt*no);
                derr = gerr + sigTrl^2/nt + sigOcc^2/no + sigTxO^2/(nt*no);
                Gtr = uni / (uni + gerr);
                Dtr = uni / (uni + derr);
                testCase.verifyEqual(out.G.pt(a), Gtr, 'AbsTol', 0.10);
                testCase.verifyEqual(out.D.pt(a), Dtr, 'AbsTol', 0.15);
            end
        end
    end
end
