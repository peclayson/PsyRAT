classdef TestDynrelNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the one-facet dynamic/
    % conditional reliability model (Rast & Clayson, analysis 11; dynrel = 2).
    %
    % This is the native-engine counterpart to TestDynrelRecovery (which uses
    % CmdStan). It drives the SAME case-11 pipeline through psyrat_computevarcomp
    % with engine = 3, so the model goes
    %   case 11 -> psyrat_run_stan -> psyrat_resolve_engine -> psyrat_run_native
    %   -> psyrat_native_hmc (matching psyrat_build_stan_dynrel) -> the case-11
    %   extraction into REL.out.
    % It requires NO CmdStan: the native HMC engine is pure MATLAB (hmcSampler +
    % automatic-differentiation gradients via dlarray). It needs the Statistics
    % and Machine Learning Toolbox (hmcSampler) and the Deep Learning Toolbox
    % (dlarray); when either is absent the test skips cleanly (Incomplete, via an
    % assumption) rather than failing. It is slow (samples several hundred draws
    % with AD gradients) and lives in tests/integration so the unit lane excludes
    % it.
    %
    % The accuracy oracle proves only "unchanged"; this proves the native engine
    % recovers known variance components and a sensible G(z)/D(z) surface.

    methods (Test)
        function testOneDimensionRecoversTruth(testCase)
            % skip cleanly unless both required toolboxes are installed
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (one standardized dimension) ----
            rng(20260626, 'twister');
            nsub = 40; ntrl = 12;
            grand = 4.0; sigLog0 = log(1.2);
            bMean = 0.8; bSigma = 0.35;
            sigP = 1.5; sigDP = 0.40; rho = 0.5; sigI = 0.7;

            zraw = randn(nsub, 1);
            zsub = (zraw - mean(zraw)) / std(zraw);     % standardize once
            S = [sigP^2, rho*sigP*sigDP; rho*sigP*sigDP, sigDP^2];
            L = chol(S, 'lower');
            e = L * randn(2, nsub);
            nu_p = e(1, :)';        % person location effect
            del_p = e(2, :)';       % person log-residual effect
            beta_t = sigI * randn(ntrl, 1);

            ids = cell(nsub * ntrl, 1);
            meas = zeros(nsub * ntrl, 1);
            dim1 = zeros(nsub * ntrl, 1);
            r = 0;
            for pp = 1:nsub
                for tt = 1:ntrl
                    r = r + 1;
                    mu = grand + bMean * zsub(pp) + nu_p(pp) + beta_t(tt);
                    sg = exp(sigLog0 + bSigma * zsub(pp) + del_p(pp));
                    ids{r} = sprintf('S%02d', pp);
                    meas(r) = mu + sg * randn;
                    dim1(r) = zsub(pp);
                end
            end
            dataTbl = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});

            % ---- run the native HMC engine through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', dataTbl, ...
                'dynrel', 2, ...
                'engine', 3, ...        % force native HMC
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            % ---- provenance + routing ----
            testCase.verifyEqual(rel.analysis, 'ic_dynrel');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws in the CmdStan-shaped layout ----
            gro = cell2mat(rel.out.gro_sds{1});      % draws x 2
            log0 = rel.out.pop_sdlog(:, 1);          % Intercept_sigma
            bsig = cell2mat(rel.out.b_sigma{1});     % draws x 1 (KDIM = 1)
            sigi = rel.out.sig_trl(:, 1);
            cholmat = cell2mat(rel.out.chol_corrmat{1}); % draws x 2 x 2
            nd = size(gro, 1);
            testCase.verifyEqual(size(gro), [nd 2]);
            testCase.verifyEqual(size(bsig), [nd 1]);
            testCase.verifyEqual(size(cholmat), [nd 2 2]);
            testCase.verifyEqual(size(cell2mat(rel.out.ind_bs{1})), [nd nsub]);

            % ---- recovery: posterior mean within a generous absolute band of
            % the population truth (the same philosophy as the CmdStan
            % TestDynrelRecovery; the bands are wide enough to absorb the small-N
            % realized-vs-population sampling gap noted above and any cross-machine
            % MCMC variation, while still catching a gross engine regression). rho
            % is weakly identified at this N and is deliberately NOT asserted.
            % Reference recovery at this fixed config/seed (bit-reproducible):
            % sigma_p 1.70, sigma_dp 0.52, Intercept_sigma 0.24, b_sigma 0.30,
            % sig_trl 1.01 (sig_trl is weakly identified with only 12 trials).
            testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 0.6);   % sigma_p
            testCase.verifyEqual(mean(gro(:, 2)), sigDP, 'AbsTol', 0.3);  % sigma_dp
            testCase.verifyEqual(mean(log0), sigLog0, 'AbsTol', 0.25);    % Intercept_sigma
            testCase.verifyEqual(mean(bsig), bSigma, 'AbsTol', 0.2);      % scale slope
            testCase.verifyEqual(mean(sigi), sigI, 'AbsTol', 0.7);        % sigma_i (trial)

            % ---- dynamic reliability surface recovers the closed-form truth ----
            zg = [-1.5 0 1.5]; nprime = ntrl;
            out = psyrat_rel_dynrel('sig_p', gro(:, 1), 'sig_log0', log0, ...
                'b_sigma', bsig, 'sig_i', sigi, 'ndim', 1, 'z1', zg, ...
                'obs', nprime, 'CI', .95);
            for a = 1:numel(zg)
                sige = exp(sigLog0 + bSigma * zg(a));
                Gtr = sigP^2 / (sigP^2 + sige^2 / nprime);
                Dtr = sigP^2 / (sigP^2 + (sige^2 + sigI^2) / nprime);
                testCase.verifyEqual(out.G.pt(a), Gtr, 'AbsTol', 0.06);
                testCase.verifyEqual(out.D.pt(a), Dtr, 'AbsTol', 0.06);
            end

            % ---- and the surface is well-formed across a finer grid ----
            zgrid = linspace(-2, 2, 9);
            outf = psyrat_rel_dynrel('sig_p', gro(:, 1), 'sig_log0', log0, ...
                'b_sigma', bsig, 'sig_i', sigi, 'ndim', 1, 'z1', zgrid, ...
                'obs', ntrl, 'CI', .95);
            testCase.verifyTrue(all(outf.D.pt >= 0 & outf.D.pt <= 1), ...
                'D(z) must lie in [0,1]');
            testCase.verifyTrue(all(outf.G.pt >= 0 & outf.G.pt <= 1), ...
                'G(z) must lie in [0,1]');
            % dependability includes sigma_i, so D(z) <= G(z) everywhere
            testCase.verifyTrue(all(outf.D.pt <= outf.G.pt + 1e-9), ...
                'D(z) must not exceed G(z)');
            % residual SD grows with z (b_sigma > 0), so reliability decreases
            testCase.verifyGreaterThan(outf.D.pt(1), outf.D.pt(end));
        end
    end
end
