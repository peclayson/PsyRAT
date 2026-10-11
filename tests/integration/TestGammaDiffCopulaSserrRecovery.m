classdef TestGammaDiffCopulaSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL CONCURRENT two-event DIFFERENCE-score design (analysis 8,
    % ic_diff_sserrvar, diffrescor=2). The concurrent gamma case-8 model is the
    % case-7 Gaussian-copula model (psyrat_build_stan_diff_copula_chisq); the
    % subject-level distinction is the per-participant CONDITIONAL residual per event
    % PLUS the per-person copula residual cross-covariance wp_cov_ss
    % (psyrat_gamma_extract_diff_copula_sserr), fed through psyrat_ssrel_diff. This
    % test simulates copula-coupled paired observations with a bivariate person block
    % (correlated event means + per-event person-VARYING dispersion) and a KNOWN
    % copula correlation rho_e, and confirms:
    %   (1) a positive residual correlation is recovered (rho_e_true > 0);
    %   (2) the recovered per-subject residual cross-covariance wp_cov_ss is positive
    %       (the concurrency signal reaches the subject level);
    %   (3) the concurrent per-subject difference reliability EXCEEDS the
    %       non-concurrent counterpart from the SAME fit (folding wp_cov_ss into the
    %       residual shrinks the difference error);
    %   (4) the person log-nu SD s_v is recovered per event (the subject-level
    %       mechanism), and the cross-event person covariance (numerator off-diagonal)
    %       is recovered; and
    %   (5) the true per-subject concurrent generalizability coefficients fall inside
    %       the recovered 95% credible intervals (coverage), with the sample-mean
    %       coefficient recovered.
    %
    % SCORE, ESTIMAND, CONDITIONING, COUNTS. The score is the difference between
    % the two events' ntrl-trial means of trial-paired (copula-coupled) gamma
    % observations. The estimand is the per-person CONCURRENT generalizability
    %   G_p = S / (S + (E_{1,p} + E_{2,p} - 2 K_p)/ntrl),
    %   S = P_1 + P_2 - 2 CovP12,   E_{e,p} = M_e + C_{e,p},   K_p = E_i Cov(Y1,Y2 | p,i),
    % relative error. S, P_e, CovP12 and M_e are POPULATION quantities at the
    % event intercepts b(e); C_{e,p} and K_p are conditioned on the person's own
    % (u_pe, v_pe). The copula preserves the gamma marginals, so each event's
    % residual DIAGONAL is the non-concurrent one (S20 pooled: M_e + C_{e,p});
    % the population mean-interaction CROSS-covariance m_1 m_2 expm1(cP12) expm1(cI12)
    % is exactly zero here only because the trial effects are independent
    % across events (cI12 = 0), and it is evaluated rather than assumed. The
    % concurrent error uses the equal-count form (E1 + E2 - 2K)/ntrl, which is
    % what psyrat_ssrel_diff forms when trls1 = trls2 = ntrl; it is not
    % generalized to unequal counts. Allocation: nsub = 60 persons x 2 events x
    % ntrl = 40 trial-paired trials, n' = 40 for both events.
    %
    % Ground truth for the signal and residual diagonals is the independent 3A
    % reference (tests/helpers/PsyRAT3ARecoveryReference.m, pinned offline in
    % tests/TestGamma3ARecoveryReference.m); no production converter is called
    % for the truth. K_p is a DIRECT brute-force simulation of E_i[Cov(Y1,Y2|p,i)]
    % (unchanged: 150000 pairs per person, seed 24680; independent of the
    % extractor's mu*W factorization, which is separately validated offline in
    % TestGammaDiffCopulaSserrReliability). Its Monte-Carlo standard error is
    % derived from the influence-function variance of the sample covariance on
    % the same draws and propagated to G_conc through dG/dK = 2S/[n (S+err)^2];
    % it is REPORTED and stored with the truth, not used to change any threshold.
    % Before this correction the truth took P_e/CovP12 from
    % psyrat_gamma_varcomps/psyrat_gamma_crosscov and omitted M_e from the
    % residual diagonals; earlier logs remain historical evidence for that stale
    % target. Design, seed, sampler settings and every assertion tolerance are
    % unchanged. CmdStan-gated (skips Incomplete when CmdStan/MatlabStan absent).
    % When PSYRAT_2C_RECOVERY_CAPTURE_DIR is set, data and truth are captured
    % BEFORE the fit, a fit failure still captures whatever CmdStan wrote, the
    % complete run is captured BEFORE the convergence assumption, and the class
    % teardown skips its scratch cleanup. LONG (copula density + per-person MC);
    % run deliberately.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            % When evidence capture is requested, leave the raw scratch untouched
            % so it survives even a failed capture or a failed test; root's
            % external artifact guard archives/restores the owned scratch.
            if ~isempty(PsyRAT2CRecoveryCapture.captureDir())
                fprintf('[diff copula sserr recovery] %s set; scratch cleanup skipped.\n', ...
                    PsyRAT2CRecoveryCapture.EnvVar);
                return;
            end
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)
        function testConcurrentSubjectLevelDifferenceReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gaminv','file')) && ...
                logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            % ---- Known truth (log expected-score scale), per event ----
            b     = [log(5.5), log(6.0)];    % per-event log-mean intercepts
            b_nu  = [log(15),  log(13)];     % per-event log-nu intercepts
            s_p   = [0.40, 0.42];            % per-event person log-mean SDs
            s_v   = [0.55, 0.55];            % per-event person log-nu SDs
            s_i   = [0.20, 0.22];            % per-event trial log-mean SDs (independent across events)
            rho_mean   = 0.5;                % cross-event person-mean correlation (log scale)
            rho_pv     = -0.3;               % within-event mean<->log-nu correlation
            rho_e_true = 0.5;                % copula residual correlation (concurrent coupling)

            % Subject-level recovery needs the per-person dispersion (nu) identified,
            % which is driven by trials-per-person: the copula + per-person-nu funnel
            % is under-identified at ntrl=22 (max-treedepth non-convergence). Use the
            % non-concurrent case-8 scale (nsub=60, ntrl=40) + extra warmup so the
            % dispersion funnel adapts (owner decision 2026-07-14). LONG fit.
            nsub = 60; ntrl = 40;
            [dataTbl, u_p, v_p] = localSimulateGammaDiffCopulaSserr( ...
                nsub, ntrl, b, b_nu, s_p, s_v, s_i, rho_mean, rho_pv, rho_e_true);

            % ---- Independent closed-form observed-scale truth (3A reference) ----
            % SIGNAL (population, estimand-invariant): per event, B_e = exp(2 b(e)
            % + s_p(e)^2 + s_i(e)^2), P_e = B_e expm1(s_p(e)^2), and the observed
            % cross-event person covariance CovP12 = m_1 m_2 expm1(rho_mean s_p(1) s_p(2))
            % with m_e = exp(b(e) + (s_p(e)^2 + s_i(e)^2)/2). S = P_1 + P_2 - 2 CovP12.
            % The dispersion parameters (s_v, rho_pv) and the copula rho_e shape the
            % planted (u_p, v_p) and the pairing only; S does not depend on them.
            P_true = zeros(1,2); I_true = zeros(1,2); M_true = zeros(1,2); m_true = zeros(1,2);
            for e = 1:2
                comp = PsyRAT3ARecoveryReference.oneFacetMeanComponents(b(e), s_p(e), s_i(e));
                P_true(e) = comp.P; I_true(e) = comp.I; M_true(e) = comp.M;
                m_true(e) = PsyRAT3ARecoveryReference.eventMean(b(e), s_p(e), s_i(e));
            end
            cov_p_obs_true = PsyRAT3ARecoveryReference.personCrossCovariance( ...
                m_true(1), m_true(2), rho_mean * s_p(1) * s_p(2));
            var_pdiff_true = P_true(1) + P_true(2) - 2*cov_p_obs_true;

            % RESIDUAL DIAGONAL per event: the POOLED own residual wp_e = M_e + C_{e,p}
            % (M_e the population person x trial mean interaction at b(e), added
            % exactly once; C_{e,p} = 2 exp(2 (b(e) + u_pe) + 2 s_i(e)^2 - b_nu(e) - v_pe)
            % the person's own trial-averaged observation variance). The copula
            % preserves the marginals, so this is the non-concurrent diagonal.
            wp1_true = PsyRAT3ARecoveryReference.pooledOwnResidual(M_true(1), ...
                PsyRAT3ARecoveryReference.ownObservationVariance(b(1), u_p(:,1), s_i(1), b_nu(1), v_p(:,1)));
            wp2_true = PsyRAT3ARecoveryReference.pooledOwnResidual(M_true(2), ...
                PsyRAT3ARecoveryReference.ownObservationVariance(b(2), u_p(:,2), s_i(2), b_nu(2), v_p(:,2)));
            % The population mean-interaction CROSS-covariance m_1 m_2 expm1(cP12)
            % expm1(cI12) is exactly zero in THIS fixture because the trial
            % effects are independent across events (cI12 = 0). Evaluate it rather
            % than assume it, so adding the diagonal M_e terms provably needs no
            % mean-interaction off-diagonal here. Not a general rule.
            piCross_true = PsyRAT3ARecoveryReference.trialMeanCrossCovariance( ...
                m_true(1), m_true(2), rho_mean * s_p(1) * s_p(2), 0);
            testCase.assertEqual(piCross_true, 0, ...
                'Population mean-interaction cross-covariance must be exactly zero at cI12 = 0.');

            % RESIDUAL OFF-DIAGONAL: the per-person CONCURRENT observation
            % cross-covariance K_p = E_i[Cov(Y1,Y2|p,i)] by DIRECT simulation
            % (unchanged distribution and design: 1.5e5 pairs per person, seed
            % 24680), plus its Monte-Carlo standard error from the influence
            % function of the sample covariance on the same draws.
            [wpcov_true, wpcov_mcse] = localTrueConcResidCrossCov( ...
                b, b_nu, s_i, u_p, v_p, rho_e_true, 1.5e5);
            % concurrent: (wp1 + wp2 - 2 K)/ntrl (equal counts; what psyrat_ssrel_diff
            % forms with trls1 = trls2 = ntrl); non-concurrent diagnostic: (wp1 + wp2)/ntrl
            trueResid_conc = PsyRAT3ARecoveryReference.differenceError( ...
                wp1_true, ntrl, wp2_true, ntrl, wpcov_true);
            trueResid_nonc = PsyRAT3ARecoveryReference.differenceError( ...
                wp1_true, ntrl, wp2_true, ntrl, 0);
            trueG_conc = PsyRAT3ARecoveryReference.differenceCoefficient( ...
                var_pdiff_true, trueResid_conc);                                % nsub x 1
            trueG_nonc = PsyRAT3ARecoveryReference.differenceCoefficient( ...
                var_pdiff_true, trueResid_nonc);
            % First-order propagation of the K Monte-Carlo error to G_conc:
            % dG/dK = 2 S / [ntrl (S + err)^2]. Reported and stored; the assertion
            % thresholds below are unchanged.
            trueG_conc_mcse = PsyRAT3ARecoveryReference.dGdK( ...
                var_pdiff_true, trueResid_conc, ntrl) .* wpcov_mcse;
            fprintf(['[diff copula sserr recovery] truth MC uncertainty: median K=%.4f ', ...
                '(median mcse %.2e, max %.2e); max G_conc mcse=%.2e vs the unchanged ', ...
                'sample-mean tolerance .06\n'], median(wpcov_true), median(wpcov_mcse), ...
                max(wpcov_mcse), max(trueG_conc_mcse));
            truth = struct('b', b, 'b_nu', b_nu, 's_p', s_p, 's_v', s_v, 's_i', s_i, ...
                'rho_mean', rho_mean, 'rho_pv', rho_pv, 'rho_e', rho_e_true, ...
                'nsub', nsub, 'ntrl', ntrl, ...
                'P', P_true, 'I', I_true, 'M', M_true, 'm', m_true, ...
                'cov_p_obs', cov_p_obs_true, 'signal', var_pdiff_true, ...
                'pi_cross_cov', piCross_true, 'u_p', u_p, 'v_p', v_p, ...
                'E_p', [wp1_true, wp2_true], 'K_p', wpcov_true, 'K_p_mcse', wpcov_mcse, ...
                'error_conc_p', trueResid_conc, 'error_nonc_p', trueResid_nonc, ...
                'G_conc_p', trueG_conc, 'G_conc_p_mcse', trueG_conc_mcse, 'G_nonc_p', trueG_nonc, ...
                'K_simulation', 'direct: 150000 (trial, copula) pairs per person, rng(24680), sample covariance', ...
                'K_mcse_method', 'influence function: sqrt(var((x-xbar).*(y-ybar))/N); dG/dK = 2S/[n (S+err)^2]', ...
                'residual_definition', ['error_conc_p = (M_1 + C_{1,p} + M_2 + C_{2,p} - 2 K_p)/ntrl ', ...
                '(population per-event mean interaction plus own observation variance, minus ', ...
                'own concurrent cross-covariance; equal counts); G_p = S/(S + error_p), ', ...
                'S = P_1 + P_2 - 2 CovP12'], ...
                'generator', 'Gaussian copula on continuous gamma marginals, seed 12345');

            % ---- Fit (case 8 gamma, CONCURRENT copula), with evidence capture ----
            % (a no-op unless the capture root is set; see localCaptureStage)
            fitTag = 'recov_gamma_diff_sserr_copula';
            outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), fitTag);
            captureMeta = struct('test', 'TestGammaDiffCopulaSserrRecovery/testConcurrentSubjectLevelDifferenceReliabilityRecovered', ...
                'analysis', 'ic_diff_sserrvar', 'diffest', 2, 'diffwpcov', 2, 'chains', 2, ...
                'warmup', 1250, 'sampling', 750, 'seed', 12345, 'nsub', nsub, 'ntrl', ntrl, ...
                'design', 'two-event concurrent (Gaussian copula) difference, persons x trials, subject-level');
            localCaptureStage('static8_lognu_diff_copula_sserr_prefit', [], truth, dataTbl, outDir, ...
                captureMeta, 'prefit');
            try
                rel = localRunDiffCopulaSserrGamma(testCase, dataTbl, fitTag);
            catch fitErr
                localCaptureStage('static8_lognu_diff_copula_sserr_fitfailed', ...
                    struct('fit_error', getReport(fitErr, 'extended', 'hyperlinks', 'off')), ...
                    truth, dataTbl, outDir, captureMeta, 'fitfailed');
                rethrow(fitErr);
            end
            localCaptureStage('static8_lognu_diff_copula_sserr', rel, truth, dataTbl, outDir, ...
                captureMeta, 'postfit');
            testCase.verifyEqual(rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % --- Convergence gate: the recovery checks below are only meaningful on a
            % CONVERGED fit (matching truth from an R-hat>=1.1 posterior is coincidental).
            % Gate via assume so a non-converged copula fit is filtered (Incomplete),
            % not a false PASS. The copula per-person-nu geometry is known not to converge
            % under this short-warmup regime (S11); the same R-hat>=1.1 + n_eff criteria
            % as the production psyrat_checkconv path decide the gate.
            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[diff copula sserr recovery] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Copula fit ', ...
                'non-converged (max R-hat=%.3f >= 1.1 / n_eff gate); recovery not validated ', ...
                'in this environment (S11 copula per-person-nu convergence limit).'], maxRhat));

            % ---- Recovered per-subject difference reliability (mirrors psyrat_relsummary) ----
            idvar  = cell2mat(rel.out.id_varcov{1});     % nd x 2 x 2
            trlvar = cell2mat(rel.out.trl_varcov{1});    % nd x 2 x 2
            bsigma = cell2mat(rel.out.b_sigma{1});       % nd x 2
            er1    = rel.out.er_var_ss1{1};              % nd x nsub
            er2    = rel.out.er_var_ss2{1};              % nd x nsub
            wpss   = rel.out.wp_cov_ss{1};               % nd x nsub (concurrent, nonzero)
            id_match = rel.out.id_matches{1};

            idtable = id_match;
            idtable.trls1 = repmat(ntrl, height(idtable), 1);
            idtable.trls2 = repmat(ntrl, height(idtable), 1);
            idtable.trls  = min(idtable.trls1, idtable.trls2);

            ndp = size(bsigma,1);
            args = {'bp',idvar,'bt',trlvar,'er_var',bsigma,'wp_cov',zeros(ndp,1),...
                'idtable',idtable,'CI',0.95,'est','gen','er_var_ss1',er1,'er_var_ss2',er2};
            tblConc = psyrat_ssrel_diff(args{:},'wp_cov_ss',wpss);
            tblNonc = psyrat_ssrel_diff(args{:},'wp_cov_ss',zeros(ndp,size(wpss,2)));

            % align recovered per-subject quantities to truth via the id string ('S###')
            nrows = height(tblConc);
            recG = tblConc.rel_pt; recLL = tblConc.rel_ll; recUL = tblConc.rel_ul;
            tGc = zeros(nrows,1); tGn = zeros(nrows,1);
            for r = 1:nrows
                idx = sscanf(char(string(tblConc.id(r))), 'S%d');
                tGc(r) = trueG_conc(idx);
                tGn(r) = trueG_nonc(idx);
            end

            rescor   = cell2mat(rel.out.rescor{1});
            svHat    = mean(rel.out.gamma_disp.s_v, 1);
            covHat   = mean(idvar(:,1,2));
            meanWpss = mean(wpss(:));
            % rescor is the OBSERVED residual correlation, which is the gamma-marginal
            % attenuation of the latent copula rho_e (=0.50 here -> observed ~0.2), so
            % it is reported against rho_e but not asserted equal to it.
            fprintf(['[diff copula sserr recovery] rescor hat=%.3f (obs; copula rho_e=%.2f)  ', ...
                'mean wp_cov_ss=%.3f  s_v hat=[%.3f %.3f]  cov_p obs hat=%.3f (true %.3f)  ', ...
                'coverage=%.2f  mean recG=%.3f (true conc %.3f, nonc %.3f)\n'], ...
                mean(rescor), rho_e_true, meanWpss, svHat(1), svHat(2), covHat, ...
                cov_p_obs_true, mean(tGc >= recLL & tGc <= recUL), mean(recG), ...
                mean(tGc), mean(tGn));

            % (1) positive residual correlation recovered (the copula effect; the
            % OBSERVED correlation is attenuated below the latent rho_e by the gamma
            % marginals, so a modest positive threshold, not rho_e itself)
            testCase.verifyGreaterThan(mean(rescor), 0.05, ...
                'A positive copula residual correlation should be recovered (rho_e_true>0).');

            % (2) the per-subject residual cross-cov reaches the subject level (positive)
            testCase.verifyGreaterThan(meanWpss, 0, ...
                'Recovered per-subject residual cross-covariance should be positive.');

            % (3) concurrency RAISES subject-level reliability vs the non-concurrent form
            testCase.verifyGreaterThan(mean(tblConc.rel_pt), mean(tblNonc.rel_pt), ...
                'Concurrent per-subject difference reliability must exceed the non-concurrent form.');

            % (4) s_v + cross-event person covariance recovered
            testCase.verifyGreaterThan(min(svHat), 0.35, sprintf( ...
                'Recovered person log-nu SD [%.3f %.3f] far below true %.2f.', ...
                svHat(1), svHat(2), s_v(1)));
            testCase.verifyLessThan(max(svHat), 0.85);
            testCase.verifyGreaterThan(covHat, 0, ...
                'Recovered cross-event person covariance should be positive.');
            testCase.verifyLessThan(abs(covHat - cov_p_obs_true), 0.5*abs(cov_p_obs_true) + 1, ...
                sprintf('Cross-event person covariance off: recovered %.3f vs true %.3f.', ...
                covHat, cov_p_obs_true));

            % (5) COVERAGE: most true CONCURRENT per-subject coefficients inside the 95%
            % CIs, and the sample-mean coefficient recovered.
            coverage = mean(tGc >= recLL & tGc <= recUL);
            testCase.verifyGreaterThanOrEqual(coverage, 0.80, sprintf( ...
                'Only %.0f%% of true per-subject concurrent coefficients inside the 95%% CIs.', ...
                100*coverage));
            testCase.verifyLessThan(abs(mean(recG) - mean(tGc)), 0.06, sprintf( ...
                'Sample-mean concurrent per-subject reliability off: recovered %.3f vs true %.3f.', ...
                mean(recG), mean(tGc)));
        end
    end
end

function [tbl, u_p, v_p] = localSimulateGammaDiffCopulaSserr(nsub, ntrl, b, b_nu, ...
    s_p, s_v, s_i, rho_mean, rho_pv, rho_e)
% Copula-coupled two-event Persons x Trials long table: a bivariate log-linked mean
% (correlated event means + independent crossed per-event trial main effects) and a
% person-VARYING per-event scaled-chi-square dispersion. Each (person, trial) yields
% a Gaussian-copula-coupled pair (Y1, Y2) at correlation rho_e, emitted in E1/E2 runs
% paired by trial position (as psyrat_build_diff_pairs re-pairs them).
rng(12345);
R = eye(4);
R(1,2) = rho_mean; R(2,1) = rho_mean;
R(1,3) = rho_pv;   R(3,1) = rho_pv;
R(2,4) = rho_pv;   R(4,2) = rho_pv;
sd = [s_p(1), s_p(2), s_v(1), s_v(2)];
L = chol(R, 'lower');
Z = randn(4, nsub);
blk = (diag(sd) * L * Z)';          % nsub x 4
u_p = blk(:,1:2);                   % person mean effects, per event
v_p = blk(:,3:4);                   % person log-nu effects, per event

u_t = [s_i(1) * randn(ntrl,1), s_i(2) * randn(ntrl,1)];   % crossed trial effects (independent)

nrow = nsub * 2 * ntrl;
ids  = cell(nrow,1); evs = cell(nrow,1); meas = zeros(nrow,1);
r = 0; evnames = {'E1','E2'};
for p = 1:nsub
    nu1 = exp(b_nu(1) + v_p(p,1));
    nu2 = exp(b_nu(2) + v_p(p,2));
    z1 = randn(ntrl,1);
    z2 = rho_e .* z1 + sqrt(1 - rho_e^2) .* randn(ntrl,1);
    mu1 = exp(b(1) + u_p(p,1) + u_t(:,1));
    mu2 = exp(b(2) + u_p(p,2) + u_t(:,2));
    Y1 = gaminv(normcdf(z1), 0.5*nu1, 2*mu1/nu1);
    Y2 = gaminv(normcdf(z2), 0.5*nu2, 2*mu2/nu2);
    for e = 1:2
        for t = 1:ntrl
            r = r + 1;
            ids{r} = sprintf('S%03d', p);
            evs{r} = evnames{e};
            if e == 1, meas(r) = Y1(t); else, meas(r) = Y2(t); end
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});
end

function [wpcov, mcse] = localTrueConcResidCrossCov(b, b_nu, s_i, u_p, v_p, rho_e, NBIG)
% Direct brute-force truth for the per-person conditional residual cross-covariance
% E_i[Cov(Y1,Y2|p,i)]. Trial effects are INDEPENDENT across events here, so the
% trial-mean cross-cov Cov_i(mu1,mu2) = 0 and the sample covariance of (Y1,Y2) over
% many (trial, copula) draws estimates E_i[Cov(Y1,Y2|p,i)] directly. Independent of
% the extractor's mu*W factorization.
%
% The simulation (distribution, NBIG, seed, RNG call order) is unchanged. The
% second output is the Monte-Carlo standard error of each person's sample
% covariance from its influence function (PsyRAT3ARecoveryReference.covarianceMcse),
% evaluated on the SAME draws; it consumes no random numbers, so wpcov is the
% historical value bit for bit.
rng(24680);
nsub = size(u_p,1);
wpcov = zeros(nsub,1);
mcse  = zeros(nsub,1);
for p = 1:nsub
    nu1 = exp(b_nu(1) + v_p(p,1));
    nu2 = exp(b_nu(2) + v_p(p,2));
    ut1 = s_i(1) * randn(NBIG,1);
    ut2 = s_i(2) * randn(NBIG,1);           % independent trial effects across events
    mu1 = exp(b(1) + u_p(p,1) + ut1);
    mu2 = exp(b(2) + u_p(p,2) + ut2);
    z1 = randn(NBIG,1);
    z2 = rho_e .* z1 + sqrt(1 - rho_e^2) .* randn(NBIG,1);
    Y1 = gaminv(normcdf(z1), 0.5*nu1, 2*mu1/nu1);
    Y2 = gaminv(normcdf(z2), 0.5*nu2, 2*mu2/nu2);
    c = cov(Y1, Y2);
    wpcov(p) = c(1,2);
    mcse(p) = PsyRAT3ARecoveryReference.covarianceMcse(Y1, Y2);
end
end

function rel = localRunDiffCopulaSserrGamma(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'diffwpcov', 2, ...       % -> diffrescor=2 (concurrent copula)
    'chains', 2, ...
    'warmup', 1250, ...       % extra warmup for the per-person-nu / copula funnel
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCaptureStage(tag, rel, truth, dataTbl, outDir, meta, stage)
% Evidence capture through the existing PsyRAT2CRecoveryCapture API (a no-op
% unless PSYRAT_2C_RECOVERY_CAPTURE_DIR is set). Three stages, three tags:
%   prefit    - data + truth only (no REL, no CmdStan output yet). The helper
%               writes the MAT/CSV/JSON first and then reports the raw capture
%               as INCOMPLETE because no chain CSV exists; that is the expected
%               outcome here and is tolerated, so data and truth are on disk
%               before CmdStan starts.
%   fitfailed - the fit threw: preserve the error report with whatever CmdStan
%               wrote (possibly partial or no chain CSVs); an incomplete raw
%               capture is again tolerated, and the caller rethrows the fit error.
%   postfit   - the complete REL + raw output, BEFORE any assertion; a capture
%               failure here is a real failure and propagates (sibling pattern).
% Any other capture error (non-absolute root, overlap with the scratch, an
% existing destination) propagates at every stage. Nothing is deleted.
if isempty(PsyRAT2CRecoveryCapture.captureDir())
    return;
end
meta.stage = stage;
try
    PsyRAT2CRecoveryCapture.capture(tag, rel, truth, dataTbl, outDir, 'meta', meta);
catch capErr
    if strcmp(stage, 'postfit') || ...
            ~strcmp(capErr.identifier, 'PsyRAT2CRecoveryCapture:rawCaptureIncomplete')
        rethrow(capErr);
    end
    fprintf(['[3A recovery capture] %s (%s): data/truth preserved; raw CmdStan capture ', ...
        'incomplete as expected at this stage (%s)\n'], tag, stage, capErr.message);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'recov_gamma_diff_sserr_copula');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            'Could not remove recovery artifact directory ''%s''.', p);
    end
end
pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
for i = 1:numel(pats)
    files = dir(fullfile(baseDir, pats{i}));
    for j = 1:numel(files)
        if ~files(j).isdir
            delete(fullfile(baseDir, files(j).name));
        end
    end
end
dirs = dir(fullfile(baseDir, 'cmdstan_*'));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
for i = 1:numel(dirs)
    try
        rmdir(fullfile(baseDir, dirs(i).name), 's');
    catch
        % best-effort cleanup
    end
end
end
