classdef TestGammaDiffSserrReliability < matlab.unittest.TestCase
    %TESTGAMMADIFFSSERRRELIABILITY Offline proof of the Gamma (scaled-chi-square)
    % SUBJECT-LEVEL two-event DIFFERENCE-score math (analysis 8, ic_diff_sserrvar,
    % family='gamma', NON-concurrent). No CmdStan.
    %
    % The gamma case-8 model is IDENTICAL to the group-level gamma difference
    % (case 7, psyrat_build_stan_diff_chisq; byte-identical golden proven by
    % TestStanModelSnapshots/ic_gamma_diff_sserr). The subject-level distinction is
    % a downstream extraction (psyrat_gamma_extract_diff_sserr, a local function of
    % psyrat_computevarcomp) that, per event e, forms each participant's CONDITIONAL
    % residual from their own log-mean (u_pe = r_id(:,:,e)) and log-nu
    % (v_pe = r_id(:,:,e+2)),
    %     sigma_pi_e2(e,p) = 2*exp(2*(b[e]+u_pe) + 2*sd_trl[e]^2 - (b_nu[e]+v_pe)),
    % encodes it as a log-SD-equivalent (er_var_ss = 0.5*log(sigma_pi_e2)), and feeds
    % the SAME subject-level difference machinery (psyrat_ssrel_diff) the Gaussian
    % case-8 path uses. These tests prove the gamma-specific math the extraction
    % implements, mirroring TestGammaSserrConversion (the case-6 single-event proof);
    % the fit->extraction plumbing (indexing, field wiring) is validated by the live
    % CmdStan recovery in tests/integration/TestGammaDiffSserrRecovery.
    %
    % All draws use base-MATLAB randn only (no Statistics Toolbox): the per-event
    % residual 2*mu^2/nu_p is a function of (mu, nu_p) only, so mean(2*mu.^2 ./ nu_p)
    % is its exact ground truth with no observation-sampling noise.

    methods (Test)

        function testPerEventResidualMatchesConditionalMonteCarlo(testCase)
            % DECISIVE (per event). For fixed persons (u_p, v_p) each event's
            % per-subject residual formula must equal the Monte-Carlo within-person
            % expected observation variance E_t[2*mu_pt^2/nu_p] over the trial
            % effects, using THAT event's own (b, b_nu, sd_trl). The two events carry
            % distinct parameters, so a per-event parameter mix-up would fail here.
            evt(1) = struct('alpha',log(5.5),'beta',log(18),'s_i',0.30);
            evt(2) = struct('alpha',log(6.2),'beta',log(14),'s_i',0.22);
            persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00;  0.25 0.35 ];
            rng(8101);
            T = 4e6;
            for e = 1:2
                for k = 1:size(persons,1)
                    u_p = persons(k,1); v_p = persons(k,2);
                    nu_p = exp(evt(e).beta + v_p);
                    u_t  = evt(e).s_i * randn(T,1);
                    mu_pt = exp(evt(e).alpha + u_p + u_t);
                    mc  = mean(2 .* mu_pt.^2 ./ nu_p);                          % ground truth
                    fml = 2 * exp(2*(evt(e).alpha + u_p) + 2*evt(e).s_i^2 ...
                        - (evt(e).beta + v_p));                                 % formula
                    testCase.verifyEqual(fml, mc, 'RelTol', 5e-3, ...
                        sprintf('Event %d per-subject residual formula wrong.', e));
                end
            end
        end

        function testLogSdEncodingRoundTrips(testCase)
            % The extractor encodes each event's per-subject residual as a log-SD-
            % equivalent er_var_ss = 0.5*log(sigma_pi_e2) so psyrat_ssrel_diff's
            % wp = exp(er_var_ss).^2 reproduces sigma_pi_e2(e,p) exactly, and the
            % typical-subject (u=v=0) population residual anchors b_sigma.
            for e = 1:2
                alpha = log(5.0)+0.4*e; beta = log(16)+0.1*e; s_i = 0.25+0.05*e;
                persons = [ 0.60 -0.50;  -0.40 0.70;  0.00 0.00 ];
                for k = 1:size(persons,1)
                    u_p = persons(k,1); v_p = persons(k,2);
                    sig_pi_e2 = 2 * exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));
                    er_var_ss = 0.5 * log(sig_pi_e2);
                    testCase.verifyEqual(exp(er_var_ss)^2, sig_pi_e2, 'RelTol', 1e-13);
                end
                sig_pop = 2 * exp(2*alpha + 2*s_i^2 - beta);   % typical subject u=v=0
                b_sigma = 0.5 * log(sig_pop);
                testCase.verifyEqual(exp(b_sigma)^2, sig_pop, 'RelTol', 1e-13);
            end
        end

        function testNonConcurrentDiffReliabilityMatchesSsrelDiff(testCase)
            % INTEGRATION. Feeding the observed-scale SD-equivalents the extractor
            % emits (per-event person var-cov = id_varcov; per-subject residual log-
            % SDs = er_var_ss1/2; wp_cov_ss = 0) into the REAL subject-level
            % difference reliability path (psyrat_ssrel_diff) must reproduce the
            % Rocha difference-score coefficient
            %   D = Var(p_diff) / (Var(p_diff) + err_diff),
            % with Var(p_diff) = sp1 + sp2 - 2*cov_p and, NON-concurrent,
            %   err_diff = wp1/o1 + wp2/o2 (+ trial terms for dependability).
            % This proves the stored contract yields correct subject-level difference
            % reliability with the residual covariance held at 0.
            nd = 200; NSUB = 3;
            % per-draw observed-scale person (co)variance across the two events
            sp1 = 4.0 + 0.2*randn(nd,1);
            sp2 = 5.0 + 0.2*randn(nd,1);
            cov_p = -0.6 + 0.1*randn(nd,1);   % negative cross-event person covariance
            si1 = 0.15 + 0.02*randn(nd,1);    % observed trial variance, event 1
            si2 = 0.12 + 0.02*randn(nd,1);    % observed trial variance, event 2

            bp = zeros(nd,2,2);
            bp(:,1,1) = sp1; bp(:,2,2) = sp2; bp(:,1,2) = cov_p; bp(:,2,1) = cov_p;
            bt = zeros(nd,2,2);
            bt(:,1,1) = si1; bt(:,2,2) = si2;   % non-concurrent trial cross-cov = 0

            % per-subject residual log-SDs (er_var_ss = 0.5*log(sigma_pi_e2))
            wp1 = 1.2 + 0.05*randn(nd,NSUB);    % sigma_pi_e2 event 1 per subject
            wp2 = 1.6 + 0.05*randn(nd,NSUB);    % sigma_pi_e2 event 2 per subject
            er_var_ss1 = 0.5*log(wp1);
            er_var_ss2 = 0.5*log(wp2);
            b_sigma = [0.5*log(mean(wp1,2)), 0.5*log(mean(wp2,2))];  % population fallback

            id_match = table((1:NSUB)', (1:NSUB)', 'VariableNames', {'id','id2'});
            idtable = id_match;
            idtable.trls  = repmat(20, NSUB, 1);
            idtable.trls1 = repmat(20, NSUB, 1);   % o1
            idtable.trls2 = repmat(25, NSUB, 1);   % o2

            out = psyrat_ssrel_diff('bp',bp,'bt',bt,'er_var',b_sigma,...
                'wp_cov',zeros(nd,1),'idtable',idtable,'CI',0.95,'est','gen',...
                'er_var_ss1',er_var_ss1,'er_var_ss2',er_var_ss2,...
                'wp_cov_ss',zeros(nd,NSUB));

            % hand reference for each subject: generalizability (relative error)
            for r = 1:NSUB
                o1 = idtable.trls1(r); o2 = idtable.trls2(r);
                var_pdiff = sp1 + sp2 - 2*cov_p;
                rel_err   = wp1(:,r)./o1 + wp2(:,r)./o2;    % wp_cov = 0
                Gdraw = var_pdiff ./ (var_pdiff + rel_err);
                testCase.verifyEqual(out.rel_pt(r), mean(Gdraw), 'RelTol', 1e-10, ...
                    sprintf('Subject %d generalizability mismatch.', r));
            end
        end

    end
end
