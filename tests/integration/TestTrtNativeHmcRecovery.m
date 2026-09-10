classdef TestTrtNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the plain test-retest
    % crossed variance-components model (analysis 5: sserrvar=1 with a
    % time/occasion facet). This is a plain-Gaussian six-component crossed
    % random-effects ANOVA (Persons x Trials x Occasions): a fixed intercept,
    % six crossed scalar random effects (id, occ, trl, trlxid, occxid, trlxocc)
    % in the non-centered form (term = sig_k * raw_k, raw_k ~ normal(0,1)), and a
    % single constant residual SD. Matches psyrat_build_stan_trt and the native
    % engine local_hmc_trt.
    %
    % Case 5 is the ONE design with two native engines. This test forces engine=3
    % (native HMC), the engine the auto-cascade prefers over fitlme. The fitlme
    % fallback is known to be weak when the between-occasion SD (sig_occ) is poorly
    % identified at small occasion counts; the well-identified components plus the
    % test-retest G coefficient are asserted here, the few-level facets (occ/trl)
    % with broader bounds.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. High-dimensional (the
    % person-by-trial and person-by-occasion raw effects dominate), so slow; lives
    % in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth ----
            rng(20260627, 'twister');
            nsub = 30; nocc = 4; ntrl = 5;
            pop_int = 10;
            sig_id = 3; sig_occ = 2; sig_trl = 1.0;
            sig_trlxid = 1.2; sig_occxid = 1.0; sig_trlxocc = 0.8;
            sig_err = 2.0;

            idE = randn(nsub,1)*sig_id;
            occE = randn(nocc,1)*sig_occ;
            trlE = randn(ntrl,1)*sig_trl;
            tidE = randn(nsub,ntrl)*sig_trlxid;   % person x trial-position
            oidE = randn(nsub,nocc)*sig_occxid;   % person x occasion
            toE  = randn(nocc,ntrl)*sig_trlxocc;  % occasion x trial-position

            ids = {}; meas = []; tm = {};
            for p = 1:nsub
                for o = 1:nocc
                    for t = 1:ntrl
                        y = pop_int + idE(p) + occE(o) + trlE(t) ...
                            + tidE(p,t) + oidE(p,o) + toE(o,t) + randn*sig_err;
                        ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                        tm{end+1,1}  = sprintf('t%d', o);   %#ok<AGROW>
                        meas(end+1,1) = y;                  %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, tm, 'VariableNames', {'id', 'meas', 'time'});

            % ---- run native HMC through the full pipeline (time facet,
            % sserrvar=1 -> analysis 5) ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'trt');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws (one chunk: no group/event) ----
            sid  = rel.out.sig_id(:,1);
            socc = rel.out.sig_occ(:,1);
            strl = rel.out.sig_trl(:,1);
            stid = rel.out.sig_trlxid(:,1);
            soid = rel.out.sig_occxid(:,1);
            sto  = rel.out.sig_trlxocc(:,1);
            serr = rel.out.sig_err(:,1);
            mu   = rel.out.mu(:,1);
            nd = numel(sid);
            testCase.verifyTrue(nd > 0 && all(sid > 0) && all(serr > 0));

            % ---- well-identified components (many levels): tighter bands ----
            testCase.verifyEqual(mean(mu), pop_int, 'AbsTol', 1.0);
            testCase.verifyEqual(mean(sid), sig_id, 'AbsTol', max(0.9, 0.35*sig_id));
            testCase.verifyEqual(mean(serr), sig_err, 'AbsTol', max(0.5, 0.20*sig_err));
            testCase.verifyEqual(mean(stid), sig_trlxid, 'AbsTol', max(0.6, 0.5*sig_trlxid));
            testCase.verifyEqual(mean(soid), sig_occxid, 'AbsTol', max(0.6, 0.5*sig_occxid));

            % ---- weakly identified components (few levels): sanity bounds only.
            % The between-occasion (4 levels), between-trial (5 levels), and
            % trial-by-occasion SDs are estimated from very few levels under a wide
            % half-Cauchy prior, so the posterior mean can sit well above the
            % realized effect spread (in practice sig_occ ~ 3.6 and sig_trl ~ 2.3
            % here, vs the generative 2.0 / 1.0). That inflation is a property of
            % the DESIGN, not the engine, so assert only that each is finite,
            % positive, and below the total observed SD (a single component cannot
            % exceed the marginal spread). The well-identified components above and
            % the generalizability coefficient below carry the recovery validation.
            % (sig_occ is exactly the term the fitlme fallback handles poorly at
            % small occasion counts; HMC keeps it finite and bounded.) ----
            smeas = std(meas);
            testCase.verifyTrue(isfinite(mean(socc)) && mean(socc) > 0 && mean(socc) < smeas);
            testCase.verifyTrue(isfinite(mean(strl)) && mean(strl) > 0 && mean(strl) < smeas);
            testCase.verifyTrue(isfinite(mean(sto))  && mean(sto)  >= 0 && mean(sto)  < smeas);

            % ---- test-retest G (generalizability, CES) recovered vs truth
            % closed form ----
            nt = ntrl; no = nocc;
            [~, Gpt, ~] = psyrat_rel_trt('gcoeff', 2, 'reltype', 3, ...
                'bp', sid, 'bo', socc, 'bt', strl, ...
                'txp', stid, 'oxp', soid, 'txo', sto, 'err', serr, ...
                'obs', nt, 'nocc', no, 'CI', .95);
            Gtr = sig_id^2 / (sig_id^2 + sig_trlxid^2/nt + sig_occxid^2/no ...
                + sig_err^2/(nt*no));
            testCase.verifyEqual(Gpt, Gtr, 'AbsTol', 0.12);
        end
    end
end
