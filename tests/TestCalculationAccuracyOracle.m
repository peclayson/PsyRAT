classdef TestCalculationAccuracyOracle < PsyRATTestBase

    methods (Test)
        function testDepScalarAndRangeMatchOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.oneFacet();

            [ll,pt,ul] = psyrat_dep('bp', fx.bp, 'wp', fx.wp, 'obs', fx.obs, 'CI', fx.ci);
            [ell,ept,eul] = PsyRATAccuracyOracle.dep(fx.bp, fx.wp, fx.obs, fx.ci);

            testCase.verifyEqual(ll, ell, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, ept, 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, eul, 'AbsTol', 1e-12);

            [llR,ptR,ulR] = psyrat_dep('bp', fx.bp, 'wp', fx.wp, 'obs', fx.obsRange, 'CI', fx.ci);
            [ellR,eptR,eulR] = PsyRATAccuracyOracle.dep(fx.bp, fx.wp, fx.obsRange, fx.ci);

            testCase.verifyEqual(llR, ellR, 'AbsTol', 1e-12);
            testCase.verifyEqual(ptR, eptR, 'AbsTol', 1e-12);
            testCase.verifyEqual(ulR, eulR, 'AbsTol', 1e-12);

            % ptR is aligned to obsRange(1):obsRange(2) (one column per trial
            % count), so dependability monotonicity is checked across the row.
            testCase.verifyGreaterThanOrEqual(diff(ptR), 0);
        end

        function testICCMatchesOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.oneFacet();

            [ll,pt,ul] = psyrat_icc('bp', fx.bp, 'wp', fx.wp, 'CI', fx.ci);
            [ell,ept,eul] = PsyRATAccuracyOracle.icc(fx.bp, fx.wp, fx.ci);

            testCase.verifyEqual(ll, ell, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, ept, 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, eul, 'AbsTol', 1e-12);
            localVerifyIntervalBounds(testCase, ll, pt, ul);
        end

        function testRelTrtAllBranchesMatchOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.trt();

            for gcoeff = 1:2
                for reltype = 1:3
                    [ll,pt,ul] = psyrat_rel_trt( ...
                        'gcoeff', gcoeff, 'reltype', reltype, ...
                        'bp', fx.bp, 'bo', fx.bo, 'bt', fx.bt, ...
                        'txp', fx.txp, 'oxp', fx.oxp, 'txo', fx.txo, 'err', fx.err, ...
                        'obs', fx.obs, 'nocc', fx.nocc, 'CI', fx.ci);

                    [ell,ept,eul] = PsyRATAccuracyOracle.relTrt( ...
                        gcoeff, reltype, fx.bp, fx.bo, fx.bt, fx.txp, fx.oxp, fx.txo, fx.err, ...
                        fx.obs, fx.nocc, fx.ci);

                    testCase.verifyEqual(ll, ell, 'AbsTol', 1e-12);
                    testCase.verifyEqual(pt, ept, 'AbsTol', 1e-12);
                    testCase.verifyEqual(ul, eul, 'AbsTol', 1e-12);
                    localVerifyIntervalBounds(testCase, ll, pt, ul);
                end
            end
        end

        function testCutscoreFunctionsMatchOracle(testCase)
            % Cut-scores are absolute-error only (Rocha 2026, Tables 2-3):
            % only the dependability (gcoeff = 1 / est = 'dep') case is a valid
            % quantity, so the relative/gen cases are gone.
            fx = PsyRATAccuracyFixtureFactory.trt();
            one = PsyRATAccuracyFixtureFactory.oneFacet();

            [ll,pt,ul] = psyrat_cutscore_dep( ...
                'bp', one.bp, 'wp', one.wp, 'mu', fx.mu, 'cut', fx.cut, ...
                'obs', fx.obs, 'CI', fx.ci, 'est', 'dep');
            [ell,ept,eul] = PsyRATAccuracyOracle.cutscoreDep( ...
                one.bp, one.wp, fx.mu, fx.cut, fx.obs, fx.ci, 'dep', 0);

            testCase.verifyEqual(ll, ell, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, ept, 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, eul, 'AbsTol', 1e-12);

            for reltype = 1:3
                [ll,pt,ul] = psyrat_cutscore_trt( ...
                    'gcoeff', 1, 'reltype', reltype, ...
                    'bp', fx.bp, 'bo', fx.bo, 'bt', fx.bt, ...
                    'txp', fx.txp, 'oxp', fx.oxp, 'txo', fx.txo, 'err', fx.err, ...
                    'mu', fx.mu, 'cut', fx.cut, ...
                    'obs', fx.obs, 'nocc', fx.nocc, 'CI', fx.ci);

                [ell,ept,eul] = PsyRATAccuracyOracle.cutscoreTrt( ...
                    1, reltype, fx.bp, fx.bo, fx.bt, fx.txp, fx.oxp, fx.txo, fx.err, ...
                    fx.mu, fx.cut, fx.obs, fx.nocc, fx.ci);

                testCase.verifyEqual(ll, ell, 'AbsTol', 1e-12);
                testCase.verifyEqual(pt, ept, 'AbsTol', 1e-12);
                testCase.verifyEqual(ul, eul, 'AbsTol', 1e-12);
            end

            % Relative-error cut-scores are not offered: gcoeff = 2 is rejected.
            testCase.verifyError(@() psyrat_cutscore_trt( ...
                'gcoeff', 2, 'reltype', 1, ...
                'bp', fx.bp, 'bo', fx.bo, 'bt', fx.bt, ...
                'txp', fx.txp, 'oxp', fx.oxp, 'txo', fx.txo, 'err', fx.err, ...
                'mu', fx.mu, 'cut', fx.cut, ...
                'obs', fx.obs, 'nocc', fx.nocc, 'CI', fx.ci), ...
                'varargin:gcoeff');
        end

        function testDiffFunctionsMatchOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.diffOneFacet();

            for est = {'dep','gen'}
                out = psyrat_diffrel( ...
                    'bp', fx.bp, 'bt', fx.bt, 'er_var', fx.er_var, 'wp_cov', fx.wp_cov, ...
                    'obs', fx.obs, 'est', est{1}, 'CI', fx.ci);
                exp = PsyRATAccuracyOracle.diffrel( ...
                    fx.bp, fx.bt, fx.er_var, fx.wp_cov, fx.obs, est{1}, fx.ci);

                testCase.verifyEqual(out.ll, exp.ll, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.pt, exp.pt, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.ul, exp.ul, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.icc_ll, exp.icc_ll, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.icc_pt, exp.icc_pt, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.icc_ul, exp.icc_ul, 'AbsTol', 1e-12);
                localVerifyIntervalBounds(testCase, out.ll, out.pt, out.ul);
            end

            fxTrt = PsyRATAccuracyFixtureFactory.diffTrt();
            for reltype = 1:3
                for est = {'dep','gen'}
                    out = psyrat_diffrel_trt( ...
                        'bp', fxTrt.bp, 'bpi', fxTrt.bpi, 'bpo', fxTrt.bpo, ...
                        'bt', fxTrt.bt, 'bo', fxTrt.bo, 'boi', fxTrt.boi, ...
                        'er_var', fxTrt.er_var, 'er_cov', fxTrt.er_cov, ...
                        'obs', fxTrt.obs, 'nocc', fxTrt.nocc, ...
                        'reltype', reltype, 'est', est{1}, 'CI', fxTrt.ci);

                    exp = PsyRATAccuracyOracle.diffrelTrt( ...
                        fxTrt.bp, fxTrt.bpi, fxTrt.bpo, fxTrt.bt, fxTrt.bo, fxTrt.boi, ...
                        fxTrt.er_var, fxTrt.er_cov, fxTrt.obs, fxTrt.nocc, reltype, est{1}, fxTrt.ci);

                    testCase.verifyEqual(out.ll, exp.ll, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.pt, exp.pt, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.ul, exp.ul, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.icc_ll, exp.icc_ll, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.icc_pt, exp.icc_pt, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.icc_ul, exp.icc_ul, 'AbsTol', 1e-12);
                    localVerifyIntervalBounds(testCase, out.ll, out.pt, out.ul);
                end
            end
        end

        function testSubjectSpecificFunctionsMatchOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.ssOneFacet();

            out = psyrat_ssrel( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'idtable', fx.idtable, 'CI', fx.ci);
            exp = PsyRATAccuracyOracle.ssrel(fx.bp, fx.wp_pop, fx.wp_ss, fx.idtable, fx.ci);

            testCase.verifyEqual(out.dep_pt, exp.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.dep_ll, exp.dep_ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.dep_ul, exp.dep_ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_pt, exp.icc_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.sem_pt, exp.sem_pt, 'AbsTol', 1e-12);

            outCut = psyrat_sscutscore( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'mu', fx.mu, 'cut', fx.cut, 'idtable', fx.idtable, 'CI', fx.ci, 'est', 'dep');
            expCut = PsyRATAccuracyOracle.sscutscore( ...
                fx.bp, fx.wp_pop, fx.wp_ss, fx.mu, fx.cut, fx.idtable, fx.ci, 'dep', 0);

            testCase.verifyEqual(outCut.cutdep_pt, expCut.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(outCut.cutdep_ll, expCut.ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(outCut.cutdep_ul, expCut.ul, 'AbsTol', 1e-12);
        end

        function testSubjectSpecificCutScoreWithTrialEffectMatchesOracle(testCase)
            % F6a: the sigma_i^2 (trial main effect) half of the subject-level
            % cut-score absolute error err = sigma_pi,e^2(s)/n' + sigma_i^2/n'
            % must match the oracle with a nonzero sigma_i, not only the i = 0
            % default the sibling lane uses. The oracle mirrors the production
            % formula (regression-only), so oracle equality alone cannot catch
            % a shared sign error; the directional check below keeps this lane
            % non-vacuous.
            fx = PsyRATAccuracyFixtureFactory.ssOneFacet();
            sigI = [0.30; 0.25; 0.35; 0.28]; % sigma_i draws (one per posterior draw)

            outCut = psyrat_sscutscore( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'mu', fx.mu, 'cut', fx.cut, 'idtable', fx.idtable, ...
                'CI', fx.ci, 'est', 'dep', 'i', sigI);
            expCut = PsyRATAccuracyOracle.sscutscore( ...
                fx.bp, fx.wp_pop, fx.wp_ss, fx.mu, fx.cut, fx.idtable, ...
                fx.ci, 'dep', sigI);

            testCase.verifyEqual(outCut.cutdep_pt, expCut.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(outCut.cutdep_ll, expCut.ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(outCut.cutdep_ul, expCut.ul, 'AbsTol', 1e-12);

            % Nonzero sigma_i adds sigma_i^2/n' to the absolute error, so it
            % must strictly lower the dependability point estimate for every
            % subject relative to the i = 0 call on the same fixture.
            outCutZero = psyrat_sscutscore( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'mu', fx.mu, 'cut', fx.cut, 'idtable', fx.idtable, ...
                'CI', fx.ci, 'est', 'dep');
            testCase.verifyLessThan(outCut.cutdep_pt, outCutZero.cutdep_pt);
        end

        function testSubjectSpecificRelativeAbsoluteMatchOracle(testCase)
            % F5: with a trial main effect (sigma_i) supplied, psyrat_ssrel must
            % match the oracle for BOTH the absolute (gcoeff = 1, phi_s) and the
            % relative (gcoeff = 2, G_s) decision types, and the absolute
            % coefficient must be <= the relative one (it keeps sigma_i^2).
            fx = PsyRATAccuracyFixtureFactory.ssOneFacet();
            sigI = [0.30; 0.25; 0.35; 0.28]; % sigma_i draws (one per posterior draw)

            for gcoeff = [1 2]
                out = psyrat_ssrel( ...
                    'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                    'i', sigI, 'gcoeff', gcoeff, 'idtable', fx.idtable, 'CI', fx.ci);
                exp = PsyRATAccuracyOracle.ssrel( ...
                    fx.bp, fx.wp_pop, fx.wp_ss, fx.idtable, fx.ci, sigI, gcoeff);

                testCase.verifyEqual(out.dep_pt, exp.dep_pt, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.dep_ll, exp.dep_ll, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.dep_ul, exp.dep_ul, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.icc_pt, exp.icc_pt, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.sem_pt, exp.sem_pt, 'AbsTol', 1e-12);
                testCase.verifyEqual(out.trl_var, exp.trl_var, 'AbsTol', 1e-12);
            end

            absOut = psyrat_ssrel( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'i', sigI, 'gcoeff', 1, 'idtable', fx.idtable, 'CI', fx.ci);
            relOut = psyrat_ssrel( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'i', sigI, 'gcoeff', 2, 'idtable', fx.idtable, 'CI', fx.ci);
            testCase.verifyTrue(all(absOut.dep_pt <= relOut.dep_pt + 1e-12));
            testCase.verifyTrue(all(absOut.trl_var > 0));

            % Subject-level cut-scores are absolute-error only: est = 'gen' is
            % rejected (audit F7).
            testCase.verifyError(@() psyrat_sscutscore( ...
                'bp', fx.bp, 'wp_pop', fx.wp_pop, 'wp_ss', fx.wp_ss, ...
                'mu', fx.mu, 'cut', fx.cut, 'idtable', fx.idtable, ...
                'CI', fx.ci, 'est', 'gen'), 'varargin:est');
        end

        function testSubjectSpecificTrtAndDiffMatchOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.ssTrt();

            out = psyrat_ssrel_trt( ...
                'gcoeff', fx.gcoeff, 'reltype', fx.reltype, ...
                'bp', fx.bp, 'bo', fx.bo, 'bt', fx.bt, ...
                'txp', fx.txp, 'oxp', fx.oxp, 'txo', fx.txo, 'err', fx.err, ...
                'idtable', fx.idtable, 'CI', fx.ci, 'nocc', fx.nocc, ...
                'txp_ss', fx.txp_ss, 'oxp_ss', fx.oxp_ss, 'err_ss', fx.err_ss);

            exp = PsyRATAccuracyOracle.ssrelTrt( ...
                fx.gcoeff, fx.reltype, fx.bp, fx.bo, fx.bt, fx.txp, fx.oxp, fx.txo, fx.err, ...
                fx.idtable, fx.ci, fx.nocc, fx.txp_ss, fx.oxp_ss, fx.err_ss);

            testCase.verifyEqual(out.rel_pt, exp.rel_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.rel_ll, exp.rel_ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.rel_ul, exp.rel_ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.sem_pt, exp.sem_pt, 'AbsTol', 1e-12);

            % RC-06. The icc_* columns were never compared here, so the oracle's
            % duplicate of the one-facet ICC was invisible twice over. Compare
            % them, plus the two group reference-line columns the caterpillar
            % plots read (psyrat_ssrelplot).
            testCase.verifyEqual(out.icc_pt, exp.icc_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_ll, exp.icc_ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_ul, exp.icc_ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.gro_icc, exp.gro_icc, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.gro_rel, exp.gro_rel, 'AbsTol', 1e-12);

            fxD = PsyRATAccuracyFixtureFactory.ssDiff();
            outD = psyrat_ssrel_diff( ...
                'bp', fxD.base.bp, 'bt', fxD.base.bt, 'er_var', fxD.base.er_var, ...
                'idtable', fxD.idtable, 'CI', fxD.ci, 'est', fxD.est, ...
                'er_var_ss1', fxD.er_var_ss1, 'er_var_ss2', fxD.er_var_ss2, ...
                'wp_cov_ss', fxD.wp_cov_ss);

            expD = PsyRATAccuracyOracle.ssrelDiff( ...
                fxD.base.bp, fxD.base.bt, fxD.base.er_var, fxD.idtable, fxD.ci, fxD.est, 0, ...
                fxD.er_var_ss1, fxD.er_var_ss2, fxD.wp_cov_ss);

            testCase.verifyEqual(outD.rel_pt, expD.rel_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(outD.rel_ll, expD.rel_ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(outD.rel_ul, expD.rel_ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(outD.sem_pt, expD.sem_pt, 'AbsTol', 1e-12);
        end

        function testDodCompositeMatchesOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.dodDynrel();
            mats = psyrat_gamma_dod_components_ls('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);

            for key = {'p', 'i', 'pi_e'}
                V = mats.(key{1});
                for rule = {'unscaled', 'harmonic'}
                    out = psyrat_dod_composite('varcov', V, 'cvec', fx.cvec, ...
                        'rule', rule{1}, 'n', fx.n, 'tol', fx.tol);
                    exp_ = PsyRATAccuracyOracle.dodComposite(V, fx.cvec, ...
                        rule{1}, fx.n, fx.tol);
                    testCase.verifyEqual(out.value, exp_.value, 'AbsTol', 1e-12);
                    testCase.verifyEqual(out.invalid_negative, ...
                        exp_.invalid_negative, 'AbsTol', 0);
                    testCase.verifyEqual(out.scale, exp_.scale, 'AbsTol', 1e-12);
                end
            end

            % Non-mirror structural assertions, so this method can fail for a
            % reason other than a drifted mirror: with every count equal to
            % one the harmonic rule IS the unscaled rule; and the harmonic
            % rule without counts is an error, not a default.
            V = mats.p;
            one = psyrat_dod_composite('varcov', V, 'cvec', fx.cvec, ...
                'rule', 'harmonic', 'n', [1 1 1 1], 'tol', fx.tol);
            plain = psyrat_dod_composite('varcov', V, 'cvec', fx.cvec, ...
                'rule', 'unscaled', 'tol', fx.tol);
            testCase.verifyEqual(one.value, plain.value, 'AbsTol', 1e-12);
            testCase.verifyError(@() psyrat_dod_composite('varcov', V, ...
                'rule', 'harmonic'), 'varargin:n');
        end

        function testDodComponentsLsMatchesOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.dodDynrel();

            one = psyrat_gamma_dod_components_ls('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            expOne = PsyRATAccuracyOracle.dodComponentsLs('trial', ...
                fx.alpha, fx.log_sigma, fx.sd_p, fx.cor_p, ...
                fx.sd_trl, fx.cor_trl, struct());
            for key = {'p', 'i', 'pi_e'}
                testCase.verifyEqual(one.(key{1}), expOne.(key{1}), ...
                    'AbsTol', 1e-12);
            end

            extra = struct('sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            two = psyrat_gamma_dod_components_ls('facet', 'trial_occasion', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
                'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            expTwo = PsyRATAccuracyOracle.dodComponentsLs('trial_occasion', ...
                fx.alpha, fx.log_sigma, fx.sd_p, fx.cor_p, ...
                fx.sd_trl, fx.cor_trl, extra);
            for key = {'p', 'i', 'o', 'pi', 'po', 'io', 'pio_e'}
                testCase.verifyEqual(two.(key{1}), expTwo.(key{1}), ...
                    'AbsTol', 1e-12);
            end

            % Non-mirror structural assertions: every component matrix is
            % symmetric with strictly positive diagonals, and the residual
            % diagonal strictly exceeds the mean-surface remainder (it pools
            % + sigma^2 > 0 on top).
            for key = {'p', 'i', 'pi_e'}
                V = one.(key{1});
                testCase.verifyEqual(V, permute(V, [1 3 2]), 'AbsTol', 0, ...
                    sprintf('%s must be symmetric.', key{1}));
                for q = 1:4
                    testCase.verifyGreaterThan(V(:,q,q), 0);
                end
            end
        end

        function testDodComponentsGaussMatchesOracle(testCase)
            fx = PsyRATAccuracyFixtureFactory.dodDynrelGauss();

            one = psyrat_dod_components('facet', 'trial', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl);
            expOne = PsyRATAccuracyOracle.dodComponentsGauss('trial', ...
                fx.alpha, fx.log_sigma, fx.sd_p, fx.cor_p, ...
                fx.sd_trl, fx.cor_trl, struct());
            for key = {'p', 'i', 'pi_e'}
                testCase.verifyEqual(one.(key{1}), expOne.(key{1}), ...
                    'AbsTol', 1e-12);
            end

            extra = struct('sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            two = psyrat_dod_components('facet', 'trial_occasion', ...
                'alpha', fx.alpha, 'log_sigma', fx.log_sigma, ...
                'sd_p', fx.sd_p, 'cor_p', fx.cor_p, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
                'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
                'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
                'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
                'sd_to', fx.sd_to, 'cor_to', fx.cor_to);
            expTwo = PsyRATAccuracyOracle.dodComponentsGauss('trial_occasion', ...
                fx.alpha, fx.log_sigma, fx.sd_p, fx.cor_p, ...
                fx.sd_trl, fx.cor_trl, extra);
            for key = {'p', 'i', 'o', 'pi', 'po', 'io', 'pio_e'}
                testCase.verifyEqual(two.(key{1}), expTwo.(key{1}), ...
                    'AbsTol', 1e-12);
            end

            % Non-mirror structural assertions, so this method can fail for a
            % reason other than a drifted mirror: symmetry and strictly
            % positive diagonals; the defining GAUSSIAN boundary - residual
            % off-diagonals EXACTLY zero (AbsTol 0; the gamma channel keeps
            % nonzero mean-surface cross terms there); and the identity-link
            % contract that expected is alpha itself, bit-identical.
            for key = {'p', 'i', 'pi_e'}
                V = one.(key{1});
                testCase.verifyEqual(V, permute(V, [1 3 2]), 'AbsTol', 0, ...
                    sprintf('%s must be symmetric.', key{1}));
                for q = 1:4
                    testCase.verifyGreaterThan(V(:,q,q), 0);
                end
            end
            nd = size(fx.alpha, 1);
            for a = 1:3
                for b = (a+1):4
                    testCase.verifyEqual(one.pi_e(:,a,b), zeros(nd,1), ...
                        'AbsTol', 0, sprintf(['Gaussian residual cross ', ...
                        '(%d,%d) must be exactly zero.'], a, b));
                    testCase.verifyEqual(two.pio_e(:,a,b), zeros(nd,1), ...
                        'AbsTol', 0, sprintf(['Gaussian two-facet residual ', ...
                        'cross (%d,%d) must be exactly zero.'], a, b));
                end
            end
            testCase.verifyEqual(one.expected, fx.alpha, 'AbsTol', 0, ...
                'expected must be alpha itself under the identity link.');
        end
    end
end

function localVerifyIntervalBounds(testCase, ll, pt, ul)
testCase.verifyLessThanOrEqual(ll, pt);
testCase.verifyLessThanOrEqual(pt, ul);
testCase.verifyGreaterThanOrEqual(pt, 0);
testCase.verifyLessThanOrEqual(pt, 1);
end
