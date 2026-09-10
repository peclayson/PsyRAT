classdef TestNonconcurrentDodGaussReliability < matlab.unittest.TestCase
    %TESTNONCONCURRENTDODGAUSSRELIABILITY Estimand tests for the GAUSSIAN arm
    % of the person-specific dynamic NONCONCURRENT difference-of-differences
    % (analysis 28, ic_dodiff_dynrel_sserrvar, family = 'gaussian'), plus the
    % summary dispatch and disclosure surfaces the arm shares with the gamma
    % family. The two-facet twin (analysis 29) is TestNonconcurrentDodGaussTrt-
    % Reliability.
    %
    % NO FROZEN FIXTURE LANE EXISTS FOR THIS ARM, and none can be minted: the
    % reference bundle is gamma end to end (SCIENTIFIC_FORMULA_AUDIT.md
    % section 26, owner ruling 3). The replacement evidence ladder here:
    %  - CONSTANT-DRAW CLOSED FORMS transcribed longhand (scalar sums, a
    %    different route from production's matrix quadratic forms through the
    %    shared composite kernel). Under Gaussian these are fully analytic -
    %    no lognormal moment identities intervene - so this two-route check is
    %    STRONGER than its gamma sibling, not weaker.
    %  - THE VETTED-KERNEL TIE, upgraded from the gamma increment's shared
    %    subspace to the FULL input space: the shipped static Gaussian DoD
    %    kernel PSYRAT_DODIFFREL is exactly this estimand at the typical
    %    person, because the Gaussian residual channel is diagonal by theorem.
    %    The tie is a TEST, not a call; the static kernel stays untouched.
    % Live CmdStan planted-truth recovery (tests/integration) is the
    % load-bearing end-to-end gate, per the same ruling.
    %
    % Tolerances: AbsTol 1e-10 for longhand closed-form transcriptions (the
    % analysis-13 precedent); AbsTol 1e-12 for same-arithmetic routes; AbsTol
    % 0 only for bit-identity claims with a stated reason.

    methods (Test)

        function testPerParticipantMatchesClosedForm(testCase)
            %Constant draws: every per-participant coefficient is an explicit
            %scalar function of the planted numbers, transcribed longhand.
            fx = localGaussConstFixture();
            ssrel = localRunGaussSsrel(fx);
            testCase.assertEqual(height(ssrel), fx.P);

            w = fx.cvec(:);
            Sp = localCovFromBlocks(fx.sd_id(1,1:4), squeeze(fx.cor_id(1,1:4,1:4)));
            Si = localCovFromBlocks(fx.sd_trl(1,:), squeeze(fx.cor_trl(1,:,:)));
            u = w' * Sp * w;

            for c = 1:fx.P
                z = fx.idtable.z1(c);
                n_act = [fx.idtable.trls1(c) fx.idtable.trls2(c) ...
                    fx.idtable.trls3(c) fx.idtable.trls4(c)];
                sig = zeros(1,4); mu_own = zeros(1,4);
                for q = 1:4
                    sig(q) = exp(fx.er_ss.(sprintf('c%d',q))(1,c) + ...
                        fx.b_sigma_dim(1,q,1) * z);
                    mu_own(q) = fx.mu_ss.(sprintf('c%d',q))(1,c) + ...
                        fx.b_dim(1,q,1) * z;
                end

                relerr = sum((w'.^2) .* sig.^2 ./ n_act);
                relerr_unscaled = sum((w'.^2) .* sig.^2);
                i_scaled = localHarmonicLonghand(Si, w, n_act);
                i_unscaled = w' * Si * w;
                abserr = relerr + i_scaled;

                gen = u / (u + relerr);
                dep = u / (u + abserr);
                iccg = u / (u + relerr_unscaled);
                iccd = u / (u + relerr_unscaled + i_unscaled);

                testCase.verifyEqual(ssrel.gen_pt(c), gen, 'AbsTol', 1e-10, ...
                    sprintf('gen_pt row %d', c));
                testCase.verifyEqual(ssrel.dep_pt(c), dep, 'AbsTol', 1e-10, ...
                    sprintf('dep_pt row %d', c));
                testCase.verifyEqual(ssrel.iccg_pt(c), iccg, 'AbsTol', 1e-10, ...
                    sprintf('iccg_pt row %d', c));
                testCase.verifyEqual(ssrel.iccd_pt(c), iccd, 'AbsTol', 1e-10, ...
                    sprintf('iccd_pt row %d', c));
                testCase.verifyEqual(ssrel.uni_pt(c), u, 'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.relerr_pt(c), relerr, ...
                    'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.abserr_pt(c), abserr, ...
                    'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.sem_gen_pt(c), sqrt(relerr), ...
                    'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.sem_dep_pt(c), sqrt(abserr), ...
                    'AbsTol', 1e-10);

                %identity-link per-person quantities: the expected score is
                %the participant's own additive predictor (no lognormal mean
                %correction), and the DoD residual variance is EXACTLY the
                %four-term sum
                testCase.verifyEqual(ssrel.expected_dod_person_pt(c), ...
                    mu_own * w, 'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.residual_var_dod_person_pt(c), ...
                    sum(sig.^2), 'AbsTol', 1e-10);
                for q = 1:4
                    testCase.verifyEqual( ...
                        ssrel.(sprintf('ressd%d_pt',q))(c), sig(q), ...
                        'AbsTol', 1e-10);
                end

                %constant draws collapse every credible interval onto the
                %point estimate
                testCase.verifyEqual(ssrel.gen_ll(c), ssrel.gen_pt(c), ...
                    'AbsTol', 1e-12);
                testCase.verifyEqual(ssrel.gen_ul(c), ssrel.gen_pt(c), ...
                    'AbsTol', 1e-12);
                testCase.verifyEqual(ssrel.invalid_frac(c), 0);
                testCase.verifyEqual(ssrel.n_floored(c), 0);
            end
        end

        function testStandardizedTwinsDiffer(testCase)
            %The standardized design must be genuinely computed at nprime, not
            %copied from the actual design.
            fx = localGaussConstFixture();
            ssrel = localRunGaussSsrel(fx);
            w = fx.cvec(:);
            Sp = localCovFromBlocks(fx.sd_id(1,1:4), squeeze(fx.cor_id(1,1:4,1:4)));
            Si = localCovFromBlocks(fx.sd_trl(1,:), squeeze(fx.cor_trl(1,:,:)));
            u = w' * Sp * w;
            for c = 1:fx.P
                z = fx.idtable.z1(c);
                sig = zeros(1,4);
                for q = 1:4
                    sig(q) = exp(fx.er_ss.(sprintf('c%d',q))(1,c) + ...
                        fx.b_sigma_dim(1,q,1) * z);
                end
                rel_s = sum((w'.^2) .* sig.^2 ./ fx.nprime);
                gen_s = u / (u + rel_s);
                testCase.verifyEqual(ssrel.gen_pt_std(c), gen_s, ...
                    'AbsTol', 1e-10, sprintf('gen_pt_std row %d', c));
                testCase.verifyNotEqual(ssrel.gen_pt_std(c), ssrel.gen_pt(c));
                testCase.verifyEqual(ssrel.relerr_pt_std(c), rel_s, ...
                    'AbsTol', 1e-10);
            end
        end

        function testTableSchemaSharedWithGammaMinusNu(testCase)
            %The shared consumers (viewer table, flag counter, plot overlay)
            %read the same columns from both families' tables. The one
            %deliberate difference: nu = 2*(mu/sigma)^2 is a scaled-chi-square
            %quantity, so the Gaussian table carries NO nu columns; ressd
            %columns stay.
            fx = localGaussConstFixture();
            ssrel = localRunGaussSsrel(fx);
            required = {'id','id2','z1','trls1','trls2','trls3','trls4', ...
                'gen_pt','gen_ll','gen_ul','dep_pt','dep_ll','dep_ul', ...
                'iccg_pt','iccg_ll','iccg_ul','iccd_pt','iccd_ll','iccd_ul', ...
                'sem_gen_pt','sem_dep_pt', ...
                'gen_pt_std','dep_pt_std','sem_gen_pt_std','sem_dep_pt_std', ...
                'uni_pt','relerr_pt','abserr_pt','relerr_pt_std', ...
                'abserr_pt_std','expected_dod_pop_pt', ...
                'expected_dod_person_pt','residual_var_dod_person_pt', ...
                'ressd1_pt','ressd2_pt','ressd3_pt','ressd4_pt', ...
                'invalid_frac','invalid_frac_std','n_floored'};
            for k = 1:numel(required)
                testCase.verifyTrue( ...
                    ismember(required{k}, ssrel.Properties.VariableNames), ...
                    sprintf('missing required column %s', required{k}));
            end
            for q = 1:4
                testCase.verifyFalse(ismember(sprintf('nu%d_pt',q), ...
                    ssrel.Properties.VariableNames), ...
                    'nu columns are scaled-chi-square quantities; drop them.');
            end
        end

        function testSurfaceMatchesZeroEffectParticipant(testCase)
            %The typical-person surface at z must equal the per-participant
            %read-out of a participant whose person effects are all zero and
            %whose counts equal the surface n'. Nonzero dimension slopes make
            %this non-vacuous in z. AbsTol 1e-12: identical draws through two
            %code paths.
            fx = localGaussRandomFixture(60);
            z = 0.7;
            obs = [10 12 9 15];

            rel = psyrat_rel_dodiffdynrel('b', fx.b, 'b_sigma', fx.b_sigma, ...
                'b_dim', fx.b_dim, 'b_sigma_dim', fx.b_sigma_dim, ...
                'sd_id', fx.sd_id, 'cor_id', fx.cor_id, ...
                'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
                'ndim', 1, 'z1', z, 'obs', obs, 'CI', 0.95);

            idt = table({'zero'}, 1, z, obs(1), obs(2), obs(3), obs(4), ...
                'VariableNames', {'id','id2','z1','trls1','trls2','trls3','trls4'});
            ssrel = psyrat_ssrel_dodiffdynrel('b', fx.b, ...
                'b_sigma', fx.b_sigma, 'b_dim', fx.b_dim, ...
                'b_sigma_dim', fx.b_sigma_dim, 'sd_id', fx.sd_id, ...
                'cor_id', fx.cor_id, 'sd_trl', fx.sd_trl, ...
                'cor_trl', fx.cor_trl, ...
                'mu_ss1', fx.b(:,1), 'mu_ss2', fx.b(:,2), ...
                'mu_ss3', fx.b(:,3), 'mu_ss4', fx.b(:,4), ...
                'er_ss1', fx.b_sigma(:,1), 'er_ss2', fx.b_sigma(:,2), ...
                'er_ss3', fx.b_sigma(:,3), 'er_ss4', fx.b_sigma(:,4), ...
                'idtable', idt, 'nprime', obs, 'ndim', 1, 'CI', 0.95);

            testCase.verifyEqual(rel.G.pt, ssrel.gen_pt, 'AbsTol', 1e-12, ...
                'Surface G must equal the zero-effect participant.');
            testCase.verifyEqual(rel.D.pt, ssrel.dep_pt, 'AbsTol', 1e-12, ...
                'Surface D must equal the zero-effect participant.');
            testCase.verifyEqual(rel.ICCg.pt, ssrel.iccg_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(rel.ICCd.pt, ssrel.iccd_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(rel.estimand, ...
                'typical person (all person effects at zero)');
        end

        function testVettedKernelTieFullSpace(testCase)
            %FULL-SPACE tie to the shipped static Gaussian DoD kernel: on
            %arbitrary seeded draws, the Gaussian stage-3 draws must equal
            %PSYRAT_DODIFFREL draw for draw. Full-space (not the gamma tie's
            %shared subspace) because the Gaussian residual channel is
            %diagonal by theorem - the static kernel's independent-residual
            %arithmetic IS this estimand at the typical person. The tie is a
            %TEST, not a call; the static kernel stays untouched. AbsTol
            %1e-12: same arithmetic assembled in a different order.
            rng(24601, 'twister');   %the gamma tie's own fixture doctrine
            nd = 500;
            obs = [18 22 15 20];
            Lp = [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22];
            Lt = [0.28 0.12; 0.24 -0.15; -0.20 0.14; 0.18 0.11];
            base_bp = Lp * Lp' + diag([0.35 0.42 0.30 0.38]);
            base_bt = Lt * Lt' + diag([0.06 0.075 0.055 0.07]);
            bp = zeros(nd, 4, 4); bt = zeros(nd, 4, 4);
            for d = 1:nd
                s = 1 + 0.3 * rand;
                bp(d, :, :) = base_bp * s;
                bt(d, :, :) = base_bt * (2 - s);
            end
            er = [log(1.5) log(1.2) log(1.8) log(1.35)] + 0.2 * randn(nd, 4);

            dep = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'dep', 'ci', .95);
            gen = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'gen', 'ci', .95);

            %decompose the target covariance draws into the calculator's
            %sd/correlation inputs
            [sd_p, cor_p] = localDecompose(bp);
            [sd_t, cor_t] = localDecompose(bt);
            sd_id = [sd_p, 0.3 * ones(nd, 4)];
            cor_id = zeros(nd, 8, 8);
            cor_id(:, 1:4, 1:4) = cor_p;
            for j = 5:8
                cor_id(:, j, j) = 1;
            end

            zeroslope = zeros(nd, 4, 1);
            idt = table({'tie'}, 1, 0, obs(1), obs(2), obs(3), obs(4), ...
                'VariableNames', {'id','id2','z1','trls1','trls2','trls3','trls4'});
            [~, dl] = psyrat_ssrel_dodiffdynrel('b', zeros(nd,4), ...
                'b_sigma', er, 'b_dim', zeroslope, 'b_sigma_dim', zeroslope, ...
                'sd_id', sd_id, 'cor_id', cor_id, 'sd_trl', sd_t, ...
                'cor_trl', cor_t, ...
                'mu_ss1', zeros(nd,1), 'mu_ss2', zeros(nd,1), ...
                'mu_ss3', zeros(nd,1), 'mu_ss4', zeros(nd,1), ...
                'er_ss1', er(:,1), 'er_ss2', er(:,2), ...
                'er_ss3', er(:,3), 'er_ss4', er(:,4), ...
                'idtable', idt, 'nprime', obs, 'ndim', 1, 'CI', 0.95);

            testCase.verifyEqual(dl.G(:,1), gen.rel_draws, 'AbsTol', 1e-12, ...
                'Generalizability draws must match the vetted static kernel.');
            testCase.verifyEqual(dl.D(:,1), dep.rel_draws, 'AbsTol', 1e-12, ...
                'Dependability draws must match the vetted static kernel.');
            testCase.verifyEqual(dl.iccg(:,1), gen.icc_draws, 'AbsTol', 1e-12, ...
                'Generalizability ICC draws must match the vetted static kernel.');
            testCase.verifyEqual(dl.iccd(:,1), dep.icc_draws, 'AbsTol', 1e-12, ...
                'Dependability ICC draws must match the vetted static kernel.');
        end

        function testDepNotAboveGenOnUnflaggedRows(testCase)
            %The absolute error adds a channel, so D <= G wherever nothing is
            %flagged (a materially negative trial composite can invert it,
            %and is flagged when it does).
            fx = localGaussRandomFixture(80);
            P = 4;
            idt = table(compose('s%d', (1:P)'), (1:P)', ...
                [-0.9; -0.1; 0.4; 1.1], ...
                [12; 8; 20; 10], [9; 15; 6; 11], [14; 10; 18; 9], ...
                [11; 13; 7; 16], 'VariableNames', ...
                {'id','id2','z1','trls1','trls2','trls3','trls4'});
            %person effects: one column per PARTICIPANT (P x 4, transposed
            %into the [ndraws x nsub] per-subject matrices)
            ulocP = 0.5 * [0.9 -0.5 0.3 -0.7; -1.1 0.6 -0.2 0.8; ...
                0.2 0.1 -0.4 0.5; -0.6 0.3 0.7 -0.2];
            ulsigP = 0.1 * [1.5 -1.0 0.5 -2.0; -0.8 1.2 -0.4 0.6; ...
                0.3 -0.2 0.9 -1.1; -1.4 0.7 -0.6 1.0];
            args = {'b', fx.b, 'b_sigma', fx.b_sigma, 'b_dim', fx.b_dim, ...
                'b_sigma_dim', fx.b_sigma_dim, 'sd_id', fx.sd_id, ...
                'cor_id', fx.cor_id, 'sd_trl', fx.sd_trl, ...
                'cor_trl', fx.cor_trl, 'idtable', idt, ...
                'nprime', [25 25 25 25], 'ndim', 1, 'CI', 0.95};
            for q = 1:4
                args = [args, {sprintf('mu_ss%d',q), ...
                    fx.b(:,q) + ulocP(:,q)', ...
                    sprintf('er_ss%d',q), ...
                    fx.b_sigma(:,q) + ulsigP(:,q)'}]; %#ok<AGROW>
            end
            ssrel = psyrat_ssrel_dodiffdynrel(args{:});
            ok = ssrel.invalid_frac == 0;
            testCase.assertTrue(any(ok), 'fixture produced no unflagged rows');
            testCase.verifyTrue(all(ssrel.dep_pt(ok) <= ssrel.gen_pt(ok) + 1e-9), ...
                'dependability exceeded generalizability on unflagged rows');
            testCase.verifyTrue(all(ssrel.gen_pt(ok) > 0 & ssrel.gen_pt(ok) < 1));
        end

        function testCouplingVocabularyRejected(testCase)
            %The nonconcurrent boundary at all four Gaussian calculator
            %signatures: coupling vocabulary errors by name, never silently
            %swallowed.
            calcs = { ...
                @psyrat_rel_dodiffdynrel, 'psyrat_rel_dodiffdynrel'; ...
                @psyrat_ssrel_dodiffdynrel, 'psyrat_ssrel_dodiffdynrel'; ...
                @psyrat_rel_dodiffdynrel_trt, 'psyrat_rel_dodiffdynrel_trt'; ...
                @psyrat_ssrel_dodiffdynrel_trt, 'psyrat_ssrel_dodiffdynrel_trt'};
            for k = 1:size(calcs,1)
                fn = calcs{k,1}; name = calcs{k,2};
                for coupled = {'rho', 'rho_e', 'cov_e_obs', ...
                        'cov_res_e_obs', 'er_cov', 'rescor', 'copula'}
                    testCase.verifyError(@() fn(coupled{1}, 0.3), ...
                        sprintf('%s:noresidualcoupling', name), ...
                        sprintf('%s must reject %s', name, coupled{1}));
                end
            end
        end

        function testZeroCountCellIsFlooredNotFatal(testCase)
            %B16 port: a participant with zero trials in one cell floors that
            %count to 1 for the composite, keeps the true zero in the table,
            %and discloses through n_floored. Corrupt counts still error.
            fx = localGaussConstFixture();
            fx.idtable.trls2(2) = 0;
            ssrel = localRunGaussSsrel(fx);
            testCase.verifyEqual(ssrel.trls2(2), 0, ...
                'the table must keep the true zero');
            testCase.verifyEqual(ssrel.n_floored(2), 1);
            testCase.verifyEqual(ssrel.n_floored([1 3]), [0; 0]);
            testCase.verifyTrue(isfinite(ssrel.gen_pt(2)), ...
                'the floored row must still carry coefficients');

            %longhand at the floored count: cell 2's residual term at n = 1
            w = fx.cvec(:);
            Sp = localCovFromBlocks(fx.sd_id(1,1:4), squeeze(fx.cor_id(1,1:4,1:4)));
            u = w' * Sp * w;
            z = fx.idtable.z1(2);
            n_eff = [fx.idtable.trls1(2) 1 fx.idtable.trls3(2) fx.idtable.trls4(2)];
            sig = zeros(1,4);
            for q = 1:4
                sig(q) = exp(fx.er_ss.(sprintf('c%d',q))(1,2) + ...
                    fx.b_sigma_dim(1,q,1) * z);
            end
            relerr = sum((w'.^2) .* sig.^2 ./ n_eff);
            testCase.verifyEqual(ssrel.gen_pt(2), u / (u + relerr), ...
                'AbsTol', 1e-10, 'the floored coefficient must use n = 1');

            %NaN and negative counts must keep erroring loudly (never floored)
            fxbad = localGaussConstFixture();
            fxbad.idtable.trls1(1) = NaN;
            testCase.verifyError(@() localRunGaussSsrel(fxbad), 'varargin:n');
            fxneg = localGaussConstFixture();
            fxneg.idtable.trls1(1) = -3;
            testCase.verifyError(@() localRunGaussSsrel(fxneg), 'varargin:n');
        end

        function testSummaryDispatchGaussian(testCase)
            %The dodgammals flip: a Gaussian-stamped DoD REL routes through
            %psyrat_dynrel_summary to the Gaussian calculators; the fail-loud
            %arms keep catching every other stamp (log-nu gamma, and a REL
            %with no family stamp at all - the damaged-file arrival the
            %documented case-13 fall-through failure mode is about).
            fx = localGaussConstFixture();
            REL = localSyntheticGaussRel(fx);
            w = warning('off', 'dynrel:dodinvalidity');
            cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
            summ = psyrat_dynrel_summary(REL);
            st = summ.strata(1);
            testCase.assertTrue(istable(st.ssrel_table));
            testCase.verifyEqual(height(st.ssrel_table), 3);
            testCase.verifyEqual(numel(st.obs_ev), 4);
            %dispatch pin: the summary's table must equal a DIRECT calculator
            %call on the same draws at the summary's resolved standardized
            %counts (the stored per-cell medians). AbsTol 1e-12: identical
            %draws through two entry points - this is what catches a
            %slot-order or slot-semantics mistake in the shared dargs/ssargs
            %assembly.
            fxd = fx;
            fxd.nprime = max(1, round(REL.out.ntrials(:,1)'));
            direct = localRunGaussSsrel(fxd);
            testCase.verifyEqual(st.ssrel_table.gen_pt, direct.gen_pt, ...
                'AbsTol', 1e-12, 'summary vs direct gen_pt diverged');
            testCase.verifyEqual(st.ssrel_table.dep_pt, direct.dep_pt, ...
                'AbsTol', 1e-12, 'summary vs direct dep_pt diverged');
            testCase.verifyEqual(st.ssrel_table.gen_pt_std, ...
                direct.gen_pt_std, 'AbsTol', 1e-12, ...
                'summary vs direct standardized twin diverged');
            testCase.verifyTrue(st.dod_flags.recorded);
            testCase.verifyEmpty(st.copula_flags, ...
                'No copula machinery may touch a nonconcurrent DoD result.');
            testCase.verifyEqual(summ.dod_labels.inference_framework, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1');
            testCase.verifyFalse(summ.isgamma);

            %fail-loud arms
            RELlognu = REL;
            RELlognu.family = 'gamma';
            RELlognu.gammascale = 1;
            testCase.verifyError(@() psyrat_dynrel_summary(RELlognu), ...
                'dynrel:dodgammals');
            RELnofam = rmfield(REL, 'family');
            testCase.verifyError(@() psyrat_dynrel_summary(RELnofam), ...
                'dynrel:dodgammals');
        end

        function testProvenanceDisclosureGaussian(testCase)
            %The nonconcurrent-DoD disclosure block must reach a GAUSSIAN run
            %(it was born inside a family=='gamma' gate), with the Gaussian
            %framework token and WITHOUT the gamma-only lines.
            REL = localSyntheticGaussRel(localGaussConstFixture());
            pd = localReportData(REL);
            lines = psyrat_provenance_lines(pd);   %the wrapper; the function unwraps .rel
            txt = strjoin(lines(:)', ' | ');
            testCase.verifyTrue(contains(txt, ...
                'NONCONCURRENT COUNTERFACTUAL ESTIMAND'), ...
                'the counterfactual disclosure must reach the Gaussian arm');
            testCase.verifyTrue(contains(txt, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1'), ...
                'the Gaussian framework token must be quoted');
            testCase.verifyTrue(contains(txt, 'Gaussian margins'), ...
                'the framework sentence must name the family');
            testCase.verifyFalse(contains(txt, 'modular_four_gamma_margins'), ...
                'a Gaussian run must not carry the gamma framework token');
            testCase.verifyFalse(contains(txt, 'Gamma scale parameterization'), ...
                'the gammascale contrast line is gamma-only');
            testCase.verifyFalse(contains(txt, 'log-nu'), ...
                'the log-nu contrast is gamma-only');
            testCase.verifyTrue(contains(txt, 'record-never-repair'), ...
                'the signed-projection validity line must be present');

            %degraded arm: the string-keyed fallback still fires without the
            %carrier, and the framework token is gone
            pdbare = pd;
            pdbare.rel.out = rmfield(pdbare.rel.out, 'dod_labels');
            linesbare = psyrat_provenance_lines(pdbare);
            txtbare = strjoin(linesbare(:)', ' | ');
            testCase.verifyTrue(contains(txtbare, 'labels missing'), ...
                'the degraded fallback must disclose the missing carrier');
            testCase.verifyFalse(contains(txtbare, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1'));

            %negative control: a non-DoD analysis emits no DoD block
            pdneg = pd;
            pdneg.rel.analysis = 'ic';
            pdneg.rel.out = rmfield(pdneg.rel.out, 'dod_labels');
            linesneg = psyrat_provenance_lines(pdneg);
            txtneg = strjoin(linesneg(:)', ' | ');
            testCase.verifyFalse(contains(txtneg, 'NONCONCURRENT'), ...
                'the DoD block leaked onto a non-DoD analysis');
        end

        function testReportEstimandNoteGaussian(testCase)
            %B15 port for the Gaussian arm: the headless report's estimand
            %note must carry the typical-person sentence on the populated
            %path, decline to assert an estimand on the degraded path, and
            %never claim the marginal wording.
            REL = localSyntheticGaussRel(localGaussConstFixture());
            pd = localReportData(REL);
            w = warning('off', 'dynrel:dodinvalidity');
            cleaner = onCleanup(@() warning(w)); %#ok<NASGU>

            rep = psyrat_report(pd);
            testCase.assertTrue(contains(rep.estimand_note, 'typical person'), ...
                'the typical-person sentence must reach the Gaussian arm');
            testCase.verifyFalse(contains(rep.estimand_note, ...
                'population average'));

            pdbare = pd;
            pdbare.rel.out = rmfield(pdbare.rel.out, 'dod_labels');
            repbare = psyrat_report(pdbare);
            testCase.verifyFalse(contains(repbare.estimand_note, ...
                'population average'));
            testCase.verifyFalse(contains(repbare.estimand_note, ...
                'typical person'));
            testCase.verifyTrue(contains(repbare.estimand_note, 'dod_labels'), ...
                'the degraded note must name the missing carrier');
        end

        function testViewerHeaderEstimandGaussian(testCase)
            %B15, viewer surface, Gaussian arm: the viewer duplicates the
            %report's estimand branching, so it needs its own pin. Drive it
            %headlessly and read the written variance-table CSV header.
            REL = localSyntheticGaussRel(localGaussConstFixture());
            pd = localReportData(REL);
            w = warning('off', 'dynrel:dodinvalidity');
            cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
            shadowCleaner = localInstallSaveShadow(); %#ok<NASGU>

            fig = localDodViewerFig(testCase, pd);
            figCleaner = onCleanup(@() close(fig, 'force'));
            csvfile = [tempname '.csv'];
            localFireDodSave(testCase, fig, csvfile);
            txt = fileread(csvfile);
            delete(csvfile);
            testCase.assertTrue(contains(txt, 'typical person'), ...
                'viewer header lost the typical-person line on the Gaussian arm');
            testCase.verifyFalse(contains(txt, 'population average'));
            testCase.verifyFalse(contains(txt, 'modular_four_gamma_margins'), ...
                'a Gaussian export must not carry the gamma framework token');

            %the ON-SCREEN variance-window disclosure note must name the
            %right family too (the adversarial pre-merge review caught it
            %hardcoding 'Gamma margins' on the Gaussian arm)
            btns = findall(fig, 'Style', 'pushbutton');
            fired = false;
            for bb = btns(:)'
                s = bb.String; if iscell(s), s = strjoin(s, ' '); end
                if contains(string(s), 'Show variance')
                    cb = bb.Callback; cb(bb, []); fired = true; break;
                end
            end
            testCase.assertTrue(fired, 'Show variance components button not found');
            vf = findall(0, 'Type', 'figure', 'Name', ...
                'Dynamic reliability: variance components');
            testCase.assertNotEmpty(vf, 'variance-components window not created');
            vfCleaner = onCleanup(@() close(vf(1), 'force')); %#ok<NASGU>
            notes = findall(vf(1), 'Style', 'text');
            notetxt = '';
            for nn = notes(:)'
                s = nn.String; if iscell(s), s = strjoin(s, ' '); end
                notetxt = [notetxt ' ' char(string(s))]; %#ok<AGROW>
            end
            testCase.verifyTrue(contains(notetxt, 'Gaussian margins'), ...
                'the on-screen DoD note must name the Gaussian family');
            testCase.verifyFalse(contains(notetxt, 'Gamma margins'), ...
                'the on-screen DoD note must not claim Gamma margins');
            clear vfCleaner
            clear figCleaner
        end

    end
end

% ---------------------------------------------------------------------------
% Fixtures. Constant-draw fixture: literal values, nd identical rows, so every
% posterior summary collapses onto the closed form. Random fixture: seeded
% draws for route-identity and invariant tests. Gaussian (identity-link)
% scales throughout, negative cell means deliberate.
% ---------------------------------------------------------------------------

function fx = localGaussConstFixture()
nd = 4;
fx.nd = nd;
fx.P = 3;
fx.cvec = [1 -1 -1 1];
fx.ndim = 1;
fx.b = repmat([-4.8 -1.9 -4.1 -3.4], nd, 1);
fx.b_sigma = repmat([0.61 1.01 0.64 0.66], nd, 1);
fx.b_dim = repmat(reshape([0.6 -0.25 0.4 -0.1], [1 4 1]), nd, 1);
fx.b_sigma_dim = repmat(reshape([-0.08 0.12 -0.05 0.07], [1 4 1]), nd, 1);
fx.sd_id = repmat([2.0 1.7 2.2 1.9 0.30 0.34 0.28 0.32], nd, 1);
fx.cor_id = repmat(reshape(localCorr8(), [1 8 8]), nd, 1);
fx.sd_trl = repmat([0.80 0.70 0.90 0.75], nd, 1);
Rt = [1 .25 .17 -.14; .25 1 -.19 .21; .17 -.19 1 .16; -.14 .21 .16 1];
fx.cor_trl = repmat(reshape(Rt, [1 4 4]), nd, 1);

uloc = [0.9 -0.5 0.3 -0.7; -1.1 0.6 -0.2 0.8; 0.2 0.1 -0.4 0.5];
ulsig = [0.15 -0.10 0.05 -0.20; -0.12 0.18 -0.06 0.09; 0.05 -0.03 0.11 -0.14];
for q = 1:4
    fx.mu_ss.(sprintf('c%d',q)) = repmat(fx.b(1,q) + uloc(:,q)', nd, 1);
    fx.er_ss.(sprintf('c%d',q)) = repmat(fx.b_sigma(1,q) + ulsig(:,q)', nd, 1);
end

fx.idtable = table({'s1';'s2';'s3'}, (1:3)', [-0.8; 0.1; 1.2], ...
    [12; 8; 20], [9; 15; 6], [14; 10; 18], [11; 13; 7], ...
    'VariableNames', {'id','id2','z1','trls1','trls2','trls3','trls4'});
fx.nprime = [30 30 30 30];
end

function fx = localGaussRandomFixture(nd)
rng(6363, 'twister');
fx.nd = nd;
fx.b = [-4.8 -1.9 -4.1 -3.4] + 0.3 * randn(nd, 4);
fx.b_sigma = [0.61 1.01 0.64 0.66] + 0.05 * randn(nd, 4);
fx.b_dim = reshape([0.6 -0.25 0.4 -0.1] + 0.05 * randn(nd, 4), [nd 4 1]);
fx.b_sigma_dim = reshape([-0.08 0.12 -0.05 0.07] + 0.02 * randn(nd, 4), ...
    [nd 4 1]);
fx.sd_id = abs([2.0 1.7 2.2 1.9 0.30 0.34 0.28 0.32] + ...
    [0.15*randn(nd,4) 0.03*randn(nd,4)]);
R8 = localCorr8();
fx.cor_id = repmat(reshape(R8, [1 8 8]), nd, 1);
fx.sd_trl = abs([0.80 0.70 0.90 0.75] + 0.06 * randn(nd, 4));
Rt = [1 .25 .17 -.14; .25 1 -.19 .21; .17 -.19 1 .16; -.14 .21 .16 1];
fx.cor_trl = repmat(reshape(Rt, [1 4 4]), nd, 1);
end

function R = localCorr8()
%A literal, valid 8x8 correlation: loadings + diagonal, normalized. Distinct
%off-diagonals with mixed signs; the [1:4,1:4] block is the person location
%correlation the calculators consume.
L = [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22; ...
    0.25 0.31; -0.28 0.24; 0.22 -0.26; 0.30 0.20];
S = L * L' + diag([0.60 0.68 0.75 0.83 0.71 0.77 0.82 0.74]);
R = diag(1 ./ sqrt(diag(S))) * S * diag(1 ./ sqrt(diag(S)));
end

function S = localCovFromBlocks(sd, R)
S = diag(sd) * R * diag(sd);
end

function v = localHarmonicLonghand(S, w, n)
%The harmonic count rule written longhand: diagonal terms carry their own
%cell's count; each covariance pair carries the pairwise harmonic mean.
v = 0;
for q = 1:4
    v = v + w(q)^2 * S(q,q) / n(q);
end
for a = 1:3
    for b = (a+1):4
        h = 2 / (1/n(a) + 1/n(b));
        v = v + 2 * w(a) * w(b) * S(a,b) / h;
    end
end
end

function [sd, cor] = localDecompose(V)
%Split covariance draws into SD and correlation draws (the calculator's input
%contract).
nd = size(V,1);
sd = zeros(nd,4);
cor = zeros(nd,4,4);
for d = 1:nd
    S = squeeze(V(d,:,:));
    s = sqrt(diag(S));
    sd(d,:) = s';
    cor(d,:,:) = diag(1./s) * S * diag(1./s);
end
end

function ssrel = localRunGaussSsrel(fx)
ssrel = psyrat_ssrel_dodiffdynrel('b', fx.b, 'b_sigma', fx.b_sigma, ...
    'b_dim', fx.b_dim, 'b_sigma_dim', fx.b_sigma_dim, ...
    'sd_id', fx.sd_id, 'cor_id', fx.cor_id, ...
    'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
    'mu_ss1', fx.mu_ss.c1, 'mu_ss2', fx.mu_ss.c2, ...
    'mu_ss3', fx.mu_ss.c3, 'mu_ss4', fx.mu_ss.c4, ...
    'er_ss1', fx.er_ss.c1, 'er_ss2', fx.er_ss.c2, ...
    'er_ss3', fx.er_ss.c3, 'er_ss4', fx.er_ss.c4, ...
    'idtable', fx.idtable, 'nprime', fx.nprime, ...
    'ndim', fx.ndim, 'CI', 0.95, 'cvec', fx.cvec);
end

function REL = localSyntheticGaussRel(fx)
%A synthetic Gaussian analysis-28 REL in the engine's storage idiom
%({num2cell(...)} per stratum), sufficient for the summary branch. The
%gammascale stamp is 1 deliberately: the engine stamps it on EVERY run and
%the Gaussian default resolution leaves it there - the summary must gate on
%family first (the documented trap its isgamma_ls comment records).
REL = struct();
REL.analysis = 'ic_dodiff_dynrel_sserrvar';
REL.family = 'gaussian';
REL.gammascale = 1;
REL.ndim = 1;
REL.dim_names = {'dim1'};
out = struct();
out.labels = {'none'};
out.b = {num2cell(fx.b)};
out.b_sigma = {num2cell(fx.b_sigma)};
out.b_dim = {num2cell(fx.b_dim)};
out.b_sigma_dim = {num2cell(fx.b_sigma_dim)};
out.sd_id = {num2cell(fx.sd_id)};
out.sd_trl = {num2cell(fx.sd_trl)};
out.cor_id = {num2cell(fx.cor_id)};
out.cor_trl = {num2cell(fx.cor_trl)};
for q = 1:4
    out.(sprintf('mu_ss%d',q)) = {num2cell(fx.mu_ss.(sprintf('c%d',q)))};
    out.(sprintf('er_var_ss%d',q)) = {num2cell(fx.er_ss.(sprintf('c%d',q)))};
end
out.ssinfo = {fx.idtable};
out.dimz = {fx.idtable.z1};
out.ntrials = [13; 10; 14; 10];   %per-cell medians of the fixture counts
out.elabels = {'E1';'E2';'E3';'E4'};
out.glabels = {'none'};
out.dod_labels = struct(...
    'inference_framework',...
    'four_gaussian_margins_all_residual_covariances_fixed_zero_v1',...
    'cut_posterior_status',...
    'no_dependence_module_exists_because_all_constituents_are_assumed_nonconcurrent',...
    'residual_dependence_structure',...
    'all six pairwise conditional residual covariances fixed exactly to zero',...
    'scale_parameterization','person_specific_log_residual_sd',...
    'estimand_version','signed_four_component_nonconcurrent_dod_v1');
REL.out = out;
end

function pd = localReportData(REL)
%Wrap the synthetic REL as the psyrat_data struct the report and viewer
%surfaces consume (the TestNonconcurrentDodReliability idiom).
pd = struct();
pd.rel = REL;
pd.ver = 'test';
pd.rel.filename = 'nonconcurrent_dod_gauss_fixture.psyrat';
pd.rel.nchains = 2;
pd.rel.niter = 600;
pd.rel.seed = 12345;
pd.rel.engine = 'cmdstan';
pd.rel.priors = psyrat_default_priors;
pd.rel.out.conv.converged = 1;
pd.rel.out.conv.ndivergent = 0;
end

function fig = localDodViewerFig(testCase, pd)
%Launch the dynamic-reliability viewer headlessly and return its figure
%(located by window name rather than by Tag).
name = 'Dynamic Reliability: View Results';
close(findall(0, 'Type', 'figure', 'Name', name));
psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);
fig = findall(0, 'Type', 'figure', 'Name', name);
testCase.assertNotEmpty(fig, 'Dynrel viewer figure not created');
fig = fig(1);
end

function cleaner = localInstallSaveShadow()
%Shadow uiputfile with a stub that auto-answers the save dialog from the
%PSYRAT_FAKE_SAVE global (the TestNonconcurrentDodReliability idiom).
shadowdir = tempname;
mkdir(shadowdir);
fid = fopen(fullfile(shadowdir, 'uiputfile.m'), 'w');
fprintf(fid, '%s\n', ...
    'function [filename,pathname,filterindex] = uiputfile(varargin)', ...
    '% TEST SHADOW: auto-answers the save dialog from a global.', ...
    'global PSYRAT_FAKE_SAVE', ...
    '[p,n,e] = fileparts(PSYRAT_FAKE_SAVE);', ...
    'filename = [n e]; pathname = [p filesep]; filterindex = 1;', ...
    'end');
fclose(fid);
addpath(shadowdir, '-begin');
rehash;
cleaner = onCleanup(@() localRemoveSaveShadow(shadowdir));
end

function localRemoveSaveShadow(shadowdir)
rmpath(shadowdir);
rehash;
if exist(shadowdir, 'dir')
    rmdir(shadowdir, 's');
end
end

function localFireDodSave(testCase, fig, outfile)
%Fire the 'Save variance table (CSV)' button through the PSYRAT_FAKE_SAVE
%test seam.
global PSYRAT_FAKE_SAVE %#ok<GVMIS> -- the uiputfile shadow's test seam
PSYRAT_FAKE_SAVE = outfile;
btns = findall(fig, 'Style', 'pushbutton');
fired = false;
for b = btns(:)'
    s = b.String; if iscell(s), s = strjoin(s, ' '); end
    if contains(string(s), 'Save variance')
        cb = b.Callback;
        cb(b, []);
        fired = true;
        break;
    end
end
PSYRAT_FAKE_SAVE = [];
testCase.assertTrue(fired, 'Save variance table button not found');
testCase.assertTrue(isfile(outfile), 'The variance CSV was not written');
end
