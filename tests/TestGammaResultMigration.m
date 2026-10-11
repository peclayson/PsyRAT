classdef TestGammaResultMigration < PsyRATTestBase
    % Offline tests for psyrat_reprocess_gamma_result: explicit, versioned
    % reprocessing of saved gamma results to the gamma_target_version 2 contract.
    %
    % EXPECTATIONS ARE INDEPENDENT of the production helper. The version-1
    % fixtures are built here from planted LOG-scale parameters through an
    % explicit forward map (support-union enumeration of the observed-scale
    % components, localSupportUnion), and the expected remainder R is taken
    % from that same enumeration, never from psyrat_gamma_mean_remainder_trt or
    % the extractor. Frozen anchors (SCIENTIFIC_CONTRACT.md "Frozen independent
    % anchors"): main-effects fixture alpha = log 5, s_p = .6, s_o = s_i = .5,
    % two-way SDs 0 gives R = 2.0652177285562118; LS old E = .04 ->
    % 2.1052177285562118; log-nu old C = 7.5507828568306812 ->
    % 9.6160005853868930; a fully removed facet gives R = 0.
    %
    % The inverse-compatibility cases and the declared error model come from
    % the frozen 2C inverse-compatibility set (local development evidence)
    % (120-digit reference, frozen before production comparison); the root
    % adversarial corruption set mirrors the 2C migration adversarial check.
    %
    % The G/D coefficient anchors and real-summary parity are not pinned here:
    % they depend on the summary path another worker owns (root lane).

    methods (Test)

        % ------------------------------------------------------------ static 25
        function testStatic25LognuTwoStrataMigrates(testCase)
            spec = localTwoStrataSpec();
            pd = localStatic25V1('lognu', spec);
            pd.relsummary = struct('stale', true);       % cached derived summary
            pd0 = pd;

            [upd, rec] = psyrat_reprocess_gamma_result(pd);

            testCase.verifyEqual(pd, pd0, 'Source struct was modified.');
            testCase.verifyEqual(rec.status, 'migrated');
            testCase.verifyEqual(rec.version_from, 1);
            testCase.verifyEqual(rec.version_to, 2);
            testCase.verifyEqual(rec.parameterization, 'lognu');
            testCase.verifyFalse(isfield(upd, 'relsummary'), 'Stale relsummary not removed.');
            testCase.verifyTrue(any(contains(rec.invalidation, 'relsummary')));

            o = upd.rel.out;
            gt = o.gamma_target;
            testCase.verifyEqual(gt.gamma_target_version, 2);
            testCase.verifyEqual(gt.surface_target, 'not_applicable');
            testCase.verifyEqual(gt.surface_count_source, 'not_applicable');
            testCase.verifyEqual(gt.reference_target, 'typical_person_zero_person_effects');
            testCase.verifyEqual(gt.subject_target, 'own_observation_variance_with_population_mean_components_at_own_z');
            testCase.verifyEqual(gt.remainder_scheme, 'population_highest_order_plus_person_observation_v2');
            testCase.verifyEqual(gt.parameterization, 'lognu');
            testCase.verifyEqual(gt.count_source, 'mean_trials_per_occupied_occasion');
            testCase.verifyEqual(gt.analysis, 'trt_sserrvar');
            testCase.verifyEqual(numel(gt.migration_history), 1);
            testCase.verifyEqual(numel(o.gamma_resid), 2);

            for s = 1:2
                st = spec(s);
                st.Ctyp = st.Ctyp_nu; st.Cp = st.Cp_nu;         % this fixture is log-nu
                Rexp = st.K.R;
                g = o.gamma_resid{s};
                testCase.verifyEqual(g.mean_remainder, Rexp, 'RelTol', 1e-10, ...
                    sprintf('Stratum %d: reconstructed R differs from the independent enumeration.', s));
                testCase.verifyEqual(g.obs_var_typ, st.Ctyp, 'RelTol', 1e-12);
                testCase.verifyEqual(g.obs_var_ss, st.Cp, 'RelTol', 1e-12);
                testCase.verifyEqual(g.parameterization, 'lognu');
                testCase.verifyEmpty(g.scale0, 'log-nu scale0 must be absent, not manufactured.');
                testCase.verifyEmpty(g.ind_scale, 'log-nu ind_scale must be absent, not manufactured.');
                popExp = 0.5 * log(Rexp + st.Ctyp);
                indExp = 0.5 * log(Rexp + st.Cp) - popExp;
                testCase.verifyEqual(o.pop_sdlog(:, s), popExp, 'RelTol', 1e-10);
                testCase.verifyEqual(cell2mat(o.ind_sdlog{s}), indExp, 'AbsTol', 1e-10);
                testCase.verifyEqual(exp(2 * (o.pop_sdlog(:, s) + cell2mat(o.ind_sdlog{s}))), Rexp + st.Cp, 'RelTol', 1e-10);
                testCase.verifyEqual(o.id_matches{s}, pd.rel.out.id_matches{s});
                testCase.verifyEqual(rec.id_mapping{s}, pd.rel.out.id_matches{s});
                % all three pair deficits are far above the error-model bound
                testCase.verifyEqual(rec.transforms(s).compatibility.n_unresolved, 0);
                testCase.verifyTrue(rec.transforms(s).reconstruction_exact_within_model);
            end
            testCase.verifyEqual(o.id_matches{2}.id2(:)', [3 1 2]);
            for f = {'mu', 'gro_sds', 'sig_occ', 'sig_trl', 'sig_trlxid', 'sig_occxid', 'sig_trlxocc', 'ind_bs', 'labels'}
                testCase.verifyEqual(o.(f{1}), pd.rel.out.(f{1}), sprintf('%s changed.', f{1}));
            end
            testCase.verifyEqual(o.gamma_disp.sigma_pot_e2_pop, spec(2).K.R + spec(2).Ctyp_nu, 'RelTol', 1e-10);
            testCase.verifyEqual(o.gamma_disp.sigma_pot_e2_pop_v1_observation_only, spec(2).Ctyp_nu);
            testCase.verifyTrue(rec.transforms(2).crosschecked_gamma_disp);
            testCase.verifyFalse(rec.transforms(1).crosschecked_gamma_disp);
        end

        function testStatic25LsMigratesWithExactScaleDraws(testCase)
            spec = localTwoStrataSpec();
            pd = localStatic25V1('ls', spec);
            [upd, rec] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(rec.status, 'migrated');
            testCase.verifyEqual(rec.parameterization, 'ls');
            o = upd.rel.out;
            testCase.verifyEqual(o.gamma_target.parameterization, 'ls');
            for s = 1:2
                st = spec(s);
                g = o.gamma_resid{s};
                testCase.verifyEqual(g.mean_remainder, st.K.R, 'RelTol', 1e-10);
                testCase.verifyEqual(g.scale0, st.ell0, 'AbsTol', 1e-14);
                testCase.verifyEqual(g.ind_scale, st.w_p, 'AbsTol', 1e-14);
                testCase.verifyEqual(g.parameterization, 'ls');
                testCase.verifyEqual(exp(2 * o.pop_sdlog(:, s)), st.K.R + exp(2 * st.ell0), 'RelTol', 1e-10);
                testCase.verifyEqual(exp(2 * (o.pop_sdlog(:, s) + cell2mat(o.ind_sdlog{s}))), ...
                    st.K.R + exp(2 * (st.ell0 + st.w_p)), 'RelTol', 1e-10);
            end
        end

        function testFrozenAnchorsMainEffectsFixture(testCase)
            st = localStratum(log(5), .6, .5, .5, 0, 0, 0, 3, 1:3, [301 302 303]);
            st.alpha(:) = log(5);
            st.K = localSupportUnion(st.alpha, st.v);
            testCase.verifyEqual(st.K.R, 2.0652177285562118 * ones(st.nd, 1), 'AbsTol', 1e-12);

            stLs = st; stLs.ell0 = 0.5 * log(0.04) * ones(st.nd, 1); stLs.w_p = zeros(st.nd, st.nsub);
            [uLs, rLs] = psyrat_reprocess_gamma_result(localStatic25V1('ls', localFinish('ls', stLs)));
            testCase.verifyEqual(exp(2 * uLs.rel.out.pop_sdlog(:, 1)), 2.1052177285562118 * ones(st.nd, 1), 'AbsTol', 1e-12);
            % zero two-way effects: deficits sit inside the bound and are recorded, not claimed exact
            testCase.verifyGreaterThan(rLs.transforms(1).compatibility.n_unresolved, 0);
            testCase.verifyFalse(rLs.transforms(1).reconstruction_exact_within_model);

            Ctarget = 7.5507828568306812;
            Vw = sum(st.v(2:6));
            stNu = st; stNu.b_nu = 2 * st.alpha + 2 * Vw - log(Ctarget / 2);
            stNu.u_p = zeros(st.nd, st.nsub); stNu.v_p = zeros(st.nd, st.nsub);
            [uNu, ~] = psyrat_reprocess_gamma_result(localStatic25V1('lognu', localFinish('lognu', stNu)));
            testCase.verifyEqual(uNu.rel.out.gamma_resid{1}.obs_var_typ, Ctarget * ones(st.nd, 1), 'RelTol', 1e-12);
            testCase.verifyEqual(exp(2 * uNu.rel.out.pop_sdlog(:, 1)), 9.6160005853868930 * ones(st.nd, 1), 'AbsTol', 1e-12);
        end

        function testZeroFacetGivesZeroRemainder(testCase)
            st = localStratum(log(4), .5, 0, .4, 0, .2, 0, 3, 1:3, [1 2 3]);
            st.ell0 = log(1.3) * ones(st.nd, 1); st.w_p = [.1 -.2 .05] .* ones(st.nd, 1);
            pd = localStatic25V1('ls', localFinish('ls', st));
            [upd, rec] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(upd.rel.out.gamma_resid{1}.mean_remainder, zeros(st.nd, 1));
            testCase.verifyEqual(upd.rel.out.pop_sdlog, pd.rel.out.pop_sdlog, 'AbsTol', 1e-14);
            testCase.verifyEqual(cell2mat(upd.rel.out.ind_sdlog{1}), cell2mat(pd.rel.out.ind_sdlog{1}), 'AbsTol', 1e-14);
            testCase.verifyEqual(rec.status, 'migrated');
        end

        function testRepeatedMigrationDoesNotPoolTwice(testCase)
            for param = {'lognu', 'ls'}
                pd = localStatic25V1(param{1}, localTwoStrataSpec());
                [u1, r1] = psyrat_reprocess_gamma_result(pd);
                [u2, r2] = psyrat_reprocess_gamma_result(u1);
                testCase.verifyEqual(r1.status, 'migrated');
                testCase.verifyEqual(r2.status, 'already_current');
                testCase.verifyEqual(r2.version_from, 2);
                testCase.verifyEqual(u2, u1, 'Second call changed the result (double pooling or restamp).');
                testCase.verifyEqual(numel(u2.rel.out.gamma_target.migration_history), 1);
            end
        end

        % ------------------------------------------------------------ frozen inverse cases
        function testFrozenInverseCompatibilityCases(testCase)
            % retained normalized components (mu = 1, B = 1) from
            % frozen_inverse_compatibility.json, order P O I PO PI OI
            valid_positive = [.2 .3 .1 .1224 .0464 .0729];
            valid_zero_pairs = [.00002 .00002 .1 4e-10 .000002 .000002];
            invalid_missing = [.00002 .00002 .1 0 .000002 .000002];
            invalid_deficit = [.001 .001 .1 9.995e-7 .0001 .0001];
            legacy_unresolved = [1e-9 1e-9 .1 0 1e-10 1e-10];
            valid_tiny = [.00002 .00002 .1 5.0000400004e-10 .000002 .000002];

            % exact compatible positive pairs: all deficits above the bound
            [u, r] = psyrat_reprocess_gamma_result(localStaticFromNormalized(valid_positive));
            testCase.verifyEqual(r.status, 'migrated');
            testCase.verifyEqual(r.transforms(1).compatibility.n_unresolved, 0);
            q = [.2 .3 .1 .04 .02 .03];
            testCase.verifyEqual(u.rel.out.gamma_resid{1}.mean_remainder, localRfromQ(q), 'RelTol', 1e-10);

            % exact zero pairs: deficit 0 inside the bound; PO recorded unresolved; bound matches the frozen model
            [u, r] = psyrat_reprocess_gamma_result(localStaticFromNormalized(valid_zero_pairs));
            testCase.verifyEqual(r.status, 'migrated');
            c = r.transforms(1).compatibility;
            testCase.verifyTrue(any(strcmp({c.unresolved.pair}, 'PO')));
            poRow = c.unresolved(strcmp({c.unresolved.pair}, 'PO'));
            testCase.verifyEqual(poRow(1).bound, 3.11326584849126220419e-14, 'RelTol', 1e-6, ...
                'PO bound does not reproduce the frozen declared error model.');
            testCase.verifyEqual(u.rel.out.gamma_resid{1}.mean_remainder, localRfromQ([.00002 .00002 .1 0 0 0]), 'RelTol', 1e-9);

            % PO stored as 0 instead of its required 4e-10: materially incompatible, REFUSED
            testCase.verifyError(@() psyrat_reprocess_gamma_result(localStaticFromNormalized(invalid_missing)), ...
                'psyrat_reprocess_gamma_result:invalid_inverse');
            % PO 5e-10 below its required 1e-6: REFUSED
            testCase.verifyError(@() psyrat_reprocess_gamma_result(localStaticFromNormalized(invalid_deficit)), ...
                'psyrat_reprocess_gamma_result:invalid_inverse');

            % legacy roundoff: PO loses 1e-18 (a = b = 1e-9); inside the bound, canonicalized and recorded, not exact
            [u, r] = psyrat_reprocess_gamma_result(localStaticFromNormalized(legacy_unresolved));
            testCase.verifyEqual(r.status, 'migrated');
            c = r.transforms(1).compatibility;
            testCase.verifyGreaterThan(c.n_unresolved_negative_canonicalized, 0);
            testCase.verifyFalse(r.transforms(1).reconstruction_exact_within_model);
            testCase.verifyGreaterThan(u.rel.out.gamma_resid{1}.mean_remainder, 0);

            % tiny positive pair 1e-10 above the bound: retained, not clamped
            [u, r] = psyrat_reprocess_gamma_result(localStaticFromNormalized(valid_tiny));
            testCase.verifyEqual(r.status, 'migrated');
            c = r.transforms(1).compatibility;
            testCase.verifyFalse(any(strcmp({c.unresolved.pair}, 'PO')), 'A positive q_PO above the bound must not be treated as unresolved.');
            testCase.verifyEqual(u.rel.out.gamma_resid{1}.mean_remainder, localRfromQ([.00002 .00002 .1 1e-10 0 0]), 'RelTol', 1e-6);

            % raw last-stratum cross-check: gamma_disp claims vPO twice the true value;
            % the implied normalized PO differs by 1e-10, beyond the two-sided bound
            pd = localStaticFromNormalized(valid_tiny);
            trueV = log1p([.00002 .00002 .1 1e-10 0 0]);
            pd.rel.out.gamma_disp = struct('alpha', -sum(trueV) / 2, ...
                's_o', sqrt(trueV(2)), 's_t', sqrt(trueV(3)), 's_po', sqrt(2 * trueV(4)), ...
                's_pt', 0, 's_ot', 0, 'log_sigma', 0, 's_sd', .3, 'rho_ps', 0, ...
                'disp_factor', 1, 'sigma_pot_e2_pop', 1, 'is_ls', true);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), ...
                'psyrat_reprocess_gamma_result:inverse_crosscheck_failed');
            % and the consistent raw claim passes the same check
            pd.rel.out.gamma_disp.s_po = sqrt(trueV(4));
            [~, r] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyTrue(r.transforms(1).crosschecked_gamma_disp);
        end

        function testInvalidInverseRefused(testCase)
            st = localStratum(log(5), .6, .5, .5, 0, 0, 0, 3, 1:3, [1 2 3]);
            st.ell0 = zeros(st.nd, 1); st.w_p = zeros(st.nd, st.nsub);
            pd = localStatic25V1('ls', localFinish('ls', st));
            pd.rel.out.sig_occxid(:) = 0;
            pd.rel.out = rmfield(pd.rel.out, 'gamma_disp');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_inverse');
        end

        function testCrossCheckRefusesInconsistentGammaDisp(testCase)
            pd = localStatic25V1('lognu', localTwoStrataSpec());
            pd.rel.out.gamma_disp.s_o = pd.rel.out.gamma_disp.s_o + 0.05;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:inverse_crosscheck_failed');
        end

        % ------------------------------------------------------------ version-2 validation
        function testVersion2EncodingMismatchRefused(testCase)
            % frozen small_residual_mismatch: retained E = 1e-12, claimed R + C = 2e-12
            [u1, ~] = psyrat_reprocess_gamma_result(localStatic25V1('ls', localTwoStrataSpec()));
            bad = u1;
            nd = size(bad.rel.out.mu, 1);
            bad.rel.out.gamma_resid{1}.mean_remainder = 1e-12 * ones(nd, 1);
            bad.rel.out.gamma_resid{1}.obs_var_typ = 1e-12 * ones(nd, 1);
            bad.rel.out.pop_sdlog(:, 1) = 0.5 * log(1e-12);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_v2_claim');
            % NaN observation component on a current result: refused before any comparison
            bad = u1; bad.rel.out.gamma_resid{1}.obs_var_ss(:) = NaN;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_source_values');
            % non-finite encoding on a current result
            bad = u1; bad.rel.out.ind_sdlog{2}{1, 1} = Inf;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_source_values');
            % duplicate id2 on a current result: refused before the current return
            bad = u1; bad.rel.out.id_matches{1}.id2(2) = bad.rel.out.id_matches{1}.id2(1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_id_matches');
            % duplicate subject label on a current result
            bad = u1; bad.rel.out.id_matches{1}.id(2) = bad.rel.out.id_matches{1}.id(1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_ids');
            % R doubled: encoding identity broken
            bad = u1; bad.rel.out.gamma_resid{1}.mean_remainder = 2 * bad.rel.out.gamma_resid{1}.mean_remainder;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_v2_claim');
            % version-2 stamp without raw pieces
            bad = localStatic25V1('lognu', localTwoStrataSpec());
            bad.rel.out.gamma_target = u1.rel.out.gamma_target;
            bad.rel.out.gamma_target.parameterization = 'lognu';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:invalid_v2_claim');
        end

        function testStampIdentifierContradictionsRefused(testCase)
            % static: every fixed identifier, an unsupported version, and parameterization contradictions
            [u1, ~] = psyrat_reprocess_gamma_result(localStatic25V1('lognu', localTwoStrataSpec()));
            for f = {'surface_target', 'reference_target', 'subject_target', 'remainder_scheme', ...
                    'count_source', 'surface_count_source', 'analysis'}
                bad = u1; bad.rel.out.gamma_target.(f{1}) = 'wrong';
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), ...
                    'psyrat_reprocess_gamma_result:contradictory_stamp', sprintf('field %s', f{1}));
                bad = u1; bad.rel.out.gamma_target = rmfield(bad.rel.out.gamma_target, f{1});
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), ...
                    'psyrat_reprocess_gamma_result:contradictory_stamp', sprintf('missing %s', f{1}));
            end
            for v = {99, 1, 2.5, NaN, '2'}
                bad = u1; bad.rel.out.gamma_target.gamma_target_version = v{1};
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), ...
                    'psyrat_reprocess_gamma_result:unsupported_target_version');
            end
            bad = u1; bad.rel.out.gamma_target.parameterization = 'ls';     % vs gammascale 1
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:contradictory_parameterization');
            bad = u1; bad.rel.out.gamma_resid{2}.parameterization = 'ls';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:contradictory_parameterization');
            for g = {3, NaN, [1 2], 'lognu'}
                bad = u1; bad.rel.gammascale = g{1};
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:unknown_parameterization');
            end
            % gammascale absent: gamma_disp (log-nu fields) vs an LS stamp contradict
            bad = u1; bad.rel = rmfield(bad.rel, 'gammascale'); bad.rel.out.gamma_target.parameterization = 'ls';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:contradictory_parameterization');
            % gammascale absent, LS gamma_disp (log_sigma/s_sd) identifies LS unambiguously
            pdLs = localStatic25V1('ls', localTwoStrataSpec()); pdLs.rel = rmfield(pdLs.rel, 'gammascale');
            [~, r] = psyrat_reprocess_gamma_result(pdLs);
            testCase.verifyEqual(r.parameterization, 'ls');

            % dynamic: fixed identifiers, then a legitimate explicit target change
            [d1, ~] = psyrat_reprocess_gamma_result(localDynFixture(1, 2, 1, 2, 2));
            for f = {'reference_target', 'subject_target', 'remainder_scheme', 'parameterization', ...
                    'count_source', 'surface_count_source', 'analysis'}
                bad = d1; bad.rel.out.gamma_target.(f{1}) = 'wrong';
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), ...
                    'psyrat_reprocess_gamma_result:contradictory_stamp', sprintf('dynamic field %s', f{1}));
            end
            bad = d1; bad.rel.out.gamma_target.surface_target = 'wrong';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:contradictory_stamp');
            bad = d1; bad.rel.out.gamma_target.gamma_target_version = 99;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:unsupported_target_version');
            [d2, r2] = psyrat_reprocess_gamma_result(d1, 'surface_target', 'typical_person');
            testCase.verifyEqual(r2.status, 'target_restamped');
            testCase.verifyEqual(d2.rel.out.gamma_target.surface_target, 'typical_person');
        end

        % ------------------------------------------------------------ root adversarial corruption set
        function testRootAdversarialCorruptionsRefused(testCase)
            % Mirrors the 2C migration adversarial check (local development evidence)
            % (17 cases; the three "current static" cases use this class's migrated
            % LS/log-nu fixtures in place of the retained capture).
            good = localDynFixture(1, 2, 1, 2, 2);
            cases = {}; labels = {}; ids = {};
            p = good; p.rel.out.mu(1, 1) = NaN;                 add('nonfinite_mu', p, 'invalid_source_values');
            p = good; p.rel.out.b{1}{1, 1} = Inf;               add('nonfinite_mean_slope', p, 'invalid_source_values');
            p = good; p.rel.out.b_nu{2}{1, 1} = NaN;            add('nonfinite_scale_slope', p, 'invalid_source_values');
            p = good; p.rel.out.sig_occ(1, 2) = NaN;            add('nonfinite_facet_sd', p, 'invalid_source_values');
            p = good; p.rel.out.sig_trl(1, 1) = -.2;            add('negative_facet_sd', p, 'invalid_source_values');
            p = good; p.rel.out.pop_lognu(1, 1) = Inf;          add('nonfinite_scale_intercept', p, 'invalid_source_values');
            p = good; p.rel.out.labels = p.rel.out.labels(1);   add('missing_stratum_label', p, 'invalid_strata');
            p = good; p.rel.out.ntrials = p.rel.out.ntrials(1); add('missing_stratum_trial_count', p, 'invalid_strata');
            p = good; p.rel.out.ssinfo{2}.z2(1) = NaN;          add('nonfinite_own_z2', p, 'invalid_ssinfo');
            p = good; p.rel.out.chol_corrmat{1} = {1};          add('malformed_correlation_draws', p, 'invalid_correlation_draws');
            p = good; p.rel.gammascale = 3;                     add('unknown_parameterization', p, 'unknown_parameterization');
            p = psyrat_reprocess_gamma_result(good); p.rel.out.gamma_target.gamma_target_version = 99;
            add('unsupported_future_target_version', p, 'unsupported_target_version');
            p = psyrat_reprocess_gamma_result(good); p.rel.out.gamma_target.subject_target = 'wrong';
            add('contradictory_subject_target', p, 'contradictory_stamp');
            p = good; p.rel.out.ssinfo{1}.id(2) = p.rel.out.ssinfo{1}.id(1);
            add('duplicate_subject_id', p, 'invalid_ids');
            cur = psyrat_reprocess_gamma_result(localStatic25V1('lognu', localTwoStrataSpec()));
            p = cur; p.rel.out.gamma_resid{1}.obs_var_ss(:) = NaN; add('current_static_nonfinite_C', p, 'invalid_source_values');
            p = cur; p.rel.out.id_matches{1}.id2(2) = p.rel.out.id_matches{1}.id2(1);
            add('current_static_duplicate_id2', p, 'invalid_id_matches');
            p = cur; p.rel.out.gamma_target.parameterization = 'ls';
            add('current_static_contradictory_parameterization', p, 'contradictory_parameterization');
            testCase.assertEqual(numel(cases), 17);
            for k = 1:numel(cases)
                testCase.verifyError(@() psyrat_reprocess_gamma_result(cases{k}), ...
                    ['psyrat_reprocess_gamma_result:' ids{k}], labels{k});
            end
            function add(label, value, id)
                labels{end+1} = label; cases{end+1} = value; ids{end+1} = id;
            end
        end

        % ------------------------------------------------------------ refusals
        function testMissingFieldRefusedWithExactName(testCase)
            pd = localStatic25V1('lognu', localTwoStrataSpec());
            pd.rel.out = rmfield(pd.rel.out, 'sig_trlxocc');
            try
                psyrat_reprocess_gamma_result(pd);
                testCase.verifyFail('Expected a missing_field error.');
            catch err
                testCase.verifyEqual(err.identifier, 'psyrat_reprocess_gamma_result:missing_field');
                testCase.verifyTrue(contains(err.message, 'sig_trlxocc'), 'Message must name the missing field.');
            end
        end

        function testMisalignedDrawsSubjectsAndIdsRefused(testCase)
            base = localStatic25V1('lognu', localTwoStrataSpec());
            pd = base; pd.rel.out.ind_sdlog{2} = pd.rel.out.ind_sdlog{2}(1:end-1, :);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:misaligned_draws');
            pd = base; pd.rel.out.ind_sdlog{2} = pd.rel.out.ind_sdlog{2}(:, 1:2);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:misaligned_subjects');
            pd = base; pd.rel.out.id_matches{2}.id2 = [1; 1; 2];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_id_matches');
            pd = base; pd.rel.out.sig_occ = pd.rel.out.sig_occ(:, 1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:misaligned_draws');
            pd = base; pd.rel.out.labels = {'only one'};
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_strata');
            pd = base; pd.rel.out.mu(2, 1) = -1;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_source_values');
            pd = base; pd.rel.out.pop_sdlog(1, 2) = NaN;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_source_values');
        end

        function testAmbiguousParameterizationRefused(testCase)
            pd = localStatic25V1('lognu', localTwoStrataSpec());
            pd.rel = rmfield(pd.rel, 'gammascale');
            pd.rel.out = rmfield(pd.rel.out, 'gamma_disp');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:ambiguous_parameterization');
        end

        function testUnsupportedFamilyAndAnalysisRefused(testCase)
            pd = localStatic25V1('lognu', localTwoStrataSpec());
            g = pd; g.rel.family = 'gaussian';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(g), 'psyrat_reprocess_gamma_result:unsupported_family');
            d = pd; d.rel.analysis = 'trt_diff';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(d), 'psyrat_reprocess_gamma_result:unsupported_analysis');
            s6 = pd; s6.rel.analysis = 'ic_sserrvar';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(s6), 'psyrat_reprocess_gamma_result:unsupported_analysis');
        end

        function testSummaryOnlyInputRefused(testCase)
            d = testCase.scratchDir();
            for ext = {'.xlsx', '.csv', '.xls', '.txt'}
                f = fullfile(d, ['summary' ext{1}]);
                testCase.verifyError(@() psyrat_reprocess_gamma_result(f, 'output', fullfile(d, 'o.psyrat')), ...
                    'psyrat_reprocess_gamma_result:summary_only_input');
            end
        end

        % ------------------------------------------------------------ files
        function testFileRoundTripPreservesSourceAndWrapper(testCase)
            d = testCase.scratchDir();
            src = fullfile(d, 'legacy.psyrat');
            out = fullfile(d, 'migrated.psyrat');
            psyrat_data = localStatic25V1('lognu', localTwoStrataSpec()); %#ok<NASGU>
            sidecar_note = 'keep me'; %#ok<NASGU>
            save(src, 'psyrat_data', 'sidecar_note');
            bytesBefore = localReadBytes(src);

            [upd, rec] = psyrat_reprocess_gamma_result(src, 'output', out);

            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'Source file bytes changed.');
            testCase.verifyTrue(isfile(out));
            S = load(out, '-mat');
            testCase.verifyTrue(all(ismember({'psyrat_data', 'sidecar_note'}, fieldnames(S))));
            testCase.verifyEqual(S.sidecar_note, 'keep me');
            testCase.verifyEqual(S.psyrat_data.rel.out.gamma_target.gamma_target_version, 2);
            testCase.verifyEqual(S.psyrat_data, upd);
            testCase.verifyEqual(rec.output, out);
            testCase.verifyEqual(S.psyrat_data.rel.out.gamma_target.migration_history(end).output, out);

            outputBefore = localReadBytes(out);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(src, 'output', out), 'psyrat_reprocess_gamma_result:output_exists');
            testCase.verifyEqual(localReadBytes(out), outputBefore, 'Refusal must preserve the existing destination.');
            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'Refusal must preserve the source.');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(src, 'output', src), 'psyrat_reprocess_gamma_result:output_is_source');
            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'In-place refusal must preserve the source.');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(src), 'psyrat_reprocess_gamma_result:output_required');
            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'Missing-output refusal must preserve the source.');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(src, 'output', fullfile(d, 'bad.xlsx')), 'psyrat_reprocess_gamma_result:bad_output');
            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'Refusals must not touch the source.');
        end

        function testRelativeSourceAndAbsoluteOutputRefusedAsSource(testCase)
            d = testCase.scratchDir();
            src = fullfile(d, 'legacy.mat');
            psyrat_data = localStatic25V1('lognu', localTwoStrataSpec()); %#ok<NASGU>
            save(src, 'psyrat_data');
            bytesBefore = localReadBytes(src);
            originalFolder = pwd;
            restoreFolder = onCleanup(@() cd(originalFolder)); %#ok<NASGU>
            cd(d);
            testCase.verifyError(@() psyrat_reprocess_gamma_result('legacy.mat', 'output', src), ...
                'psyrat_reprocess_gamma_result:output_is_source');
            testCase.verifyEqual(localReadBytes(src), bytesBefore, 'Refusal must preserve the source.');
        end

        function testFileWithoutPsyratDataRefused(testCase)
            d = testCase.scratchDir();
            f = fullfile(d, 'other.mat');
            something = 1; %#ok<NASGU>
            save(f, 'something');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(f, 'output', fullfile(d, 'o.psyrat')), ...
                'psyrat_reprocess_gamma_result:no_psyrat_data');
        end

        function testBareRelInputCannotBeWritten(testCase)
            pd = localStatic25V1('ls', localTwoStrataSpec());
            [updRel, rec] = psyrat_reprocess_gamma_result(pd.rel);
            testCase.verifyTrue(isfield(updRel, 'out') && ~isfield(updRel, 'rel'));
            testCase.verifyEqual(rec.source_kind, 'rel');
            testCase.verifyEmpty(rec.output);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd.rel, 'output', fullfile(testCase.scratchDir(), 'x.psyrat')), ...
                'psyrat_reprocess_gamma_result:output_not_allowed');
        end

        % ------------------------------------------------------------ dynamic
        function testDynamicRoutesStampWithoutNumericChange(testCase)
            expectScheme = struct('dyn11', 'population_mean_remainder_plus_selected_surface_observation_v1', ...
                'dyn14', 'population_mean_remainder_plus_selected_surface_observation_v1', ...
                'dyn26', 'population_mean_remainder_plus_person_observation_v1', ...
                'dyn27', 'population_highest_order_plus_person_observation_v2');
            for scale = 1:2
                for facets = 1:2
                    for subjects = 0:1
                        for ndim = 1:2
                            for nstrata = 1:2
                                pd = localDynFixture(scale, facets, subjects, ndim, nstrata);
                                pd0 = pd;
                                label = sprintf('scale%d facet%d subj%d dim%d strata%d', scale, facets, subjects, ndim, nstrata);
                                [u1, r1] = psyrat_reprocess_gamma_result(pd);
                                testCase.verifyEqual(pd, pd0, [label ': source modified']);
                                testCase.verifyEqual(r1.status, 'migrated', label);
                                gt = u1.rel.out.gamma_target;
                                testCase.verifyEqual(gt.gamma_target_version, 2, label);
                                testCase.verifyEqual(gt.surface_target, 'population_marginal', label);
                                testCase.verifyEqual(gt.surface_count_source, 'default_stratum_median', label);
                                testCase.verifyEqual(gt.reference_target, 'not_applicable', label);
                                testCase.verifyEqual(gt.remainder_scheme, expectScheme.(r1.route), label);
                                testCase.verifyEqual(gt.analysis, pd.rel.analysis, label);
                                if scale == 2, testCase.verifyEqual(gt.parameterization, 'ls', label);
                                else, testCase.verifyEqual(gt.parameterization, 'lognu', label); end
                                if subjects
                                    testCase.verifyEqual(gt.subject_target, 'own_observation_variance_with_population_mean_components_at_own_z', label);
                                    testCase.verifyEqual(gt.count_source, 'median_trials_per_occupied_occasion', label);
                                    for s = 1:nstrata
                                        testCase.verifyEqual(r1.id_mapping{s}, pd.rel.out.ssinfo{s}(:, {'id', 'id2'}), label);
                                    end
                                else
                                    testCase.verifyEqual(gt.subject_target, 'not_applicable', label);
                                    testCase.verifyEqual(gt.count_source, 'not_applicable', label);
                                end
                                testCase.verifyEqual(rmfield(u1.rel.out, 'gamma_target'), pd.rel.out, [label ': numerics changed']);
                                [u2, r2] = psyrat_reprocess_gamma_result(u1);
                                testCase.verifyEqual(r2.status, 'already_current', label);
                                testCase.verifyEqual(u2, u1, [label ': second call changed the result']);
                                [u3, r3] = psyrat_reprocess_gamma_result(u1, 'surface_target', 'typical_person');
                                testCase.verifyEqual(r3.status, 'target_restamped', label);
                                testCase.verifyEqual(u3.rel.out.gamma_target.surface_target, 'typical_person', label);
                                testCase.verifyEqual(rmfield(u3.rel.out, 'gamma_target'), pd.rel.out, [label ': restamp changed numerics']);
                                testCase.verifyEqual(numel(u3.rel.out.gamma_target.migration_history), 2, label);
                            end
                        end
                    end
                end
            end
        end

        function testDynamicSubjectRouteRequiresOwnEffects(testCase)
            pd = localDynFixture(1, 2, 1, 1, 2);
            pd.rel.out = rmfield(pd.rel.out, 'ind_nu');
            try
                psyrat_reprocess_gamma_result(pd);
                testCase.verifyFail('Expected missing_field.');
            catch err
                testCase.verifyEqual(err.identifier, 'psyrat_reprocess_gamma_result:missing_field');
                testCase.verifyTrue(contains(err.message, 'ind_nu'));
            end
            pd = localDynFixture(2, 1, 1, 1, 1);
            pd.rel.out.ssinfo{1}.id2 = [1; 1; 3];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_id_matches');
            pd = localDynFixture(1, 1, 1, 1, 1);
            pd.rel.out.ind_bs{1} = pd.rel.out.ind_bs{1}(:, 1:2);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:misaligned_subjects');
            pd = localDynFixture(1, 2, 0, 2, 1);
            pd.rel.out.b{1} = pd.rel.out.b{1}(:, 1:2);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:misaligned_draws');
            pd = localDynFixture(1, 1, 1, 1, 1);
            pd.rel.out.ssinfo{1}.trls(2) = 0;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ssinfo');
            pd = localDynFixture(1, 1, 1, 1, 1);
            pd.rel.out.ssinfo{1}.z1(3) = Inf;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ssinfo');
        end

        function testDynamicDimensionAndCountConsistencyRefused(testCase)
            pd = localDynFixture(1, 2, 1, 2, 2); pd.rel.ndim = 3;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_dimensions');
            pd = localDynFixture(1, 2, 1, 2, 2); pd.rel.dim_names = pd.rel.dim_names(1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_dimensions');
            pd = localDynFixture(1, 2, 1, 2, 2); pd.rel.out.dimz{2}(1, 2) = NaN;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_dimensions');
            pd = localDynFixture(1, 2, 1, 2, 2); pd.rel.out.dimz{2} = pd.rel.out.dimz{2}(:, 1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_dimensions');
            pd = localDynFixture(1, 2, 1, 1, 2); pd.rel.out.nocc = pd.rel.out.nocc(1);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_strata');
            pd = localDynFixture(1, 2, 1, 1, 2); pd.rel.out.ntrials(2) = 0;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_strata');
            pd = localDynFixture(1, 2, 1, 1, 2); pd.rel.out.labels{2} = '';
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_strata');
            pd = localDynFixture(2, 1, 0, 1, 1); pd.rel.out.gro_sds{1}{1, 2} = -0.1;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_source_values');
            % a valid two-strata result with its own id2 permutation remains usable
            pd = localDynFixture(2, 2, 1, 2, 2);
            [~, r] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(r.status, 'migrated');
            testCase.verifyEqual(r.id_mapping{2}.id2(:)', [3 1 2]);
        end

        function testCorrelationDrawsRequiredForBothScalesAndTargets(testCase)
            % complete draw-aligned chol_corrmat is required under log-nu AND
            % location-scale, for population_marginal AND typical_person
            for scale = 1:2
                for target = {'population_marginal', 'typical_person'}
                    pd = localDynFixture(scale, 1, 0, 1, 1);
                    pd.rel.out = rmfield(pd.rel.out, 'chol_corrmat');
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd, 'surface_target', target{1}), ...
                        'psyrat_reprocess_gamma_result:missing_field', sprintf('scale %d %s', scale, target{1}));
                    pd = localDynFixture(scale, 2, 1, 1, 2);
                    pd.rel.out.chol_corrmat{2} = {1};
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd, 'surface_target', target{1}), ...
                        'psyrat_reprocess_gamma_result:invalid_correlation_draws');
                    pd = localDynFixture(scale, 1, 1, 1, 1);
                    pd.rel.out.chol_corrmat{1}{1, 2, 2} = 0.5;          % row norm broken
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd, 'surface_target', target{1}), ...
                        'psyrat_reprocess_gamma_result:invalid_correlation_draws');
                    pd = localDynFixture(scale, 1, 1, 1, 1);
                    pd.rel.out.chol_corrmat{1}{1, 2, 1} = NaN;
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd, 'surface_target', target{1}), ...
                        'psyrat_reprocess_gamma_result:invalid_correlation_draws');
                end
            end
            % rho = +-1 is a valid boundary (L22 = 0) in the estimator's num2cell form
            for rho = [-1 1]
                pd = localDynFixture(1, 1, 0, 1, 1);
                nd = size(pd.rel.out.mu, 1);
                L = zeros(nd, 2, 2); L(:, 1, 1) = 1; L(:, 2, 1) = rho; L(:, 2, 2) = 0;
                pd.rel.out.chol_corrmat{1} = num2cell(L);
                [~, r] = psyrat_reprocess_gamma_result(pd);
                testCase.verifyEqual(r.status, 'migrated');
                testCase.verifyLessThan(r.transforms(1).chol_max_deviation, 1e-12);
                % the same block stored numerically breaks the consumers' cell2mat contract
                pd.rel.out.chol_corrmat{1} = L;
                testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_storage');
            end
        end

        function testNumericInnerStorageRefused(testCase)
            % every per-stratum block the consumers read through cell2mat must
            % keep the estimator's num2cell form; a numeric inner block is
            % refused as unsupported storage (both scales, dynamic and static)
            for scale = 1:2
                if scale == 2, slopeScale = 'b_sigma'; own = {'ind_sdlog'}; else, slopeScale = 'b_nu'; own = {'ind_bs', 'ind_nu'}; end
                for f = [{'gro_sds', 'b', slopeScale, 'chol_corrmat'}, own]
                    pd = localDynFixture(scale, 2, 1, 1, 2);
                    pd.rel.out.(f{1}){2} = localCellToNumeric(pd.rel.out.(f{1}){2});
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), ...
                        'psyrat_reprocess_gamma_result:unsupported_storage', sprintf('scale %d %s', scale, f{1}));
                end
            end
            for param = {'lognu', 'ls'}
                for f = {'gro_sds', 'ind_sdlog', 'ind_bs'}
                    pd = localStatic25V1(param{1}, localTwoStrataSpec());
                    pd.rel.out.(f{1}){1} = cell2mat(pd.rel.out.(f{1}){1});
                    testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), ...
                        'psyrat_reprocess_gamma_result:unsupported_storage', sprintf('%s %s', param{1}, f{1}));
                end
                % and on an already-current result the same contract holds
                [u1, ~] = psyrat_reprocess_gamma_result(localStatic25V1(param{1}, localTwoStrataSpec()));
                u1.rel.out.ind_sdlog{2} = cell2mat(u1.rel.out.ind_sdlog{2});
                testCase.verifyError(@() psyrat_reprocess_gamma_result(u1), 'psyrat_reprocess_gamma_result:unsupported_storage');
            end
        end

        function testEmptyOrPartialProvenanceRefused(testCase)
            % a present-but-empty marker is not "unstamped": refusing it is what
            % prevents a second addition of R on an already-migrated result
            [u1, ~] = psyrat_reprocess_gamma_result(localStatic25V1('ls', localTwoStrataSpec()));
            for marker = {[], {}, '', struct(), struct('surface_target', 'not_applicable')}
                bad = u1; bad.rel.out.gamma_target = marker{1};
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), ...
                    'psyrat_reprocess_gamma_result:unsupported_target_version', class(marker{1}));
            end
            % stamp removed but the version-2 raw pieces retained: inconsistent provenance,
            % never reinterpreted as an observation-only encoding
            bad = u1; bad.rel.out = rmfield(bad.rel.out, 'gamma_target');
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:inconsistent_provenance');
            % dynamic: present-but-empty marker refused as well
            [d1, ~] = psyrat_reprocess_gamma_result(localDynFixture(1, 2, 1, 1, 2));
            for marker = {[], struct()}
                bad = d1; bad.rel.out.gamma_target = marker{1};
                testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:unsupported_target_version');
            end
            bad = localDynFixture(2, 1, 1, 1, 1); bad.rel.out.gamma_resid = {struct()};
            testCase.verifyError(@() psyrat_reprocess_gamma_result(bad), 'psyrat_reprocess_gamma_result:inconsistent_provenance');
        end

        function testSsinfoAndIdShapesRefused(testCase)
            % ssinfo columns must be exactly one value per row; a NaN/negative first
            % column beside a valid second column must not pass a flattened check
            pd = localDynFixture(1, 2, 1, 1, 1);
            z = pd.rel.out.ssinfo{1}.z1; pd.rel.out.ssinfo{1}.z1 = [[NaN; z(2:3)], z];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ssinfo');
            pd = localDynFixture(1, 2, 1, 1, 1);
            t = pd.rel.out.ssinfo{1}.trls; pd.rel.out.ssinfo{1}.trls = [[-1; t(2:3)], t];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ssinfo');
            pd = localDynFixture(1, 2, 1, 2, 1);
            z2 = pd.rel.out.ssinfo{1}.z2; pd.rel.out.ssinfo{1}.z2 = [z2, z2];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ssinfo');
            % IDs: one per row; multi-column numeric IDs and char matrices are refused
            pd = localDynFixture(1, 2, 1, 1, 1);
            id = pd.rel.out.ssinfo{1}.id; pd.rel.out.ssinfo{1}.id = [id, id + 100];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ids');
            pd = localDynFixture(1, 2, 1, 1, 1);
            pd.rel.out.ssinfo{1}.id = ['S1'; 'S2'; 'S3'];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ids');
            pd = localDynFixture(1, 2, 1, 1, 1);
            id2 = pd.rel.out.ssinfo{1}.id2; pd.rel.out.ssinfo{1}.id2 = [id2, id2];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_id_matches');
            % static id_matches: same shape contract
            pd = localStatic25V1('ls', localTwoStrataSpec());
            id = pd.rel.out.id_matches{2}.id; pd.rel.out.id_matches{2}.id = [id, id];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ids');
            pd = localStatic25V1('ls', localTwoStrataSpec());
            id2 = pd.rel.out.id_matches{2}.id2; pd.rel.out.id_matches{2}.id2 = [id2, id2];
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_id_matches');
            % supported ID types remain usable: cellstr, string and categorical columns
            pd = localDynFixture(1, 2, 1, 1, 1);
            pd.rel.out.ssinfo{1}.id = {'S401'; 'S402'; 'S403'};
            [~, r] = psyrat_reprocess_gamma_result(pd); testCase.verifyEqual(r.status, 'migrated');
            pd.rel.out.ssinfo{1}.id = string({'S401'; 'S402'; 'S403'});
            [~, r] = psyrat_reprocess_gamma_result(pd); testCase.verifyEqual(r.status, 'migrated');
            pd.rel.out.ssinfo{1}.id = categorical({'S401'; 'S402'; 'S403'});
            [~, r] = psyrat_reprocess_gamma_result(pd); testCase.verifyEqual(r.status, 'migrated');
            pd.rel.out.ssinfo{1}.id = {'S401'; 'S401'; 'S403'};
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:invalid_ids');
        end

        function testStrongCancellationRepeatIsNoOp(testCase)
            % Astra-confirmed valid round trip: all six observed SDs equal mu =
            % sqrt(expm1(.1)), so a = b = c = 1 and all pair q are 0; R = mu^2 =
            % .10517091807564755. Legacy pop = 100, ind = -100 (C_typ = e^200,
            % C_p = 1). After migration ind_v2 = .5 log(R + 1) - 100 cancels
            % against pop = 100 in the validator, which the encoding bound must
            % cover; the repeat must be a byte-equal no-op.
            m = sqrt(expm1(.1));
            o = struct();
            o.mu = m; o.pop_sdlog = 100;
            o.sig_occ = m; o.sig_trl = m; o.sig_occxid = m; o.sig_trlxid = m; o.sig_trlxocc = m;
            o.gro_sds = {num2cell([m, .3])};
            o.ind_sdlog = {num2cell([-100 -100])};
            o.ind_bs = {num2cell([0 0])};
            o.id_matches = {table([1; 2], [1; 2], 'VariableNames', {'id', 'id2'})};
            o.labels = {'cancellation'};
            o.conv = struct('data', {{{'name', 'n_eff', 'r_hat'; 'Intercept', 400, 1.00}}});
            rel = struct('analysis', 'trt_sserrvar', 'family', 'gamma', 'gammascale', 2, ...
                'engine', 'cmdstan', 'filename', 'cancel.psyrat', 'out', o);
            pd = struct('rel', rel, 'proc', struct('savename', 'cancel.psyrat'), 'ver', 'test');
            [u1, r1] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(r1.status, 'migrated');
            testCase.verifyEqual(u1.rel.out.gamma_resid{1}.mean_remainder, .10517091807564755, 'RelTol', 1e-12);
            testCase.verifyEqual(u1.rel.out.pop_sdlog, 100, 'AbsTol', 1e-12);
            testCase.verifyEqual(cell2mat(u1.rel.out.ind_sdlog{1}), (0.5 * log(.10517091807564755 + 1) - 100) * [1 1], 'AbsTol', 1e-12);
            [u2, r2] = psyrat_reprocess_gamma_result(u1);
            testCase.verifyEqual(r2.status, 'already_current');
            testCase.verifyEqual(u2, u1);
        end

        function testUnsupportedNumericRegimeRefused(testCase)
            % mu = 1e50 with all six observed SDs 1e-40: q ~ 1e-180 are finite
            % normal, but every positive remainder term (products of two or
            % three q) underflows before B = 1e100 restores it; the true R is a
            % normal number (~6e-260) that the product evaluation cannot produce
            pd = localStaticFromNormalizedRaw(1e50, 1e-40 * ones(1, 6));
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % subnormal stored SD; subnormal mean; overflowing mu^2
            pd = localStaticFromNormalizedRaw(1, [.3 .3 .3 .3 .3 1e-310]);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            pd = localStaticFromNormalizedRaw(1e-310, .3 * ones(1, 6));
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            pd = localStaticFromNormalizedRaw(1e200, .3 * ones(1, 6));
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % a subnormal encoding-decoded C is refused too
            pd = localStaticFromNormalizedRaw(1, .3 * ones(1, 6)); pd.rel.out.pop_sdlog = -360;
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % nonvacuous moderate limits stay accepted: small mean, small components, normal R
            pd = localStaticFromNormalizedRaw(1e-100, 1e-101 * ones(1, 6));
            [u, r] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(r.status, 'migrated');
            testCase.verifyGreaterThan(u.rel.out.gamma_resid{1}.mean_remainder, realmin);
            % and a complete zero facet at extreme magnitude remains a legitimate R = 0
            pd = localStaticFromNormalizedRaw(1e50, [1e-40 0 1e-40 0 1e-40 0]);
            [u, r] = psyrat_reprocess_gamma_result(pd);
            testCase.verifyEqual(r.status, 'migrated');
            testCase.verifyEqual(u.rel.out.gamma_resid{1}.mean_remainder, 0);
        end

        function testPositiveToZeroTransformRefused(testCase)
            % Astra xhigh case: mu = 1e100, sd_P = 1e-100, sd_OI = 1e150, the other
            % four SDs 0, C_typ = C_p = 1e-100. The source squares (1e-200, 1e300)
            % are normal, but the normalized q_P = (1e-100/1e100)^2 = 1e-400
            % underflows to 0; carrying that zero would make the remainder
            % preflight call every P-supported term structurally absent and
            % report R = 0 although the true R = 1e-100 (E = 2e-100). A positive
            % source mapping to zero must be refused with a precise message.
            pd = localStaticFromNormalizedRaw(1e100, [1e-100 0 0 0 0 1e150]);
            pd.rel.out.pop_sdlog = 0.5 * log(1e-100);
            try
                psyrat_reprocess_gamma_result(pd);
                testCase.verifyFail('Expected an unsupported_numeric_regime refusal.');
            catch err
                testCase.verifyEqual(err.identifier, 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
                testCase.verifyTrue(contains(err.message, '(sd_P/mu)^2'), 'Message must name the underflowing transform.');
                testCase.verifyTrue(contains(err.message, 'maps to zero'), 'Message must state positive-to-zero underflow.');
            end
            % a strictly positive pair deficit whose quotient is subnormal is refused
            % (a = 1, b = 0 -> divisor 2; h_PO = 3e-308 normal -> q_PO = 1.5e-308 < realmin)
            pd = localStaticFromNormalizedRaw(1, [1 0 0 sqrt(3e-308) 0 0]);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % a strictly positive source SD whose square underflows is refused
            pd = localStaticFromNormalizedRaw(1, [.3 .3 .3 .3 .3 1e-160]);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % a positive source whose square overflows is refused
            pd = localStaticFromNormalizedRaw(1, [.3 .3 .3 .3 .3 1e160]);
            testCase.verifyError(@() psyrat_reprocess_gamma_result(pd), 'psyrat_reprocess_gamma_result:unsupported_numeric_regime');
            % controls: an actual source zero stays a structural zero (zero facet, R = 0)
            st = localStratum(log(4), .5, 0, .4, 0, .2, 0, 3, 1:3, [1 2 3]);
            st.ell0 = zeros(st.nd, 1); st.w_p = zeros(st.nd, st.nsub);
            [u, r] = psyrat_reprocess_gamma_result(localStatic25V1('ls', localFinish('ls', st)));
            testCase.verifyEqual(r.status, 'migrated');
            testCase.verifyEqual(u.rel.out.gamma_resid{1}.mean_remainder, zeros(st.nd, 1));
            % and a tiny positive pair deficit (1e-16, inside the bound) is retained, not zeroed
            base = [.00002 .00002 .1 4e-10 .000002 .000002];
            tiny = base; tiny(4) = 4e-10 + 1e-16;
            [u, r] = psyrat_reprocess_gamma_result(localStaticFromNormalized(tiny));
            testCase.verifyEqual(r.status, 'migrated');
            Rzero = localRfromQ([.00002 .00002 .1 0 0 0]);
            testCase.verifyGreaterThan(u.rel.out.gamma_resid{1}.mean_remainder, Rzero * (1 + 1e-7), ...
                'A positive pair deficit inside the bound must be retained, not canonicalized to zero.');
            c = r.transforms(1).compatibility;
            testCase.verifyTrue(any(strcmp({c.unresolved.pair}, 'PO')), 'The 1e-16 PO deficit lies inside the bound and must be recorded.');
            poRow = c.unresolved(strcmp({c.unresolved.pair}, 'PO'));
            testCase.verifyGreaterThan(poRow(1).deficit, 0, 'The recorded PO deficit must be the retained positive value.');
        end

        function testHistoricalMarginalFlagIsNotATarget(testCase)
            for flag = [0 1]
                pd = localDynFixture(1, 2, 1, 1, 1);
                pd.rel.out.marginal = flag;
                [u, r] = psyrat_reprocess_gamma_result(pd);
                testCase.verifyEqual(u.rel.out.gamma_target.surface_target, 'population_marginal');
                testCase.verifyTrue(any(contains(r.notes, sprintf('marginal=%d', flag))));
                testCase.verifyEqual(u.rel.out.marginal, flag, 'Historical flag must be preserved, not rewritten.');
            end
        end

    end

    methods
        function d = scratchDir(testCase)
            d = fullfile(tempdir, ['psyrat_migration_' char(java.util.UUID.randomUUID().toString())]);
            mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
        end
    end
end

% =========================================================================
% Independent forward map: support-union products over all 63 subsets
% =========================================================================
function K = localSupportUnion(alpha, v)
% alpha [nd x 1]; v [1 x 6] = [vP vO vI vPO vPI vOI]. Returns [nd x 1] fields
% P, O, I, PO, PI, OI, R, B and mu = exp(alpha + V/2).
q = expm1(v);
B = exp(2 * alpha + sum(v));
acc = localUnionAccumulate(q);
K.B = B; K.mu = exp(alpha + sum(v) / 2);
K.P = B * acc(2); K.O = B * acc(3); K.PO = B * acc(4); K.I = B * acc(5);
K.PI = B * acc(6); K.OI = B * acc(7); K.R = B * acc(8);
end

function acc = localUnionAccumulate(q)
supp = logical([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]);
acc = zeros(1, 8);
for m = 1:63
    inT = logical(bitget(m, 1:6));
    u = any(supp(inT, :), 1);
    idx = 1 + u(1) + 2 * u(2) + 4 * u(3);
    acc(idx) = acc(idx) + prod(q(inT));
end
end

function R = localRfromQ(q)
% remainder from normalized q at B = 1 (mu = 1), independent enumeration
acc = localUnionAccumulate(q);
R = acc(8);
end

function st = localStratum(alpha0, sP, sO, sI, sPO, sPI, sOI, nsub, id2, ids)
st = struct();
st.nd = 4;
st.alpha = alpha0 + [0; 0.05; -0.05; 0.1];
st.v = [sP sO sI sPO sPI sOI].^2;
st.nsub = nsub;
st.id2 = id2(:);
st.ids = ids(:);
st.K = localSupportUnion(st.alpha, st.v);
st.s_scale = 0.3;
end

function st = localFinish(param, st)
if strcmp(param, 'ls')
    st.Ctyp = exp(2 * st.ell0);
    st.Cp   = exp(2 * (st.ell0 + st.w_p));
else
    Vw = sum(st.v(2:6));
    st.Ctyp = 2 * exp(2 * st.alpha + 2 * Vw - st.b_nu);
    st.Cp   = 2 * exp(2 * (st.alpha + st.u_p) + 2 * Vw - (st.b_nu + st.v_p));
end
end

function spec = localTwoStrataSpec()
% Two strata with distinct means, scales, counts and a permuted id2 in the second.
s1 = localStratum(log(5.5), .35, .22, .15, .12, .10, .08, 3, [1 2 3], [701 702 703]);
s2 = localStratum(log(7.0), .40, .18, .20, .09, .14, .06, 3, [3 1 2], [703 701 702]);
s1.ell0 = log(2.0) * ones(s1.nd, 1);   s1.w_p = [0.10 -0.25 0.05] .* ones(s1.nd, 1);
s2.ell0 = log(1.6) * ones(s2.nd, 1);   s2.w_p = [-0.15 0.30 0.00] .* ones(s2.nd, 1);
s1.b_nu = log(15) * ones(s1.nd, 1);    s1.u_p = [0.25 0.00 -0.30] .* ones(s1.nd, 1); s1.v_p = [-0.20 0.10 0.35] .* ones(s1.nd, 1);
s2.b_nu = log(22) * ones(s2.nd, 1);    s2.u_p = [-0.10 0.20 0.05] .* ones(s2.nd, 1); s2.v_p = [0.15 -0.05 -0.25] .* ones(s2.nd, 1);
spec = [localFinish('ls', s1), localFinish('ls', s2)];
specNu = [localFinish('lognu', s1), localFinish('lognu', s2)];
for s = 1:2
    spec(s).Ctyp_nu = specNu(s).Ctyp; spec(s).Cp_nu = specNu(s).Cp;
end
end

function pd = localStatic25V1(param, spec)
% Version-1 analysis-25 psyrat_data from the planted spec: observed-scale SDs
% from the independent forward map and the OBSERVATION-ONLY encoding.
nstrata = numel(spec);
o = struct();
o.mu = []; o.pop_sdlog = []; o.sig_occ = []; o.sig_trl = []; o.sig_trlxid = [];
o.sig_occxid = []; o.sig_trlxocc = [];
o.ind_bs = {}; o.gro_sds = {}; o.ind_sdlog = {}; o.id_matches = {}; o.conv = struct('data', {{}});
o.labels = cell(1, nstrata);
CtypLast = [];
for s = 1:nstrata
    st = spec(s);
    if strcmp(param, 'ls') || ~isfield(st, 'Ctyp_nu')
        Ctyp = st.Ctyp; Cp = st.Cp;
    else
        Ctyp = st.Ctyp_nu; Cp = st.Cp_nu;
    end
    CtypLast = Ctyp;
    pop = 0.5 * log(Ctyp);
    ind = 0.5 * log(Cp) - pop;
    o.mu(:, s)          = st.K.mu;
    o.pop_sdlog(:, s)   = pop;
    o.sig_occ(:, s)     = sqrt(st.K.O);
    o.sig_trl(:, s)     = sqrt(st.K.I);
    o.sig_trlxid(:, s)  = sqrt(st.K.PI);
    o.sig_occxid(:, s)  = sqrt(st.K.PO);
    o.sig_trlxocc(:, s) = sqrt(st.K.OI);
    o.gro_sds{s}   = num2cell([sqrt(st.K.P), st.s_scale * ones(st.nd, 1)]);
    o.ind_sdlog{s} = num2cell(ind);
    if isfield(st, 'u_p'), ib = st.u_p; else, ib = zeros(st.nd, st.nsub); end
    o.ind_bs{s}    = num2cell(ib);
    o.id_matches{s} = table(st.ids, st.id2, 'VariableNames', {'id', 'id2'});
    o.labels{s} = sprintf('stratum %d', s);
    o.conv.data{s} = {'name', 'n_eff', 'r_hat'; 'Intercept', 400, 1.00};
end
last = spec(nstrata);
if strcmp(param, 'ls')
    o.gamma_disp = struct('s_sd', last.s_scale * ones(last.nd, 1), 'rho_ps', zeros(last.nd, 1), ...
        'disp_factor', ones(last.nd, 1), 'alpha', last.alpha, ...
        's_o', sqrt(last.v(2)) * ones(last.nd, 1), 's_t', sqrt(last.v(3)) * ones(last.nd, 1), ...
        's_po', sqrt(last.v(4)) * ones(last.nd, 1), 's_pt', sqrt(last.v(5)) * ones(last.nd, 1), ...
        's_ot', sqrt(last.v(6)) * ones(last.nd, 1), 'log_sigma', last.ell0, ...
        'sigma_pot_e2_pop', CtypLast, 'is_ls', true);
    gammascale = 2;
else
    o.gamma_disp = struct('s_v', last.s_scale * ones(last.nd, 1), 'rho_pv', zeros(last.nd, 1), ...
        'disp_factor', ones(last.nd, 1), 'V_w', sum(last.v(2:6)) * ones(last.nd, 1), 'alpha', last.alpha, ...
        's_o', sqrt(last.v(2)) * ones(last.nd, 1), 's_t', sqrt(last.v(3)) * ones(last.nd, 1), ...
        's_po', sqrt(last.v(4)) * ones(last.nd, 1), 's_pt', sqrt(last.v(5)) * ones(last.nd, 1), ...
        's_ot', sqrt(last.v(6)) * ones(last.nd, 1), 'nu', exp(last.b_nu), ...
        'sigma_pot_e2_pop', CtypLast);
    gammascale = 1;
end
rel = struct('analysis', 'trt_sserrvar', 'family', 'gamma', 'gammascale', gammascale, ...
    'engine', 'cmdstan', 'filename', 'legacy_gamma25_fixture.psyrat', 'out', o);
pd = struct('rel', rel, 'proc', struct('savename', 'legacy_gamma25_fixture.psyrat'), 'ver', 'test');
end

function pd = localStaticFromNormalized(k)
% Single-draw legacy LS analysis-25 fixture with mu = 1 (B = 1), observed SDs
% sqrt(k) in order P O I PO PI OI, C_typ = 1 (pop_sdlog = 0), zero offsets;
% no gamma_disp (the frozen inverse cases have no raw provenance).
pd = localStaticFromNormalizedRaw(1, sqrt(k));
end

function pd = localStaticFromNormalizedRaw(mu, sds)
% Single-draw legacy LS analysis-25 fixture from a raw observed mean and the
% six raw observed SDs (order P O I PO PI OI); C_typ = 1, zero offsets.
o = struct();
o.mu = mu; o.pop_sdlog = 0;
o.sig_occ = sds(2); o.sig_trl = sds(3); o.sig_occxid = sds(4);
o.sig_trlxid = sds(5); o.sig_trlxocc = sds(6);
o.gro_sds = {num2cell([sds(1), .3])};
o.ind_sdlog = {num2cell([0 0])};
o.ind_bs = {num2cell([0 0])};
o.id_matches = {table([11; 12], [2; 1], 'VariableNames', {'id', 'id2'})};
o.labels = {'frozen'};
o.conv = struct('data', {{{'name', 'n_eff', 'r_hat'; 'Intercept', 400, 1.00}}});
rel = struct('analysis', 'trt_sserrvar', 'family', 'gamma', 'gammascale', 2, ...
    'engine', 'cmdstan', 'filename', 'frozen_inverse_case.psyrat', 'out', o);
pd = struct('rel', rel, 'proc', struct('savename', 'frozen_inverse_case.psyrat'), 'ver', 'test');
end

function pd = localDynFixture(scale, facets, subjects, ndim, nstrata)
% Historical (unstamped) dynamic gamma raw-draw result, mirroring the task
% fixture of the 2C review (local evidence) without adding it to the path.
if facets == 1
    pd = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
else
    pd = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffTrtSummaryData();
end
o = pd.rel.out; nd = size(o.mu, 1);
if scale == 2
    o.pop_sdlog = repmat(log(1.4), nd, 1);
    o.b_sigma = o.b_nu;
    o.ind_sdlog = o.ind_nu;
    o = rmfield(o, {'pop_lognu', 'b_nu', 'ind_nu'});
    if facets == 2, o = rmfield(o, 'ind_bs'); end
end
if ndim == 2
    z2 = [-.7; .8; .2];
    o.dimz{1} = [o.dimz{1} z2];
    o.ssinfo{1}.z2 = z2;
    o.b{1} = num2cell([cell2mat(o.b{1}) repmat([.07 -.05], nd, 1)]);
    if scale == 2, fld = 'b_sigma'; else, fld = 'b_nu'; end
    o.(fld){1} = num2cell([cell2mat(o.(fld){1}) repmat([.1 .04], nd, 1)]);
end
if ~subjects
    if facets == 1, pd.rel.analysis = 'ic_dynrel'; else, pd.rel.analysis = 'ic_dynrel_trt'; end
    for fld = {'ssinfo', 'ind_bs', 'ind_nu', 'ind_sdlog'}
        if isfield(o, fld{1}), o = rmfield(o, fld{1}); end
    end
end
if nstrata == 2
    names = fieldnames(o);
    for i = 1:numel(names)
        x = o.(names{i});
        if iscell(x) || isnumeric(x), o.(names{i}) = [x x]; end
    end
    o.labels = {'stratum A', 'stratum B'};
    o.mu(:, 2) = o.mu(:, 2) + log(1.3);
    o.ntrials(2) = 11;
    if subjects
        o.ssinfo{2} = o.ssinfo{1}([3 1 2], :);
        o.ssinfo{2}.id = [703; 701; 702];
        o.ssinfo{2}.trls = [7; 19; 12];
    end
end
pd.rel.out = o; pd.rel.gammascale = scale; pd.rel.ndim = ndim;
pd.rel.dim_names = {'Dimension 1'};
if ndim == 2, pd.rel.dim_names{2} = 'Dimension 2'; end
end

function m = localCellToNumeric(x)
% numeric twin of a num2cell block (2-D or [n x 2 x 2]), for storage-contract tests
if ndims(x) == 3
    m = zeros(size(x));
    for i = 1:size(x, 2)
        for j = 1:size(x, 3)
            m(:, i, j) = cell2mat(x(:, i, j));
        end
    end
else
    m = cell2mat(x);
end
end

function b = localReadBytes(f)
fid = fopen(f, 'r'); c = onCleanup(@() fclose(fid));
b = fread(fid, Inf, '*uint8');
end
