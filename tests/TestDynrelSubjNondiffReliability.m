classdef TestDynrelSubjNondiffReliability < PsyRATTestBase
    % Deterministic unit tests for the NON-difference SUBJECT-LEVEL dynamic/
    % conditional reliability variants (Rast & Clayson): the one-facet
    % ic_dynrel_sserrvar (analysis 26) and the trial + occasion two-facet
    % ic_dynrel_sserrvar_trt (analysis 27). These variants reuse the same
    % location-scale model as ic_dynrel / ic_dynrel_trt (the model already
    % estimates a per-subject residual scale, ind_sdlog) but ADD a per-participant
    % table reporting each subject's conditional reliability phi_s(z) at their own
    % standardized dimension value, while still reporting the typical-person
    % population surface as a reference.
    %
    % No CmdStan/MatlabStan required: the fixtures carry constant posterior draws
    % (PsyRATTestDataFactory.makeDynrelSubjNondiff*SummaryData) so the point
    % estimates are exact. The one-facet per-participant coefficients are checked
    % against the closed-form subject-level formula; the two-facet coefficients
    % are checked against a direct psyrat_rel_trt reconstruction. The live
    % known-truth CmdStan recovery is a separate (slow, gated) concern.

    methods (Test)

        % ---- one-facet (ic_dynrel_sserrvar, analysis 26) --------------------

        function testOneFacetProducesSurfaceAndSsrel(testCase)
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);
            % the typical-person population surface is still reported
            testCase.verifyNotEmpty(st.G.pt);
            testCase.verifyNotEmpty(st.D.pt);
            % and the per-participant table is present with the expected columns
            T = st.ssrel_table;
            testCase.verifyNotEmpty(T);
            expcols = {'id','id2','z1','trls','gen_pt','gen_ll','gen_ul',...
                'dep_pt','dep_ll','dep_ul','iccg_pt','iccd_pt',...
                'sem_gen_pt','sem_dep_pt'};
            testCase.verifyTrue(all(ismember(expcols, T.Properties.VariableNames)));
            testCase.verifyEqual(height(T), 3);
        end

        function testOneFacetPerSubjectMatchesClosedForm(testCase)
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            o = data.rel.out;
            % pull the (constant) generative inputs from the fixture
            gro   = cell2mat(o.gro_sds{1}); sig_p = gro(1,1);
            pop   = o.pop_sdlog(1);
            bsig  = cell2mat(o.b_sigma{1}); bsig = bsig(1);   % KDIM = 1
            sigi  = o.sig_trl(1);
            indsd = cell2mat(o.ind_sdlog{1}); indsd = indsd(1,:); % 1 x NSUB
            ss    = o.ssinfo{1};

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;

            for r = 1:height(ss)
                col = ss.id2(r); zr = ss.z1(r); nr = ss.trls(r);
                % per-subject relative residual variance at this subject's own z
                rv = exp(pop + bsig*zr + indsd(col))^2;
                % generalizability excludes sigma_i^2; dependability includes it
                expG   = sig_p^2 / (sig_p^2 + rv/nr);
                expPhi = sig_p^2 / (sig_p^2 + (rv + sigi^2)/nr);
                expICCg = sig_p^2 / (sig_p^2 + rv);
                expICCd = sig_p^2 / (sig_p^2 + rv + sigi^2);
                expSEMg = sqrt(rv/nr);
                expSEMd = sqrt((rv + sigi^2)/nr);
                testCase.verifyEqual(T.gen_pt(r), expG,   'AbsTol', 1e-10);
                testCase.verifyEqual(T.dep_pt(r), expPhi, 'AbsTol', 1e-10);
                testCase.verifyEqual(T.iccg_pt(r), expICCg, 'AbsTol', 1e-10);
                testCase.verifyEqual(T.iccd_pt(r), expICCd, 'AbsTol', 1e-10);
                testCase.verifyEqual(T.sem_gen_pt(r), expSEMg, 'AbsTol', 1e-10);
                testCase.verifyEqual(T.sem_dep_pt(r), expSEMd, 'AbsTol', 1e-10);
                % constant draws: the credible interval collapses to the point
                testCase.verifyEqual(T.gen_ll(r), T.gen_pt(r), 'AbsTol', 1e-10);
                testCase.verifyEqual(T.gen_ul(r), T.gen_pt(r), 'AbsTol', 1e-10);
                % dependability never exceeds generalizability
                testCase.verifyLessThanOrEqual(T.dep_pt(r), T.gen_pt(r) + 1e-12);
            end
        end

        function testOneFacetMarginalEstimandAllowed(testCase)
            % the population-average (marginal) estimand IS defined for the
            % non-difference variants; requesting it must not warn or force it off
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            summ = testCase.verifyWarningFree(@() ...
                psyrat_dynrel_summary(data.rel, 'CI', 0.95, 'marginal', 1));
            testCase.verifyEqual(summ.marginal, 1);
        end

        % ---- two-facet (ic_dynrel_sserrvar_trt, analysis 27) ----------------

        function testTwoFacetProducesSurfaceAndSsrel(testCase)
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffTrtSummaryData();
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);
            testCase.verifyNotEmpty(st.G.pt);
            testCase.verifyNotEmpty(st.D.pt);
            T = st.ssrel_table;
            testCase.verifyNotEmpty(T);
            % the two-facet per-participant table also carries the occasion count
            expcols = {'id','id2','z1','trls','occs','gen_pt','dep_pt',...
                'iccg_pt','iccd_pt','sem_gen_pt','sem_dep_pt'};
            testCase.verifyTrue(all(ismember(expcols, T.Properties.VariableNames)));
            testCase.verifyEqual(height(T), 3);

            % RC-06. The generalizability and dependability ICC columns must
            % DIFFER. They were identical while the subject-level ICC routed
            % through psyrat_icc, which sees neither gcoeff nor the facet-only
            % components. Non-vacuous here: the fixture plants sig_occ = 1.5,
            % sig_trl = 2 and sig_trlxocc = 0.5, all nonzero, so absolute error
            % strictly exceeds relative error for every reltype.
            testCase.verifyNotEqual(T.iccd_pt, T.iccg_pt);
            testCase.verifyTrue(all(T.iccd_pt < T.iccg_pt), ...
                'absolute-error ICC must be strictly below the relative-error ICC');
        end

        function testTwoFacetPerSubjectMatchesRelTrt(testCase)
            % verify the summary wires the z-conditioned per-subject residual and
            % the population facet SDs into psyrat_rel_trt correctly: reconstruct
            % each subject's coefficient directly and compare.
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffTrtSummaryData();
            o = data.rel.out;
            pop   = o.pop_sdlog;                     % nd x 1
            gro   = cell2mat(o.gro_sds{1});
            sig_p = gro(:,1);                         % nd x 1 (constant)
            bsig  = cell2mat(o.b_sigma{1});           % nd x 1 (KDIM = 1)
            indsd = cell2mat(o.ind_sdlog{1});         % nd x NSUB
            ss    = o.ssinfo{1};
            nocc  = 1; % local_resolve_nocc default when 'nocc' unset

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;

            for r = 1:height(ss)
                col = ss.id2(r); zr = ss.z1(r); nr = ss.trls(r);
                err_r = exp(pop + bsig*zr + indsd(:,col)); % nd x 1 residual SD at z_s
                [~, Gpt, ~] = psyrat_rel_trt('gcoeff', 2, 'reltype', 3, ...
                    'bp', sig_p, 'bo', o.sig_occ, 'bt', o.sig_trl, ...
                    'txp', o.sig_trlxid, 'oxp', o.sig_occxid, 'txo', o.sig_trlxocc, ...
                    'err', err_r, 'obs', nr, 'nocc', nocc, 'CI', 0.95);
                [~, Dpt, ~] = psyrat_rel_trt('gcoeff', 1, 'reltype', 3, ...
                    'bp', sig_p, 'bo', o.sig_occ, 'bt', o.sig_trl, ...
                    'txp', o.sig_trlxid, 'oxp', o.sig_occxid, 'txo', o.sig_trlxocc, ...
                    'err', err_r, 'obs', nr, 'nocc', nocc, 'CI', 0.95);
                testCase.verifyEqual(T.gen_pt(r), Gpt, 'AbsTol', 1e-9);
                testCase.verifyEqual(T.dep_pt(r), Dpt, 'AbsTol', 1e-9);
                testCase.verifyLessThanOrEqual(T.dep_pt(r), T.gen_pt(r) + 1e-12);
            end
        end

        % ---- GAMMA siblings (family = 'gamma') -----------------------------
        %
        % These pin the psyrat_dynrel_summary WIRING for the gamma subject-level
        % designs: that the isgamma branch now produces an ssrel_table at all (it
        % returned [] before, which surfaced as a dead "Show per-participant
        % table" button in the viewer), that the table carries the shared schema,
        % and that its numbers come from the conditional calculators. The
        % calculators' own mathematics is proven separately and in more depth by
        % TestGammaDynrelSserrReliability; the live fit -> extraction plumbing by
        % tests/integration/TestGammaDynrelSserrRecovery.

        function testGammaOneFacetProducesSurfaceAndSsrel(testCase)
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);
            % the population surface is still reported (estimand #2, marginal)
            testCase.verifyNotEmpty(st.G.pt);
            testCase.verifyNotEmpty(st.D.pt);
            % and the per-participant CONDITIONAL table is now present
            T = st.ssrel_table;
            testCase.verifyNotEmpty(T, ...
                ['The gamma subject-level branch produced no per-participant ', ...
                'table. Analysis 26 has silently degraded to analysis 11.']);
            expcols = {'id','id2','z1','trls','gen_pt','gen_ll','gen_ul',...
                'dep_pt','dep_ll','dep_ul','iccg_pt','iccd_pt',...
                'sem_gen_pt','sem_dep_pt'};
            testCase.verifyTrue(all(ismember(expcols, T.Properties.VariableNames)));
            testCase.verifyEqual(height(T), 3);
            testCase.verifyTrue(all(isfinite(T.gen_pt)) && all(isfinite(T.dep_pt)));
            % coefficients are genuinely per-participant, not one repeated value.
            % max-min, not range(): range is Statistics Toolbox, and CI installs
            % base MATLAB only (see tests/TestBaseMatlabOnly.m).
            testCase.verifyGreaterThan(max(T.gen_pt) - min(T.gen_pt), 1e-6);
        end

        function testGammaOneFacetPerSubjectMatchesCalculator(testCase)
            % The summary must pass exactly the draws the calculator expects; any
            % mis-wiring (wrong slot, wrong stratum unwrap, transposed slopes)
            % changes the numbers. Reconstruct by calling the calculator directly.
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            o = data.rel.out;
            gro = cell2mat(o.gro_sds{1});

            [gen, dep] = psyrat_ssrel_dynrel_gamma( ...
                'alpha0', o.mu, 'b', cell2mat(o.b{1}), ...
                'lognu0', o.pop_lognu, 'b_nu', cell2mat(o.b_nu{1}), ...
                'sig_p', gro(:,1), 'sig_i', o.sig_trl, ...
                'ind_bs', cell2mat(o.ind_bs{1}), 'ind_nu', cell2mat(o.ind_nu{1}), ...
                'idtable', o.ssinfo{1}, 'ndim', 1, 'CI', 0.95);

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;

            testCase.verifyEqual(T.gen_pt, gen.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.dep_pt, dep.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.iccg_pt, gen.icc_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.iccd_pt, dep.icc_pt, 'AbsTol', 1e-12);
            for r = 1:height(T)
                testCase.verifyLessThanOrEqual(T.dep_pt(r), T.gen_pt(r) + 1e-12);
            end
        end

        function testGammaTwoFacetProducesSurfaceAndSsrel(testCase)
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffTrtSummaryData();
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);
            testCase.verifyNotEmpty(st.G.pt);
            T = st.ssrel_table;
            testCase.verifyNotEmpty(T, ...
                ['The gamma two-facet subject-level branch produced no ', ...
                'per-participant table. Analysis 27 has degraded to analysis 14.']);
            expcols = {'id','id2','z1','trls','occs','gen_pt','dep_pt', ...
                'iccg_pt','iccd_pt','sem_gen_pt','sem_dep_pt'};
            testCase.verifyTrue(all(ismember(expcols, T.Properties.VariableNames)));
            testCase.verifyEqual(height(T), 3);
            % coefficients are genuinely per-participant, not one repeated value
            testCase.verifyGreaterThan(max(T.gen_pt) - min(T.gen_pt), 1e-6);
        end

        function testGammaTwoFacetPerSubjectMatchesCalculator(testCase)
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffTrtSummaryData();
            o = data.rel.out;
            gro = cell2mat(o.gro_sds{1});
            % local_resolve_nocc default when 'nocc' is unset is 1 (the owner
            % default), NOT REL.out.nocc - that is only used for 'observed'. Same
            % assumption the Gaussian two-facet test above makes.
            nocc = 1;

            [gen, dep] = psyrat_ssrel_dynrel_trt_gamma( ...
                'alpha0', o.mu, 'b', cell2mat(o.b{1}), ...
                'lognu0', o.pop_lognu, 'b_nu', cell2mat(o.b_nu{1}), ...
                'sig_p', gro(:,1), 'sig_occ', o.sig_occ, 'sig_trl', o.sig_trl, ...
                'sig_occxid', o.sig_occxid, 'sig_trlxid', o.sig_trlxid, ...
                'sig_trlxocc', o.sig_trlxocc, ...
                'ind_bs', cell2mat(o.ind_bs{1}), 'ind_nu', cell2mat(o.ind_nu{1}), ...
                'idtable', o.ssinfo{1}, 'ndim', 1, 'reltype', 3, ...
                'nocc', nocc, 'CI', 0.95);

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;

            testCase.verifyEqual(T.gen_pt, gen.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.dep_pt, dep.dep_pt, 'AbsTol', 1e-12);
            % RC-06 regression guard, gamma side: the absolute-error ICC must not
            % exceed the relative-error one.
            for r = 1:height(T)
                testCase.verifyLessThanOrEqual(T.iccd_pt(r), T.iccg_pt(r) + 1e-12);
            end
        end

    end
end
