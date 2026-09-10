classdef TestGammaLsDynrelSummaryDispatch < PsyRATTestBase
    %TESTGAMMALSDYNRELSUMMARYDISPATCH Pins psyrat_dynrel_summary's routing
    % between the two gamma scale parameterizations, and the labels of the
    % variance-component table it builds for each.
    %
    % WHY ROUTING NEEDS ITS OWN TEST. gammascale = 1 and gammascale = 2 are
    % different estimands that share the same family name, most of the same
    % REL.out slot names, and - critically - the same gro_sds column layout with
    % a DIFFERENT quantity in column 2 (s_v under 1, s_sd under 2). Routing a
    % location-scale fit into the log-nu branch does not error. It reads
    % pop_lognu (absent) or, worse where the slots do line up, reports s_sd under
    % the label s_v. The failure mode is a wrong answer with a plausible face, so
    % the routing has to be asserted rather than assumed.
    %
    % The calculators' own mathematics is proven separately by
    % TestGammaDynrelLsReliability and TestGammaDynrelLsSserrReliability; this
    % file only checks that the summary hands them the right draws and labels the
    % output for the right parameterization.

    methods (Test)

        function testLocationScaleRunRoutesToTheLocationScaleCalculators(testCase)
            % End-to-end through the summary: surface present, per-participant
            % table present, and both numerically identical to calling the
            % location-scale calculators directly on the same draws. Any
            % mis-wiring - wrong slot, wrong stratum unwrap, transposed slopes -
            % moves the numbers.
            data = PsyRATTestDataFactory.makeGammaLsDynrelSubjNondiffSummaryData();
            o = data.rel.out;
            gro = cell2mat(o.gro_sds{1});

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);

            testCase.assertNotEmpty(st.G.pt, 'No population surface was produced.');
            T = st.ssrel_table;
            testCase.assertNotEmpty(T, ...
                ['The location-scale gamma subject-level branch produced no ', ...
                'per-participant table. Analysis 26 has degraded to 11.']);
            testCase.verifyEqual(height(T), 3);

            [gen, dep] = psyrat_ssrel_dynrel_gamma_ls( ...
                'alpha0', o.mu, 'b', cell2mat(o.b{1}), ...
                'logsig0', o.pop_sdlog, 'b_sigma', cell2mat(o.b_sigma{1}), ...
                'sig_p', gro(:,1), 'sig_i', o.sig_trl, ...
                'ind_sd', cell2mat(o.ind_sdlog{1}), ...
                'idtable', o.ssinfo{1}, 'ndim', 1, 'CI', 0.95);

            testCase.verifyEqual(T.gen_pt,  gen.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.dep_pt,  dep.dep_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.iccg_pt, gen.icc_pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(T.iccd_pt, dep.icc_pt, 'AbsTol', 1e-12);

            % coefficients are genuinely per-participant, not one repeated value.
            % max-min, not range(): range is Statistics Toolbox and CI installs
            % base MATLAB only (tests/TestBaseMatlabOnly.m).
            testCase.verifyGreaterThan(max(T.gen_pt) - min(T.gen_pt), 1e-6);
        end

        function testSurfaceMatchesTheLocationScaleSurfaceCalculator(testCase)
            % The same check for the population surface, so a routing bug that
            % only affected the grid path could not hide behind a correct
            % per-participant table.
            data = PsyRATTestDataFactory.makeGammaLsDynrelSubjNondiffSummaryData();
            o = data.rel.out;
            gro = cell2mat(o.gro_sds{1});

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            st = summ.strata(1);

            rel = psyrat_rel_dynrel_gamma_ls('alpha0', o.mu, ...
                'b', cell2mat(o.b{1}), 'logsig0', o.pop_sdlog, ...
                'b_sigma', cell2mat(o.b_sigma{1}), 'sig_p', gro(:,1), ...
                'sig_sd', gro(:,2), 'sig_i', o.sig_trl, 'ndim', 1, ...
                'z1', st.z1, 'obs', st.obs, 'CI', 0.95);

            testCase.verifyEqual(st.G.pt, rel.G.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(st.D.pt, rel.D.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(st.sige_pt, rel.sige_pt, 'AbsTol', 1e-12);
        end

        function testVarcompTableNamesTheScaleSubmodelCorrectly(testCase)
            % THE SILENT-WRONG-ANSWER GUARD. The log-nu builder emits rows named
            % beta(0), nu(0), s_v, rho_pv and b_nu_k. None of those exist under
            % gammascale = 2, but every one of them would still PRINT, because
            % the underlying draws occupy the same slots. Assert the
            % location-scale names are there and the log-nu names are not.
            data = PsyRATTestDataFactory.makeGammaLsDynrelSubjNondiffSummaryData();
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            sym = string(summ.strata(1).varcomp{:,2});

            for want = ["log sigma(0)", "sigma(0)", "s_sd", "rho_p,sigma", "b_sigma_1"]
                testCase.verifyTrue(any(sym == want), ...
                    sprintf('Location-scale component "%s" is missing.', want));
            end
            for unwanted = ["beta(0)", "nu(0)", "s_v", "rho_pv", "b_nu_1"]
                testCase.verifyFalse(any(sym == unwanted), ...
                    sprintf(['Component "%s" belongs to the log-nu ' ...
                    'parameterization and must not appear on a gammascale=2 ' ...
                    'run - the draws in that slot are a different quantity.'], ...
                    unwanted));
            end

            % s_sd must carry the value from gro_sds column 2, not the trial SD
            % from the neighbouring read. The fixture keeps them distinct (0.27
            % vs 0.30) precisely so this can be checked.
            row = summ.strata(1).varcomp(sym == "s_sd", :);
            testCase.verifyEqual(row{1,3}, 0.27, 'AbsTol', 1e-12, ...
                's_sd was not sourced from gro_sds column 2.');
        end

        function testLogNuRunIsUnchangedByTheNewBranch(testCase)
            % REGRESSION. The gammascale = 1 fixture carries no gammascale field
            % at all, which is also how every .psyrat saved before the option
            % existed arrives. It must still route to the log-nu branch and still
            % emit the log-nu component names.
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            testCase.assertFalse(isfield(data.rel, 'gammascale'), ...
                'Fixture drifted: the log-nu fixture must carry no gammascale.');

            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            sym = string(summ.strata(1).varcomp{:,2});
            testCase.verifyTrue(any(sym == "s_v"), ...
                'A log-nu run must still report s_v.');
            testCase.verifyTrue(any(sym == "nu(0)"), ...
                'A log-nu run must still report nu(0).');
            testCase.verifyFalse(any(sym == "s_sd"), ...
                'A log-nu run must not report the location-scale s_sd.');
            testCase.verifyFalse(summ.isgamma_ls, ...
                'summary.isgamma_ls must be false for a log-nu run.');
            testCase.verifyTrue(summ.isgamma);
        end

        function testGammascaleOneIsRoutedAsLogNuEvenWhenStamped(testCase)
            % psyrat_computevarcomp stamps REL.gammascale on EVERY run and
            % defaults it to 1, so the common case is a log-nu run that DOES
            % carry the field. Gating on the field's presence rather than its
            % value would send every one of them down the wrong branch.
            data = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            data.rel.gammascale = 1;
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            testCase.verifyFalse(summ.isgamma_ls);
            testCase.verifyTrue(any(string(summ.strata(1).varcomp{:,2}) == "s_v"));
        end

        function testReportEstimandNoteNamesTheRightConditioningSet(testCase)
            % The headless report's estimand note tells the reader what the
            % per-participant table is conditional on. Under gammascale = 2 that
            % is each person's own log RESIDUAL SD; their own log-mean does not
            % enter the coefficient at all. Carrying the log-nu wording over
            % would be wrong in both halves, and it is prose, so no numeric gate
            % would catch it.
            data = PsyRATTestDataFactory.makeGammaLsDynrelSubjNondiffSummaryData();
            report = psyrat_report(data);
            note = string(report.estimand_note);

            testCase.verifyTrue(contains(note, "log residual SD"), ...
                'The location-scale note must name the log residual SD.');
            testCase.verifyFalse(contains(note, "log-nu"), ...
                'The location-scale note must not mention log-nu.');
            testCase.verifyFalse(contains(note, "delta_p"), ...
                'The Gaussian typical-person label must not appear on a gamma run.');

            % the log-nu run keeps its own wording
            d1 = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            note1 = string(psyrat_report(d1).estimand_note);
            testCase.verifyTrue(contains(note1, "log-nu"), ...
                'A log-nu run must still be described as conditional on log-nu.');
        end

        function testProvenanceReportsTheScaleSubmodelItActuallyFit(testCase)
            % psyrat_provenance_lines pools gro_sds column 2 for its dispersion
            % header. Under gammascale = 2 that column is s_sd, so the log-nu
            % lines would report a different quantity under the label s_v - a
            % silently wrong header rather than a crash. Two things are pinned:
            % that the log-nu lines are gone, and that a location-scale summary
            % appears in their place.
            data = PsyRATTestDataFactory.makeGammaLsDynrelSubjNondiffSummaryData();
            L = string(psyrat_provenance_lines(data));

            testCase.verifyFalse(any(contains(L, "s_v")), ...
                'A gammascale=2 run must not report a person log-nu SD.');
            testCase.verifyFalse(any(contains(L, "2*mu(z)^2/nu(z)")), ...
                ['The mean-coupled-residual line describes the log-nu model ', ...
                'and is false here.']);
            testCase.verifyFalse(any(contains(L, "cancels in the coefficient")), ...
                ['The mean-slope cancellation does NOT hold under the ', ...
                'location-scale parameterization.']);

            testCase.verifyTrue(any(contains(L, "LOCATION-SCALE")), ...
                'The scale parameterization must be stated on the run.');
            testCase.verifyTrue(any(contains(L, "s_sd): median 0.270")), ...
                's_sd must be reported, sourced from gro_sds column 2.');
            testCase.verifyTrue(any(contains(L, "not comparable")), ...
                ['The header must warn that the two parameterizations are ', ...
                'different estimands.']);

            % the log-nu run is untouched
            d1 = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            L1 = string(psyrat_provenance_lines(d1));
            testCase.verifyTrue(any(contains(L1, "person log-nu SD, s_v")), ...
                'A log-nu run must still report s_v.');
            testCase.verifyFalse(any(contains(L1, "LOCATION-SCALE")), ...
                'A log-nu run must not claim the location-scale parameterization.');
        end

        function testGaussianDynrelGetsNoGammaProvenanceLines(testCase)
            % The location-scale branch is keyed on pop_sdlog + b_sigma, which the
            % GAUSSIAN dynrel also writes. It is safe only because the caller
            % gates on family = 'gamma'. Pin that gate: a Gaussian run must emit
            % no gamma dispersion lines at all.
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            L = string(psyrat_provenance_lines(data));
            testCase.verifyFalse(any(contains(L, "LOCATION-SCALE")));
            testCase.verifyFalse(any(contains(L, "Gamma dispersion")));
        end

        function testGaussianRunNeverPicksUpTheLocationScaleFlag(testCase)
            % The other half of the same guard: gammascale is stamped on Gaussian
            % runs too. isgamma_ls must be gated on the family as well, or a
            % Gaussian analysis would inherit gamma labeling.
            data = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            data.rel.gammascale = 2;   % as a stamped-on-every-run field would be
            summ = psyrat_dynrel_summary(data.rel, 'CI', 0.95);
            testCase.verifyFalse(summ.isgamma, ...
                'A Gaussian run must not be flagged as gamma.');
            testCase.verifyFalse(summ.isgamma_ls, ...
                ['A Gaussian run must not be flagged as location-scale gamma ', ...
                'just because the gammascale field is present.']);
        end

    end
end
