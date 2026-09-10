classdef TestAnalysisFamily < PsyRATTestBase
    %Tests for psyrat_analysis_family, the single-source-of-truth classifier that
    %maps a REL.analysis string to its reliability family (sing/trt/dod/splits/
    %dynrel). Pure lookup: no CmdStan, no figures. These pin the exact mapping the
    %GUI router, the dynrel/splits summary guards, and the headless report
    %dispatch all rely on, so a future variant cannot silently change one site's
    %behavior without this test catching it.

    methods (Test)

        function testSingFamily(testCase)
            % One-facet (single-session) strings, including the legacy default.
            for a = {'ic','ic_sserrvar','ic_diff','ic_diff_sserrvar'}
                testCase.verifyEqual(psyrat_analysis_family(a{1}), 'sing', ...
                    sprintf('%s should map to sing', a{1}));
            end
        end

        function testTrtFamily(testCase)
            for a = {'trt','trt_sserrvar','trt_diff'}
                testCase.verifyEqual(psyrat_analysis_family(a{1}), 'trt', ...
                    sprintf('%s should map to trt', a{1}));
            end
        end

        function testDodFamily(testCase)
            testCase.verifyEqual(psyrat_analysis_family('ic_dodiff'), 'dod');
        end

        function testSplitsFamily(testCase)
            for a = {'ic_splits','trt_splits'}
                testCase.verifyEqual(psyrat_analysis_family(a{1}), 'splits', ...
                    sprintf('%s should map to splits', a{1}));
            end
        end

        function testDynrelFamily(testCase)
            % The canonical 12-string dynamic-reliability set (matches the
            % accepted set in psyrat_dynrel_summary and the GUI router).
            dynrel = {'ic_dynrel','ic_dynrel_sserrvar',...
                'ic_diff_dynrel','ic_diff_dynrel_sserrvar',...
                'ic_dynrel_trt','ic_dynrel_sserrvar_trt','ic_diff_dynrel_trt',...
                'ic_diff_dynrel_sserrvar_trt',...
                'ic_diff_dynrel_rescor','ic_diff_dynrel_trt_rescor',...
                'ic_diff_dynrel_sserrvar_rescor',...
                'ic_diff_dynrel_sserrvar_trt_rescor',...
                'ic_diff_dynrel_sserrvar_rho','ic_diff_dynrel_sserrvar_trt_rho'};
            for a = dynrel
                testCase.verifyEqual(psyrat_analysis_family(a{1}), 'dynrel', ...
                    sprintf('%s should map to dynrel', a{1}));
            end
        end

        function testEmptyAndMissingDefaultToSing(testCase)
            % Missing/empty -> the legacy one-facet default ('ic') -> sing, so the
            % GUI router opens the single-session viewer for pre-analysis-field
            % structs, matching its historical behavior.
            testCase.verifyEqual(psyrat_analysis_family(''), 'sing');
            testCase.verifyEqual(psyrat_analysis_family(), 'sing');
        end

        function testUnknownReturnsEmpty(testCase)
            % An unrecognized non-empty string returns '' so each caller keeps its
            % own fallback (router opens no viewer; guards/report error clearly).
            testCase.verifyEqual(psyrat_analysis_family('not_a_real_analysis'), '');
            testCase.verifyEqual(psyrat_analysis_family('ic_dynrel_nonexistent'), '');
        end

        function testAcceptsStringScalar(testCase)
            % string-typed input (not just char) is accepted and char-normalized.
            testCase.verifyEqual(psyrat_analysis_family("trt_diff"), 'trt');
            testCase.verifyEqual(psyrat_analysis_family("ic_splits"), 'splits');
        end

    end
end
