classdef TestViewerEditFieldFallback < PsyRATTestBase
    % Regression guard: complex-literal text in the data-splits viewer's
    % "Splits n'" field must fall back to the field default, exactly like any
    % other invalid input.
    %
    % local_parse_pos (duplicated in psyrat_startview_splits and
    % psyrat_startview_dynrel) documented "return default on empty/invalid" but
    % accepted complex literals: str2double('3-4i') returns 3-4i, isfinite is
    % true for it, and MATLAB's ordering comparison (x > 0) uses only the REAL
    % part, so the 'isfinite(x) && x > 0' accept condition passed. Nothing
    % downstream caught it either -- psyrat_splits_summary's obs screen is
    % isnumeric+isscalar (a complex scalar passes both) -- so the complex n'
    % reached psyrat_rel_sing and errored in quantile, which delegates to
    % prctile: "Invalid data type. First argument must be an array of real
    % values." The user-visible symptom was an unexplained crash instead of the
    % documented fallback. Fixed by adding isreal to the accept condition.
    %
    % The dynamic-reliability viewer's copy of the same helper is pinned by
    % TestExportCallbacks/testDynrelVarcompSaveRejectsComplexNPrime, which
    % asserts on that viewer's resolved-n' export header (the observable there),
    % alongside its existing default/override cases.
    %
    % ASSERTING ON THE COEFFICIENTS TABLE IS DELIBERATE. The splits viewer's
    % VARIANCE-COMPONENT table is invariant to n' (measured: blank vs n_s = 9
    % renders byte-identical), so an equality check on it cannot distinguish
    % "the parse fell back to the default" from "the parse accepted the value"
    % -- it would pin only the absence of the crash. The coefficients table does
    % move with n_s, and testValidOverrideMovesTheCoefficients below asserts
    % that it does, so the fallback check that follows is discriminating rather
    % than vacuously true.
    %
    % Strategy (as in TestSingViewerCoefficientControl): `matlab -batch` builds
    % real figures, so the viewer is opened on a synthetic factory fixture and
    % its real button callbacks are fired.

    methods (Test)

        function testValidOverrideMovesTheCoefficients(testCase)
            % Non-vacuity guard for the test below: prove the observable it
            % compares actually responds to the parsed n'. The fixture is 3
            % persons x 4 splits, so its observed median n_s is 4; 9 is chosen
            % to differ from that AND from any value a fallback could produce.
            base = localSplitsCoefficients(testCase, '');
            moved = localSplitsCoefficients(testCase, '9');
            testCase.verifyNotEqual(moved, base, ...
                ['the coefficients table must respond to the splits n'' ' ...
                'override, otherwise the fallback assertion below is vacuous']);
        end

        function testComplexSplitsTextFallsBackToTheDefault(testCase)
            % The regression itself. Pre-fix this errored in quantile/prctile
            % rather than returning a table, so the assertion doubles as a
            % crash guard; post-fix the parse must yield the same default the
            % blank field does.
            base = localSplitsCoefficients(testCase, '');
            got = localSplitsCoefficients(testCase, '3-4i');
            testCase.verifyEqual(got, base, ...
                ['complex edit-field text must fall back to the observed ' ...
                'median n_s and render the identical coefficients table']);
        end

    end
end


function data = localSplitsCoefficients(testCase, nsplitstr)
%Open the data-splits viewer, type nsplitstr into "Splits n'", fire the real
%"Show coefficients" callback, and return the rendered uitable Data.
%
%Values, not handles: the cleanup below closes the figures when this helper
%returns, so a returned handle would already be deleted at the call site.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force'));

pd = PsyRATTestDataFactory.makeSplitsSummaryData();
psyrat_startview_splits('psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd);

viewer = findall(groot, 'Type', 'figure', ...
    'Name', 'Data-Splits Reliability: View Results');
testCase.assertNumElements(viewer, 1, ...
    'expected exactly one data-splits viewer figure.');

localSetEditByLabel(testCase, viewer, 'Splits n', nsplitstr);

btn = findall(viewer, 'Style', 'pushbutton', 'String', 'Show coefficients');
testCase.assertNumElements(btn, 1, 'Show coefficients button not found.');
btn.Callback(btn, []);

outfig = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNumElements(outfig, 1, ...
    'expected exactly one coefficients output figure.');
tbl = findall(outfig, 'Type', 'uitable');
testCase.assertNumElements(tbl, 1, 'expected exactly one uitable in the output.');
data = tbl.Data;
end


function localSetEditByLabel(testCase, viewer, labelprefix, str)
%Type str into the edit uicontrol sharing a layout row with the static text
%starting with labelprefix. The viewer tags no individual controls, so row
%pairing by position is how the existing viewer tests identify them (see
%localSetDynrelTrials in TestExportCallbacks).
lbl = [];
for h = findall(viewer, 'Style', 'text')'
    s = h.String;
    if iscell(s); s = strjoin(s, ' '); end
    if startsWith(string(strtrim(char(s))), labelprefix)
        lbl = h;
        break;
    end
end
testCase.assertNotEmpty(lbl, ...
    sprintf('no "%s" label found on the viewer.', labelprefix));

p = get(lbl, 'Position');
eds = findall(viewer, 'Style', 'edit');
match = arrayfun(@(h) abs(h.Position(2) - p(2)) < 5, eds);
testCase.assertEqual(nnz(match), 1, ...
    sprintf('"%s" edit box not uniquely identified by its row.', labelprefix));
set(eds(match), 'String', str);
end
