classdef TestViewerMissingMeasheader < PsyRATTestBase
    % G112 regression: a saved result that lacks psyrat_data.proc.measheader
    % must open in the setup viewers with the measurement labeled
    % "(not recorded)", not stop with a raw "Unrecognized field name" error.
    %
    % The single-session, test-retest, DoD and criterion setup viewers print
    % ['Measurement:  ' psyrat_data.proc.measheader] unguarded
    % (psyrat_startview_sing.m, _trt.m, _dodiff.m, _criterion.m). Shipped
    % fitting always writes proc.measheader (psyrat_startproc, psyrat_run),
    % but a hand-built or externally produced result can omit it, as the
    % synthetic g1 fixture of the 2026-09-24 native sitting did. The fix
    % backfills the label once in the loader, psyrat_startview, beside the
    % existing 'ver' backfill, so every test here goes through that loader
    % rather than calling a viewer directly.
    %
    % Each variant sets proc itself, whatever the factory fixture carries
    % (the DoD fixture has one; the others do not): removed, a proc struct
    % without measheader, or an empty measheader.
    % A result that does record its measurement is the control: its label
    % must be shown unchanged.

    properties (TestParameter)
        % one test point per viewer family and per missing-label variant, so
        % a failure (or the RED run on the unfixed loader) names each case
        family = {'sing', 'trt', 'dod'};
        variant = {'noproc', 'nomeasheader', 'emptymeasheader'};
    end

    methods (Test)

        function testViewerOpensWithUnrecordedMeasurement(testCase, family, variant)
            % On the unfixed loader every point errors before the viewer is
            % built.
            pd = localResult(testCase, family, variant);
            label = localOpenAndReadMeasurement(testCase, pd, ...
                sprintf('%s/%s', family, variant));
            testCase.verifyEqual(label, 'Measurement:  (not recorded)', ...
                sprintf(['%s viewer, %s: the measurement should read ' ...
                '"(not recorded)", got "%s".'], family, variant, label));
        end

        function testRecordedMeasurementIsShownUnchanged(testCase, family)
            % Control: the backfill must not touch a result that records its
            % measurement.
            pd = localResult(testCase, family, 'recorded');
            label = localOpenAndReadMeasurement(testCase, pd, family);
            testCase.verifyEqual(label, 'Measurement:  meas', sprintf( ...
                '%s viewer: a recorded measurement must be shown as-is, got "%s".', ...
                family, label));
        end

        function testCriterionViewerReceivesTheBackfilledLabel(testCase)
            % The criterion viewer is reached from the single-session viewer's
            % "Criterion Score Outputs" button, whose callback carries the
            % psyrat_data the loader handed to that viewer. Invoking that
            % callback (not a click) checks the pass-through.
            close('all', 'force');
            cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

            pd = localResult(testCase, 'sing', 'noproc');
            psyrat_startview('psyrat_prefs', psyrat_defaults, 'psyrat_data', pd);

            f = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
            testCase.assertNumElements(f, 1, 'the single-session viewer should open');
            btn = localFindButton(f, 'Criterion');
            testCase.assertNotEmpty(btn, ...
                'fixture precondition: the plain one-facet screen offers criterion scoring');

            cb = btn.Callback;
            feval(cb{1}, btn, [], cb{2:end});

            crit = findall(groot, 'Type', 'figure', 'Name', 'Criterion Score Outputs Setup');
            testCase.assertNumElements(crit, 1, 'the criterion viewer should open');
            testCase.verifyEqual(localMeasurementLabel(testCase, crit), ...
                'Measurement:  (not recorded)', ...
                'the criterion viewer should show the backfilled label');
        end

    end
end


function pd = localResult(testCase, family, variant)
%A summarized fixture for one viewer family, with the measurement label set
%up as the variant asks.
switch family
    case 'sing'
        pd = PsyRATTestDataFactory.makeICSummaryData();
    case 'trt'
        pd = PsyRATTestDataFactory.makeTRTSummaryData();
    case 'dod'
        pd = PsyRATTestDataFactory.makeDoDRelData();
end
if ~isfield(pd.rel, 'filename'); pd.rel.filename = 'g112_fixture.psyrat'; end

%fixture preconditions: the right viewer family, and no notice that would
%block on uiwait before the viewer is built (a hang, not a failure, in batch)
testCase.assertEqual(psyrat_analysis_family(pd.rel.analysis), family, ...
    sprintf('fixture precondition: %s routes to the %s viewer', pd.rel.analysis, family));
testCase.assertEmpty(psyrat_engine_notice(pd.rel), ...
    sprintf('fixture precondition: the %s fixture raises no engine notice', family));
[~, requires_ack] = psyrat_convergence_notice(pd.rel);
testCase.assertFalse(logical(requires_ack), ...
    sprintf('fixture precondition: the %s fixture needs no diagnostics acknowledgment', family));

switch variant
    case 'noproc'
        if isfield(pd, 'proc'); pd = rmfield(pd, 'proc'); end
    case 'nomeasheader'
        pd.proc = struct('idheader', 'id');
    case 'emptymeasheader'
        pd.proc = struct('idheader', 'id', 'measheader', '');
    case 'recorded'
        pd.proc = struct('idheader', 'id', 'measheader', 'meas');
end
end


function label = localOpenAndReadMeasurement(testCase, pd, what)
%Open the result through the real loader and return the viewer's
%measurement line. Read inside the helper because the onCleanup below
%closes the figure when this returns.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

psyrat_startview('psyrat_prefs', psyrat_defaults, 'psyrat_data', pd);

f = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
testCase.assertNumElements(f, 1, sprintf('%s: exactly one setup viewer should open', what));
label = localMeasurementLabel(testCase, f);
end


function label = localMeasurementLabel(testCase, fig)
%The one text control whose string starts with "Measurement:".
txt = findall(fig, 'Style', 'text');
label = '';
n = 0;
for k = 1:numel(txt)
    s = txt(k).String;
    if iscell(s); s = strjoin(s, ' '); end
    if startsWith(string(s), 'Measurement:')
        label = char(s);
        n = n + 1;
    end
end
testCase.assertEqual(n, 1, 'the viewer should show exactly one measurement line');
end


function btn = localFindButton(fig, text)
%The push button whose (possibly multi-line) label contains the text.
btn = [];
btns = findall(fig, 'Style', 'pushbutton');
for k = 1:numel(btns)
    s = btns(k).String;
    if iscell(s); s = strjoin(s, ' '); end
    if contains(string(s), text, 'IgnoreCase', true)
        btn = btns(k);
        return;
    end
end
end
