classdef TestViewerComplexLiteralValidation < PsyRATTestBase
    % Regression guard: complex-literal text typed into a numeric GUI field must
    % be rejected exactly like any other invalid text.
    %
    % WHY THESE GUARDS ALL FAILED THE SAME WAY. str2double parses complex
    % literals -- str2double('3-4i') is 3-4i, and so are str2double('3-4j') and
    % str2double('3 - 4i') -- and every idiom these guards were built from is
    % blind to one (verified on R2025b, 2026-08-11):
    %
    %   isnan(3-4i)              -> false      (only NaN parts count)
    %   isfinite(3-4i)           -> true
    %   3-4i > 0, 3-4i <= 1      -> compares the REAL part only
    %   round(3-4i)              -> 3-4i, so abs(round(x)-x) is 0
    %   max(2,3-4i), min(n,3-4i) -> compare by MAGNITUDE, so clamping and
    %                               rounding do not restore a real value
    %   str2double('3-0i')       -> prints as 3 but isreal is still false
    %
    % so a complex value passed the field's own check and every downstream
    % re-check, then died somewhere far from the edit box -- or, worse, did not
    % die at all. Confirmed live before the fix, one failure mode per site:
    %
    %   dodiff trials    -> psyrat_dodiffrel, "First argument must be an array
    %                       of real values" (MATLAB:prctile:InvalidData)
    %   dodiff gridmax   -> "Colon operands must be real scalars" at 0:gridmax
    %   criterion cutoff -> same prctile error via offsq = (mu - cut).^2, and
    %                       the exported "Criterion Cutoff: %0.4f" header records
    %                       only the real part, so the number reported would not
    %                       have been the number used
    %   trt occasions k  -> mod(val,1) is undefined for complex ("Argument must
    %                       be real"), an unhandled callback error
    %   sing/trt ntrials -> "Colon operands must be real scalars" at
    %                       x = trials(1):trials(2) in the vs-trials plots
    %   reliability cut  -> NO crash: depcheck/relcheck accepted it and it then
    %                       behaved as its real part in find(llrel >= cutoff,1)
    %                       and in the subject-level include/exclude split
    %   trial browser    -> "Complex inputs are not supported" setting the
    %                       slider Value, which aborts the callback before the
    %                       indexing of signals(t,...) that would follow
    %
    % STRATEGY. Every check below is differential rather than message-matching:
    % driving a field with complex text must produce the SAME observable outcome
    % as driving it with plainly invalid text ('abc'). That is the contract each
    % of these guards already documents ("return default on empty/invalid"), and
    % it survives rewording of any dialog.
    %
    % Each differential pair is backed by a POSITIVE CONTROL in the same test: a
    % valid entry must produce the OPPOSITE observable. Without it a differential
    % assertion can pass for the wrong reason -- an earlier draft of this file
    % asserted on "some figure other than this one exists", which is true on the
    % accept path too (the preferences screen closes and the viewer reopens), so
    % those assertions could not fail. The positive control is what makes each
    % observable provably able to distinguish accept from reject.
    %
    % `matlab -batch` builds real figures, so the callbacks can be fired
    % directly (same approach as TestSingViewerCoefficientControl).
    %
    % NOT COVERED HERE. psyrat_startproc's six sampler fields (nchains, nwarmup,
    % nsampling, seed, adapt_delta, max_treedepth) were fixed in the same sweep
    % but cannot be driven from a unit test: psyrat_prefs_save is a subfunction
    % reachable only through the full processing GUI with a loaded dataset. Its
    % adapt_delta hole was the one silent path in the set -- a complex value
    % cleared the GUI check, psyrat_computevarcomp's adapt_delta guard AND
    % MatlabStan's validateattributes, reaching CmdStan as the token
    % "delta=0.5+0.5i". (A complex nchains, by contrast, did NOT reach CmdStan;
    % it died in psyrat_computevarcomp's own mod-based nchains guard.) So of the
    % fourteen conditions the sweep changed, the six here are unpinned.

    methods (Test)

        function testDodTrialCountRejectsComplexLikeGarbage(testCase)
            % local_nonnegint: complex trial counts reached psyrat_dodiffrel.
            localVerifyDodField(testCase, 'Trials for ERP_1:', '50-4i');
        end

        function testDodGridMaxRejectsComplexLikeGarbage(testCase)
            % Same helper, different consumer: gridmax feeds 0:gridmax.
            localVerifyDodField(testCase, 'Heatmap trial grid max:', '20-4i');
        end

        function testCriterionCutoffRejectsComplexLikeGarbage(testCase)
            % psyrat_criterion_cutcheck was isnan-only.
            [complexShown, complexOut] = localRunCriterion(testCase, '2-3i');
            [garbageShown, garbageOut] = localRunCriterion(testCase, 'abc');
            [validShown,   validOut]   = localRunCriterion(testCase, '0.5');

            %positive control first: the observables must be able to say "accepted"
            testCase.assertFalse(validShown, ...
                'positive control: a valid cutoff must NOT raise the error label.');
            testCase.assertGreaterThan(validOut, 0, ...
                'positive control: a valid cutoff must produce criterion output.');

            testCase.verifyEqual(complexShown, garbageShown, ...
                ['a complex criterion cutoff must raise the same inline ' ...
                '"must be numeric" label as any other invalid text.']);
            testCase.verifyTrue(complexShown, ...
                'the criterion cutoff error label should be visible for ''2-3i''.');
            testCase.verifyEqual(complexOut, garbageOut, ...
                'a rejected cutoff must not produce criterion output figures.');
            testCase.verifyEqual(complexOut, 0, ...
                'no criterion output should be generated from ''2-3i''.');
        end

        function testSingReliabilityCutoffRejectsComplexLikeGarbage(testCase)
            % depcheck: isnan + a range test that reads only the real part, so
            % '0.5+0.5i' was accepted and silently behaved as 0.5.
            localVerifyCutoffField(testCase, @psyrat_startview_sing, ...
                PsyRATTestDataFactory.makeICSummaryData(), 'sing/depcheck');
        end

        function testTrtReliabilityCutoffRejectsComplexLikeGarbage(testCase)
            % relcheck: byte-for-byte the same logic as depcheck.
            localVerifyCutoffField(testCase, @psyrat_startview_trt, ...
                PsyRATTestDataFactory.makeTRTSummaryData(), 'trt/relcheck');
        end

        function testTrtOccasionCountRejectsComplexLikeGarbage(testCase)
            % psyrat_get_nocc reached mod(val,1), which errors on complex input.
            %
            % This one is deliberately a THROW detector, and cannot be more than
            % that from this screen: an invalid k falls back to [] (observed
            % occasions), which is also what a blank field gives, and neither the
            % criterion screen nor psyrat_prefs is re-read here, so an accepted-
            % but-wrong k would be indistinguishable from the fallback. What it
            % does pin is the actual pre-fix defect, which was an unhandled
            % MATLAB:math:mustBeReal out of mod() in a pushbutton callback.
            complexOk = localRunTrtNocc(testCase, '2+3i');
            garbageOk = localRunTrtNocc(testCase, 'abc');

            testCase.verifyEqual(complexOk, garbageOk, ...
                ['a complex occasion count k must fall back exactly like any ' ...
                'other invalid text instead of erroring in mod().']);
            testCase.verifyTrue(complexOk, ...
                'the criterion screen should still open with k = ''2+3i''.');
        end

        function testSingPlotTrialsRejectsComplexLikeGarbage(testCase)
            % prefs.view.ntrials flows into x = trials(1):trials(2).
            localVerifyNtrials(testCase, @psyrat_startview_sing, ...
                PsyRATTestDataFactory.makeICSummaryData(), 'sing');
        end

        function testTrtPlotTrialsRejectsComplexLikeGarbage(testCase)
            localVerifyNtrials(testCase, @psyrat_startview_trt, ...
                PsyRATTestDataFactory.makeTRTSummaryData(), 'trt');
        end

        function testTrialBrowserJumpRejectsComplexLikeGarbage(testCase)
            % onJumpToTrial guarded with ~isfinite only, so a complex trial
            % number survived round() and the min/max clamp and was then handed
            % to the slider (which throws on complex) inside localRefresh.
            complexTrial = localRunTrialBrowserJump(testCase, '3-4i');
            garbageTrial = localRunTrialBrowserJump(testCase, 'abc');
            validTrial   = localRunTrialBrowserJump(testCase, '7');

            testCase.assertEqual(validTrial, 7, ...
                'positive control: a valid trial number must actually move currTrial.');

            testCase.verifyEqual(complexTrial, garbageTrial, ...
                ['a complex trial number must leave currTrial untouched, ' ...
                'exactly as unparseable text does.']);
            testCase.verifyTrue(isreal(complexTrial), ...
                'currTrial must never become complex.');
            testCase.verifyEqual(complexTrial, 1, ...
                'currTrial should still be the initial trial after a rejected jump.');
        end

    end
end


% ---------------------------------------------------------------------------
% Per-viewer drivers. Each returns plain values, never handles: the onCleanup
% inside closes the figures before the caller sees the result.
% ---------------------------------------------------------------------------

function localVerifyDodField(testCase, labelprefix, complextext)
%Type complex text into one DoD numeric field, click Generate, and require the
%same outcome as typing 'abc' there: an error dialog and no output window.
[complexDlg, complexOut] = localRunDod(testCase, labelprefix, complextext);
[garbageDlg, garbageOut] = localRunDod(testCase, labelprefix, 'abc');
[validDlg,   validOut]   = localRunDod(testCase, labelprefix, '');

%positive control: with the field left at its shipped default the same click
%must produce output and no dialog, so both observables can say "accepted".
testCase.assertFalse(validDlg, ...
    'positive control: valid DoD inputs must not raise a validation dialog.');
testCase.assertGreaterThan(validOut, 0, ...
    'positive control: valid DoD inputs must produce an output window.');

testCase.verifyEqual(complexDlg, garbageDlg, sprintf( ...
    ['"%s" = ''%s'' must raise the same validation dialog as invalid text; ' ...
    'before the isreal fix it raised an uncaught error instead.'], ...
    labelprefix, complextext));
testCase.verifyTrue(complexDlg, sprintf( ...
    '"%s" = ''%s'' should be rejected with a dialog.', labelprefix, complextext));
testCase.verifyEqual(complexOut, garbageOut, ...
    'a rejected DoD input must not produce output figures.');
testCase.verifyEqual(complexOut, 0, ...
    'no DoD output should be generated from a rejected trial input.');
end


function [sawDialog, nOutputs] = localRunDod(testCase, labelprefix, text)
%Open the DoD setup screen, set one field, and fire "Generate DoD Outputs".
%An empty `text` leaves the field at its shipped default (the positive control).
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

psyrat_startview_dodiff('psyrat_prefs', psyrat_defaults, ...
    'psyrat_data', localStamp(PsyRATTestDataFactory.makeDoDRelData()));

fig = localViewerFigure(testCase, 'DoD Reliability Setup');
if ~isempty(text)
    set(localEditByLabel(testCase, fig, labelprefix), 'String', text);
end

btn = localButton(testCase, fig, 'Generate DoD Outputs');
localFire(btn);

%errordlg builds a normal figure; anything that is neither the setup screen nor
%a DoD output window is the validation dialog. Unlike the preferences screens,
%this callback leaves the setup figure open on BOTH paths, so this really does
%isolate the dialog (the positive control above proves it).
sawDialog = ~isempty(localOtherFigures(fig, 'DoD Reliability Outputs'));
nOutputs = numel(findall(groot, 'Type', 'figure', 'Name', 'DoD Reliability Outputs'));
end


function [labelVisible, nOutputs] = localRunCriterion(testCase, text)
%Open the criterion screen, set the cutoff, and fire "Generate Criterion Outputs".
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

psyrat_startview_criterion('psyrat_prefs', psyrat_defaults, ...
    'psyrat_data', localStamp(PsyRATTestDataFactory.makeICSummaryData()));

fig = localViewerFigure(testCase, 'Criterion Score Outputs Setup');
set(localEditByLabel(testCase, fig, 'Criterion Cutoff:'), 'String', text);

btn = localButton(testCase, fig, 'Generate Criterion Outputs');
localFire(btn);

%the inline warning label is created hidden next to the cutoff field
lbl = localTextByPrefix(fig, 'Criterion cutoff must be numeric');
testCase.assertNotEmpty(lbl, 'the criterion cutoff warning label was not found.');
labelVisible = strcmp(lbl.Visible, 'on');
nOutputs = numel(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
end


function localVerifyCutoffField(testCase, viewerfcn, pd, tag)
%The reliability-cutoff box on the main viewer, checked by depcheck (sing) or
%relcheck (trt). A complex value here never crashed -- it was silently accepted
%and used as its real part -- so the observable has to distinguish accepted from
%rejected, not merely "did it throw".
[complexShown, complexOpen] = localRunCutoff(testCase, viewerfcn, pd, '0.5+0.5i');
[garbageShown, garbageOpen] = localRunCutoff(testCase, viewerfcn, pd, 'abc');
[validShown,   validOpen]   = localRunCutoff(testCase, viewerfcn, pd, '0.5');

testCase.assertFalse(validShown, sprintf( ...
    '%s positive control: a valid cutoff must not raise the error label.', tag));
testCase.assertTrue(validOpen, sprintf( ...
    '%s positive control: a valid cutoff must open the criterion screen.', tag));

testCase.verifyEqual(complexShown, garbageShown, sprintf( ...
    ['%s: a complex reliability cutoff must raise the same inline error as ' ...
    'any other invalid text. It was silently ACCEPTED before the fix and ' ...
    'then behaved as its real part in find(llrel >= cutoff,1).'], tag));
testCase.verifyTrue(complexShown, sprintf( ...
    '%s: ''0.5+0.5i'' should raise the cutoff error label.', tag));
testCase.verifyEqual(complexOpen, garbageOpen, sprintf( ...
    '%s: a rejected cutoff must not proceed to the criterion screen.', tag));
testCase.verifyFalse(complexOpen, sprintf( ...
    '%s: ''0.5+0.5i'' must not be treated as a usable cutoff.', tag));
end


function [labelVisible, criterionOpened] = localRunCutoff(testCase, viewerfcn, pd, text)
%Set the reliability cutoff on the main viewer and fire "Criterion Score
%Outputs", which is the callback that runs depcheck/relcheck.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

viewerfcn('psyrat_prefs', psyrat_defaults, 'psyrat_data', localStamp(pd));

fig = localViewerFigure(testCase, 'View Results Setup');
set(localEditByLabel(testCase, fig, 'Reliability Cutoff:'), 'String', text);

%grab the hidden warning label BEFORE firing: on the accept path the criterion
%screen closes the whole viewer (psyrat_startview_criterion closes any figure
%tagged psyrat_gui), so the label cannot be looked up afterwards.
lbl = localTextByPrefix(fig, 'Reliability cutoff must be numeric');
testCase.assertNotEmpty(lbl, 'the reliability cutoff warning label was not found.');

btn = localButton(testCase, fig, 'Criterion Score Outputs');
localFire(btn);

%a deleted label means the viewer was torn down, i.e. the cutoff was accepted
labelVisible = isgraphics(lbl) && strcmp(lbl.Visible, 'on');
criterionOpened = ~isempty(findall(groot, 'Type', 'figure', ...
    'Name', 'Criterion Score Outputs Setup'));
end


function opened = localRunTrtNocc(testCase, text)
%Set the occasion-composite k field and fire "Criterion Score Outputs", which
%routes through psyrat_get_nocc. True when the next screen opened without error.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

psyrat_startview_trt('psyrat_prefs', psyrat_defaults, ...
    'psyrat_data', localStamp(PsyRATTestDataFactory.makeTRTSummaryData()));

fig = localViewerFigure(testCase, 'View Results Setup');

%select the multi-occasion composite so the k field is the live one
pops = findall(fig, 'Style', 'popupmenu');
for k = 1:numel(pops)
    s = pops(k).String;
    if iscell(s) && numel(s) > 1 && contains(lower(strjoin(s, ' ')), 'composite')
        pops(k).Value = 2;
        if ~isempty(pops(k).Callback)
            localFire(pops(k));
        end
    end
end

set(localEditByLabel(testCase, fig, '# Occasions in Composite'), 'String', text);

btn = localButton(testCase, fig, 'Criterion Score Outputs');
localFire(btn);

opened = ~isempty(findall(groot, 'Type', 'figure', ...
    'Name', 'Criterion Score Outputs Setup'));
end


function localVerifyNtrials(testCase, viewerfcn, pd, tag)
%The plot-trials field lives on the "Output Prefs" sub-screen; saving it is what
%validates it. A refused save keeps that screen open, an accepted one closes it
%and returns to the viewer -- that is the observable, and the positive control
%below is what proves it discriminates.
complexOpen = localRunNtrials(testCase, viewerfcn, pd, '50-4i');
garbageOpen = localRunNtrials(testCase, viewerfcn, pd, 'abc');
validOpen   = localRunNtrials(testCase, viewerfcn, pd, '25');

testCase.assertFalse(validOpen, sprintf( ...
    ['%s positive control: a valid trials count must be accepted and close ' ...
    'the preferences screen.'], tag));

testCase.verifyEqual(complexOpen, garbageOpen, sprintf( ...
    ['%s: a complex plot-trials count must be refused exactly like invalid ' ...
    'text rather than being stored and reaching the vs-trials plot.'], tag));
testCase.verifyTrue(complexOpen, sprintf( ...
    ['%s: ''50-4i'' should be refused, leaving the preferences screen open.'], tag));
end


function prefsStillOpen = localRunNtrials(testCase, viewerfcn, pd, text)
%Open a viewer, go to "Output Prefs", set the trials field, and click Save.
%Returns whether the preferences screen survived the click.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

viewerfcn('psyrat_prefs', psyrat_defaults, 'psyrat_data', localStamp(pd));

main = localViewerFigure(testCase, 'View Results Setup');
btn = localButton(testCase, main, 'Output Prefs');
localFire(btn);

prefsfig = localPrefsFigure(testCase);
edit = localEditByLabel(testCase, prefsfig, 'Number of');
set(edit, 'String', text);

savebtn = localSaveButton(testCase, prefsfig);
localFire(savebtn);

prefsStillOpen = isgraphics(prefsfig);
end


function currTrial = localRunTrialBrowserJump(testCase, text)
%Fire the trial browser's jump-to-trial callback and read currTrial back.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force')); %#ok<NASGU>

psyrat_import_trialbrowser('project', localSyntheticProject());

figs = findall(groot, 'Type', 'figure');
testCase.assertNumElements(figs, 1, 'the trial browser should open one figure.');
fig = figs(1);

%by label, not by position in the findall list, so this cannot silently
%retarget if the browser gains another edit control
hedit = localEditByLabel(testCase, fig, 'Trial #');

set(hedit, 'String', text);
localFire(hedit);

st = guidata(fig);
currTrial = st.currTrial;
end


% ---------------------------------------------------------------------------
% Shared helpers
% ---------------------------------------------------------------------------

function pd = localStamp(pd)
%The viewers print a dataset header the bare fixtures do not carry (the same
%stamping TestSingViewerCoefficientControl does).
if ~isfield(pd, 'ver'); pd.ver = '0-test'; end
if ~isfield(pd.rel, 'filename'); pd.rel.filename = 'test'; end
if ~isfield(pd, 'proc') || ~isfield(pd.proc, 'measheader')
    pd.proc.measheader = 'meas';
end
end


function localFire(h)
%Invoke a uicontrol's Callback. These viewers wire buttons as cell arrays
%({@fcn, 'psyrat_prefs', prefs, ...}) rather than bare function handles, so the
%trailing entries have to be unpacked as extra arguments; the trial browser uses
%a bare handle. Handle both.
cb = h.Callback;
if iscell(cb)
    fcn = cb{1};
    fcn(h, [], cb{2:end});
else
    cb(h, []);
end
end


function fig = localViewerFigure(testCase, name)
fig = findall(groot, 'Type', 'figure', 'Name', name);
testCase.assertNumElements(fig, 1, ...
    sprintf('expected exactly one "%s" figure.', name));
end


function fig = localPrefsFigure(testCase)
%The output-preferences screen is the figure that is not the main viewer and
%carries an edit field; it is named differently across viewers.
figs = findall(groot, 'Type', 'figure');
fig = [];
for k = 1:numel(figs)
    if ~strcmp(figs(k).Name, 'View Results Setup') && ...
            ~isempty(findall(figs(k), 'Style', 'edit'))
        fig = figs(k);
        break;
    end
end
testCase.assertNotEmpty(fig, 'the output-preferences screen did not open.');
end


function btn = localButton(testCase, fig, label)
btns = findall(fig, 'Style', 'pushbutton');
btn = [];
for k = 1:numel(btns)
    s = btns(k).String;
    if iscell(s); s = strjoin(s, ' '); end
    if strcmp(strtrim(regexprep(char(s), '\s+', ' ')), label)
        btn = btns(k);
        return;
    end
end
testCase.assertNotEmpty(btn, sprintf('button "%s" not found.', label));
end


function btn = localSaveButton(testCase, fig)
%The preferences screens label their commit button with "Save"; take whichever
%pushbutton mentions it.
btns = findall(fig, 'Style', 'pushbutton');
btn = [];
for k = 1:numel(btns)
    s = btns(k).String;
    if iscell(s); s = strjoin(s, ' '); end
    if contains(lower(char(s)), 'save')
        btn = btns(k);
        return;
    end
end
testCase.assertNotEmpty(btn, 'no Save button on the preferences screen.');
end


function hedit = localEditByLabel(testCase, fig, labelprefix)
%Return the edit uicontrol on the same row as, and to the right of, the static
%text starting with labelprefix. The viewers tag no individual controls, so row
%pairing by position is the established way to identify them.
%
%Pairing is by vertical OVERLAP of the two control rectangles rather than by
%their bottom edges: the labels are taller than the edit boxes and the two are
%not bottom-aligned (on the output-preferences screens the bottoms differ by 16
%pixels), so a bottom-edge tolerance silently finds nothing. Among the
%overlapping candidates to the label's right, take the one sharing the most
%vertical extent, then the nearest.
txts = findall(fig, 'Style', 'text');
edits = findall(fig, 'Style', 'edit');
hedit = [];
bestoverlap = 0;
bestdx = inf;
for j = 1:numel(txts)
    s = txts(j).String;
    if iscell(s); s = strjoin(s, ' '); end
    s = strtrim(char(s));
    if ~strncmp(s, labelprefix, numel(labelprefix)); continue; end
    tp = txts(j).Position;
    for k = 1:numel(edits)
        ep = edits(k).Position;
        dx = ep(1) - tp(1);
        overlap = min(tp(2) + tp(4), ep(2) + ep(4)) - max(tp(2), ep(2));
        if dx <= 0 || overlap <= 0; continue; end
        if overlap > bestoverlap || (overlap == bestoverlap && dx < bestdx)
            bestoverlap = overlap;
            bestdx = dx;
            hedit = edits(k);
        end
    end
end
testCase.assertNotEmpty(hedit, ...
    sprintf('no edit field found on the "%s" row.', labelprefix));
end


function lbl = localTextByPrefix(fig, prefix)
lbl = [];
txts = findall(fig, 'Style', 'text');
for j = 1:numel(txts)
    s = txts(j).String;
    if iscell(s); s = strjoin(s, ' '); end
    if strncmp(strtrim(char(s)), prefix, numel(prefix))
        lbl = txts(j);
        return;
    end
end
end


function figs = localOtherFigures(known, alsoKnownName)
%Every open figure that is neither `known` nor named alsoKnownName -- i.e. the
%dialog a validation failure puts up. Only sound where the callback leaves
%`known` open on both the accept and reject paths (see localRunDod).
all_figs = findall(groot, 'Type', 'figure');
keep = false(size(all_figs));
for k = 1:numel(all_figs)
    if isequal(all_figs(k), known); continue; end
    if ~isempty(alsoKnownName) && strcmp(all_figs(k).Name, alsoKnownName)
        continue;
    end
    keep(k) = true;
end
figs = all_figs(keep);
end


function project = localSyntheticProject()
%Minimal canonical project the trial browser accepts: it only needs single-trial
%signals, a time axis, channel labels and a matching trial table.
nTrials = 10;
nChans = 2;
nSamp = 25;
signals = zeros(nTrials, nChans, nSamp);
for t = 1:nTrials
    for c = 1:nChans
        signals(t, c, :) = sin(linspace(0, pi, nSamp)) * (t + c);
    end
end

project = struct();
project.canonical = struct();
project.canonical.has_single_trial = true;
project.canonical.single_trial = struct( ...
    'signals', signals, ...
    'time_ms', linspace(-100, 400, nSamp)', ...
    'channel_labels', {{'Cz', 'Fz'}});
%localRefresh prints a metadata line, so the table needs every column it reads:
%subject_id, event_id, occasion_id, group_id and trial_uid.
project.canonical.trial_table = table( ...
    repmat({'S01'}, nTrials, 1), ...
    repmat({'ERN'}, nTrials, 1), ...
    repmat({'1'}, nTrials, 1), ...
    repmat({'control'}, nTrials, 1), ...
    arrayfun(@(k) sprintf('uid%02d', k), (1:nTrials)', 'UniformOutput', false), ...
    'VariableNames', ...
    {'subject_id', 'event_id', 'occasion_id', 'group_id', 'trial_uid'});
project.scoring_specs = [];
end
