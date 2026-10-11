classdef TestStartprocPreflightWarningDialog < matlab.unittest.TestCase
    %TESTSTARTPROCPREFLIGHTWARNINGDIALOG The GUI's estimation callback must
    %raise the preflight's non-fatal warnings as a modal dialog before the
    %fit starts.
    %
    % Council finding G8 (2026-08-24), confirmed 2026-09-05:
    % psyrat_preflight_validate computed few-participants, low-median-trials
    % and incomplete-cell warnings that both callers dropped, reading only
    % errors(1). The scripted route is pinned behaviorally in
    % TestInstallationAndWorkflowValidation (the warp prints them). The GUI
    % dialog cannot be driven headlessly (the callback needs a built Specify
    % Inputs screen and a save panel), so this class pins the source the way
    % TestModalRejectionDialogs pins the rejection dialogs: the block between
    % the preflight call and the warp call must read preflight.warnings and
    % raise uiwait(warndlg(..., 'modal')), so a multi-hour run does not start
    % until the user has seen the warning.

    methods (Test)
        function testEstimationCallbackRaisesPreflightWarningsModally(testCase)
            src = localReadSource('subroutines/guis/psyrat_startproc.m');
            span = localPreflightToWarpSpan(testCase, src);
            %control: the span is the right block (it already reads errors)
            testCase.assertTrue(contains(span, 'preflight.errors(1)'), ...
                'span control: the preflight-to-warp block must be the one that reads errors(1)');
            testCase.verifyTrue(contains(span, 'preflight.warnings'), ...
                'the estimation callback must read the preflight warnings');
            testCase.verifyTrue(~isempty(regexp(span, 'uiwait\s*\(\s*warndlg\s*\(', 'once')), ...
                ['the warnings must be raised through uiwait(warndlg(...)) so the ', ...
                'run waits for acknowledgment']);
            testCase.verifyTrue(~isempty(regexp(span, 'warndlg\s*\([^;]*''modal''', 'once')), ...
                'the warning dialog must be modal');
        end

        function testEstimationCallbackRaisesRestrictedNoticeBeforeTheSaveDialog(testCase)
            %Queue row 3a (owner rulings 2-4, 2026-09-21): the RESTRICTED
            %design-maturity notice for the gamma two-facet designs is
            %raised ONCE, through uiwait(warndlg(..., 'modal')) like the
            %preflight warnings, and BEFORE the uiputfile save dialog, so
            %the user can still back out with Cancel. The behavior is
            %driven through the real callback in TestStartprocAnalyzeDrive;
            %this pins the placement, which a drive cannot see.
            src = localReadSource('subroutines/guis/psyrat_startproc.m');
            i0 = strfind(src, 'psyrat_gamma_twofacet_restriction(');
            testCase.assertEqual(numel(i0), 1, ...
                'expected exactly one restricted-notice site in psyrat_startproc');
            isave = strfind(src, 'uiputfile(');
            testCase.assertNotEmpty(isave, 'no save dialog in psyrat_startproc');
            testCase.verifyLessThan(i0, min(isave), ...
                'the restricted-design notice must precede the save dialog');
            span = src(i0:min(isave));
            testCase.verifyTrue(~isempty(regexp(span, 'uiwait\s*\(\s*warndlg\s*\(', 'once')), ...
                'the notice must be raised through uiwait(warndlg(...)) so the run waits for acknowledgment');
            testCase.verifyTrue(~isempty(regexp(span, 'warndlg\s*\([^;]*''modal''', 'once')), ...
                'the notice dialog must be modal');
        end
    end
end

function src = localReadSource(relpath)
%absolute path from THIS file's location, never from pwd or the MATLAB path
%(worktree copies have shadowed the toolbox before)
here = fileparts(mfilename('fullpath'));
src = fileread(fullfile(fileparts(here), relpath));
end

function span = localPreflightToWarpSpan(testCase, src)
%The GUI runs preflight once per measurement and then calls the warp; the
%warnings must be raised in between.
i0 = strfind(src, 'preflight = psyrat_preflight_validate(');
testCase.assertEqual(numel(i0), 1, ...
    'expected exactly one preflight call in psyrat_startproc');
rest = src(i0:end);
i1 = strfind(rest, 'psyrat_computevarcompwarp(');
testCase.assertNotEmpty(i1, 'no warp call after the preflight call');
span = rest(1:i1(1));
end
