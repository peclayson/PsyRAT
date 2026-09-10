classdef TestRerunDialogDismissal < matlab.unittest.TestCase
    %TESTRERUNDIALOGDISMISSAL Closing the "Chains did not converge" dialog
    %with the window's close box must count as "Do Not Rerun", not become an
    %error after a finished fit.
    %
    % Council finding G13 (2026-08-24), confirmed structurally 2026-09-05:
    % psyrat_reruncheck built its figure with no CloseRequestFcn, so the
    % window's close box ran MATLAB's default closereq, deleted the figure and
    % released uiwait; psyrat_computevarcompwarp then called guidata on the
    % empty handle findobj returned and errored, discarding the estimation
    % that had just completed. The fix registers a close handler that stores
    % the "do not rerun" answer (0) and resumes, leaving the figure for the
    % warp to close and delete as it already does after the two buttons.
    %
    % The close box is exercised the way a user's click reaches it: close(h)
    % runs the figure's CloseRequestFcn. A timer fires it while
    % psyrat_reruncheck is blocked in uiwait (uiwait keeps processing timer
    % callbacks). A second, later timer force-deletes any new figure so a
    % handler that fails to resume cannot hang the test runner.

    methods (Test)
        function testCloseBoxAnswersDoNotRerunAndKeepsTheFigure(testCase)
            before = findall(groot, 'Type', 'figure');
            [t, guard] = localArmTimers(@() localCloseNewDialog(before), before);
            cleanupTimers = onCleanup(@() localStopTimers(t, guard)); %#ok<NASGU>
            psyrat_reruncheck;   % blocks in uiwait until the timer closes it
            figs = setdiff(findall(groot, 'Type', 'figure'), before);
            cleanupFigs = onCleanup(@() localDeleteFigs(figs)); %#ok<NASGU>
            testCase.assertNotEmpty(figs, ...
                ['the dialog must survive its close request so the warp can ', ...
                'read the answer (it was deleted outright before the fix)']);
            h = findobj(figs, 'Tag', 'psyrat_gui');
            testCase.assertNotEmpty(h, 'the surviving dialog must keep its psyrat_gui tag');
            testCase.verifyEqual(guidata(h(1)), 0, ...
                'closing the dialog with the close box must answer "Do Not Rerun" (0)');
        end

        function testRerunButtonStillAnswersOnePositiveControl(testCase)
            before = findall(groot, 'Type', 'figure');
            [t, guard] = localArmTimers(@() localPressButton(before, 'Rerun'), before);
            cleanupTimers = onCleanup(@() localStopTimers(t, guard)); %#ok<NASGU>
            psyrat_reruncheck;
            figs = setdiff(findall(groot, 'Type', 'figure'), before);
            cleanupFigs = onCleanup(@() localDeleteFigs(figs)); %#ok<NASGU>
            h = findobj(figs, 'Tag', 'psyrat_gui');
            testCase.assertNotEmpty(h, 'the dialog must survive the Rerun button');
            testCase.verifyEqual(guidata(h(1)), 1, 'the Rerun button must answer 1');
        end
    end
end

function [t, guard] = localArmTimers(action, before)
%ACTION fires one second in (the dialog is built in milliseconds); GUARD
%force-deletes any new figure ten seconds in so uiwait can never hang the run.
t = timer('StartDelay', 1, 'TimerFcn', @(~,~) action());
guard = timer('StartDelay', 10, 'TimerFcn', @(~,~) localForceDelete(before));
start(guard);
start(t);
end

function localCloseNewDialog(before)
figs = setdiff(findall(groot, 'Type', 'figure'), before);
for k = 1:numel(figs)
    close(figs(k));   % runs the CloseRequestFcn, exactly as the close box does
end
end

function localPressButton(before, label)
figs = setdiff(findall(groot, 'Type', 'figure'), before);
btns = findall(figs, 'Style', 'pushbutton');
for k = 1:numel(btns)
    s = btns(k).String;
    if iscell(s); s = strjoin(s, ' '); end
    if strcmp(strtrim(s), label)
        cb = btns(k).Callback;
        feval(cb{1}, btns(k), [], cb{2:end});
        return;
    end
end
error('test:buttonNotFound', 'button "%s" not found', label);
end

function localForceDelete(before)
figs = setdiff(findall(groot, 'Type', 'figure'), before);
localDeleteFigs(figs);
end

function localDeleteFigs(figs)
if ~isempty(figs)
    delete(figs(isvalid(figs)));
end
end

function localStopTimers(varargin)
for k = 1:numel(varargin)
    t = varargin{k};
    if isvalid(t)
        stop(t);
        delete(t);
    end
end
end
