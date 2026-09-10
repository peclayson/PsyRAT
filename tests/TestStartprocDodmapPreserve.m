classdef TestStartprocDodmapPreserve < PsyRATTestBase
    %Source pins for finding B20 -- the Select Events Save callback must not
    %discard a configured DoD contrast when the saved event SET is unchanged,
    %and a real wipe must be visible.
    %
    %WHY SOURCE SCANS. psyrat_whichevents_save_call is a LOCAL function of
    %psyrat_startproc.m, so it cannot be called from a test (the standing
    %local-functions limitation), and exercising it for real means clicking
    %through the live GUI -- which is exactly the deferred verification this
    %class does not claim. These pins hold the fix's structure in place until
    %the click-through session: the wipe is guarded by a set comparison, the
    %comparison is the shared normalizer applied to BOTH selections, and the
    %real-wipe arm warns modally. Pins were watched to FAIL against the
    %pre-fix source (verified by substituting the HEAD copy, 2026-08-16).
    %
    %WHAT THE FIX IS FOR. Every Save used to run
    %    psyrat_prefs.proc.dodset = 0;
    %    psyrat_prefs.proc.dodmap = {};
    %unconditionally; the only signal was a tooltip flip, and re-opening the
    %configurator seeded a DIFFERENT valid mapping (data first-appearance
    %order, made alphabetical by the Save itself). Observed live 2026-08-13;
    %B20 in PSYRAT_AUDIT_FINDINGS.md.

    methods (Test)

        function testWipeIsGuardedBySetComparison(testCase)
            body = testCase.saveCallbackBody();

            %the guard exists and the wipe lines sit AFTER it
            guard = strfind(body, 'if selectionchanged');
            testCase.assertNotEmpty(guard, ...
                ['the Save callback must guard the dodset/dodmap wipe with ' ...
                'a selection-change test; an unconditional wipe restores B20']);
            wipe1 = strfind(body, 'psyrat_prefs.proc.dodset = 0;');
            wipe2 = strfind(body, 'psyrat_prefs.proc.dodmap = {};');
            testCase.assertTrue(~isempty(wipe1) && ~isempty(wipe2), ...
                'the wipe lines are gone entirely; the changed-set arm must keep them');
            testCase.verifyTrue(wipe1(1) > guard(1) && wipe2(1) > guard(1), ...
                'the wipe must be INSIDE the selection-changed arm');

            %and the pre-fix unconditional sequence is absent
            prefix = ['psyrat_prefs.proc.inp.whichevents = events;' newline ...
                'psyrat_prefs.proc.dodset = 0;'];
            testCase.verifyEmpty(strfind(body, prefix), ...
                'the unconditional assignment-then-wipe sequence has returned');
        end

        function testComparisonIsSetBasedOnBothSelections(testCase)
            body = testCase.saveCallbackBody();

            %both the prior and the new selection go through the shared
            %normalizer, so the comparison is order- and type-insensitive
            hits = strfind(body, 'psyrat_eventset_key(');
            testCase.verifyGreaterThanOrEqual(numel(hits), 2, ...
                ['both selections must be normalized through ' ...
                'psyrat_eventset_key so the comparison is a SET comparison']);

            %the helper itself exists in the file
            src = fileread(testCase.startprocSource());
            testCase.verifyNotEmpty( ...
                strfind(src, 'function keys = psyrat_eventset_key(v)'), ...
                'the psyrat_eventset_key normalizer is missing');
        end

        function testPreserveRequiresMapConfiguredForCurrentColumn(testCase)
            %Refuter-caught hole in the first cut of the fix: switching the
            %EVENT COLUMN rewrites whichevents outside the Save callback
            %without wiping the map (psyrat_store_inputchoices writes no dod
            %fields), so set-equality alone could preserve a dodset=1 mapping
            %whose names come from ANOTHER column -- inert downstream thanks
            %to the dodmapcol triple check, but silently retained and then
            %silently lost. The preserve arm must therefore also require the
            %stored map to be configured for the CURRENT event column.
            body = testCase.saveCallbackBody();
            testCase.verifyNotEmpty(strfind(body, 'mapcurrent'), ...
                ['the preserve arm must verify the map was configured for ' ...
                'the current event column, not rely on set-equality alone']);
            testCase.verifyNotEmpty(strfind(body, 'dodmapcol'), ...
                'the current-column check must key on proc.dodmapcol');
        end

        function testRealWipeWarnsModally(testCase)
            body = testCase.saveCallbackBody();

            warnpos = strfind(body, 'warndlg');
            testCase.assertNotEmpty(warnpos, ...
                ['a real wipe of a configured contrast must be VISIBLE; ' ...
                'the tooltip flip was the whole of B20''s silent half']);
            %modal, per the B21 class ruling, and gated on a configured map
            %(wiping nothing needs no dialog)
            testCase.verifyNotEmpty(strfind(body, '''DoD contrast cleared'',''modal'''), ...
                'the wipe warning must be modal and carry its title');
            hadmap = strfind(body, 'hadmap');
            testCase.verifyTrue(~isempty(hadmap) && hadmap(1) < warnpos(1), ...
                'the warning must fire only when a contrast was configured');
        end

        function testWipeWarningFiresAfterTheRelaunch(testCase)
            %Live click-through catch (2026-08-16): MATLAB modality binds only
            %the windows that exist when the dialog is created, and a figure
            %created afterwards stacks ABOVE it. The first cut fired the
            %warndlg BEFORE the callback relaunched psyrat_startproc_gui, so
            %the relaunched setup screen both covered the warning and escaped
            %its modality -- observed on screen: the dialog was buried and the
            %Event popup still opened. The dialog must therefore be created
            %LAST, after the relaunch.
            body = testCase.saveCallbackBody();
            warnpos = strfind(body, 'warndlg');
            relaunch = strfind(body, 'psyrat_startproc_gui(');
            testCase.assertTrue(~isempty(warnpos) && ~isempty(relaunch), ...
                'the warning or the relaunch is missing from the Save callback');
            testCase.verifyTrue(warnpos(1) > relaunch(1), ...
                ['the wipe warning must be created AFTER the relaunched setup ' ...
                'screen, or the new figure buries it and escapes its modality']);
        end

    end

    methods (Access = private)

        function p = startprocSource(~)
            %absolute path derived from THIS file's location, never from pwd
            %or the MATLAB path (worktree copies have shadowed the toolbox
            %before)
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here), 'subroutines', 'guis', ...
                'psyrat_startproc.m');
        end

        function body = saveCallbackBody(testCase)
            %extract psyrat_whichevents_save_call's body: from its function
            %line to the next top-level function declaration
            src = fileread(testCase.startprocSource());
            a = strfind(src, 'function psyrat_whichevents_save_call');
            testCase.assertNotEmpty(a, ...
                'psyrat_whichevents_save_call not found; the anchor has drifted');
            rest = src(a(1):end);
            nextfun = regexp(rest(2:end), '\nfunction ', 'once');
            testCase.assertNotEmpty(nextfun, ...
                'could not delimit the callback body');
            body = rest(1:nextfun);
        end

    end
end
