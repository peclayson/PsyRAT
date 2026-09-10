classdef TestRelsummaryNoneIdSentinel < PsyRATTestBase
    %Behavioral and source pins for finding B29 -- a PARTICIPANT whose id is
    %literally 'none' must not collide with psyrat_relsummary's no-cutoff-found
    %sentinel.
    %
    %THE DEFECT. The good/bad arm tests were spelled
    %"if ~strcmp(relsummary.group(gloc).goodids,'none')", where the carrier is
    %either the CHAR 'none' (no trial cutoff was found) or a cellstr of retained
    %ids. When the retained ids include a participant named 'none', strcmp
    %returns a logical ARRAY, MATLAB's "if" requires every element true, so
    %~strcmp is false and control falls to the bad-ids arm. On the eventgoodids
    %half that skips the badids append, badids stays empty, and the
    %empty-goodids guard then FIRES -- the user is told the data do not reach
    %the reliability threshold on a run where a cutoff WAS found and
    %participants DID clear it. A silent wrong result presented as a legitimate
    %rejection, not a crash.
    %
    %THE CARRIER HAS FOUR STATES, which is why neither of the file's other two
    %idioms is a correct fix on its own:
    %  1. char 'none'          -- no cutoff found                  -> bad arm
    %  2. empty cell           -- cutoff found, excluded everyone   -> bad arm
    %     (RC-38; strcmp on an empty cell returns [], and "if ~[]" is false)
    %  3. non-empty cellstr    -- normal retained set              -> good arm
    %  4. non-empty cellstr containing 'none'                      -> good arm
    %                                                                 (was bad)
    %"~isempty(x)" alone (:1019, :1080) sends state 1 to the GOOD arm;
    %"~(ischar(x) && strcmp(x,'none'))" alone (:3140) sends state 2 to the GOOD
    %arm. The shipped predicate "iscell(x) && ~isempty(x)" is correct on all
    %four: the good arm is taken only for a genuine non-empty cellstr of ids.
    %B29 in PSYRAT_AUDIT_FINDINGS.md.

    methods (Test)

        function testRenamingAParticipantToNoneDoesNotChangeTheResult(testCase)
            %The behavioral pin, with its positive control built in as the
            %comparison arm. A participant's NAME is not data: renaming one
            %participant to the literal 'none' must leave every reported
            %quantity byte-identical. Asserting the two runs agree is what
            %makes this a differential rather than a single observation --
            %without the control arm a fixture that never reached the site
            %would pass and prove nothing.
            control = testCase.fixtureWithFirstIdRenamed('s99');
            collide = testCase.fixtureWithFirstIdRenamed('none');

            [ctrlOut, ctrlErr] = testCase.runSummary(control);
            [collOut, collErr] = testCase.runSummary(collide);

            %the control must actually reach the good arm, or the comparison
            %below is vacuous (both could be "rejected" and still agree)
            testCase.assertEqual(ctrlErr.nogooddata, 0, ...
                ['The control fixture must find a cutoff and retain ' ...
                'participants, or this test cannot detect B29.']);

            %assert, not verify: on the wrong arm psyrat_relsummary returns
            %early WITHOUT publishing relsummary, so continuing would raise a
            %bare "Unrecognized field name" that buries the real diagnosis.
            testCase.assertEqual(collErr.nogooddata, 0, ...
                ['A participant named ''none'' drove psyrat_relsummary into ' ...
                'the no-cutoff-found arm: the run was reported as failing the ' ...
                'reliability threshold even though a cutoff was found and ' ...
                'participants cleared it (B29).']);

            %every per-cell retained count must match the control
            for gloc = 1:numel(ctrlOut.relsummary.group)
                for eloc = 1:numel(ctrlOut.relsummary.group(gloc).event)
                    testCase.verifyEqual( ...
                        collOut.relsummary.group(gloc).event(eloc).goodn, ...
                        ctrlOut.relsummary.group(gloc).event(eloc).goodn, ...
                        sprintf(['group %d event %d retained a different ' ...
                        'number of participants purely because one was ' ...
                        'named ''none'' (B29).'], gloc, eloc));
                end
            end
        end

        function testTheParticipantNamedNoneIsItselfRetained(testCase)
            %The collision is only interesting because 'none' is a REAL
            %participant. Pin that it survives into the retained set rather
            %than being silently treated as the absent-cutoff sentinel.
            collide = testCase.fixtureWithFirstIdRenamed('none');
            [out, relerr] = testCase.runSummary(collide);

            testCase.assertEqual(relerr.nogooddata, 0, ...
                ['psyrat_relsummary took the no-cutoff-found arm, so no ' ...
                'summary was published to inspect (B29).']);

            retained = out.relsummary.group(1).event(1).eventgoodids;
            testCase.assertTrue(iscell(retained), ...
                ['eventgoodids must be a cellstr of ids here, not the char ' ...
                'sentinel; the fixture did not reach the good arm.']);
            testCase.verifyTrue(any(strcmp(retained, 'none')), ...
                ['The participant literally named ''none'' was dropped from ' ...
                'the retained set (B29).']);
        end

        function testNoSentinelValueKeyedArmTestRemains(testCase)
            %The sweeping source pin: zero value-keyed arm tests anywhere in
            %the file (the pre-fix source carried 22). A future site added with
            %the old idiom reintroduces the defect silently, and no behavioral
            %test covers all 22 routes.
            src = fileread(testCase.relsummarySource());
            pat = ['if\s*~strcmp\(\s*relsummary\.group\(gloc\)' ...
                '(\.event\(eloc\))?\.(event)?goodids\s*,\s*''none''\s*\)'];
            hits = regexp(src, pat, 'match');
            testCase.verifyEmpty(hits, sprintf( ...
                ['%d value-keyed good/bad arm test(s) remain in ' ...
                'psyrat_relsummary.m. The arm must be selected by TYPE and ' ...
                'EMPTINESS -- iscell(x) && ~isempty(x) -- never by comparing ' ...
                'the carrier against the literal ''none'', or a participant ' ...
                'with that id takes the wrong arm (B29).'], numel(hits)));
        end

        function testArmTestsUseTheEmptinessKeyedForm(testCase)
            %Count pin. This is what distinguishes the conversion from mere
            %DELETION of the guards: the same 22 decisions must still be made,
            %in the corrected form.
            src = fileread(testCase.relsummarySource());

            ngood = numel(regexp(src, ...
                ['iscell\(relsummary\.group\(gloc\)\.goodids\)\s*&&\s*' ...
                '~isempty\(relsummary\.group\(gloc\)\.goodids\)'], 'match'));
            nevent = numel(regexp(src, ...
                ['iscell\(relsummary\.group\(gloc\)\.event\(eloc\)\.eventgoodids\)' ...
                '\s*&&\s*~isempty\(relsummary\.group\(gloc\)\.event\(eloc\)' ...
                '\.eventgoodids\)'], 'match'));

            %Counts re-derived 2026-08-29 for the G39 strict fold. The nine
            %cross-event intersection loops no longer test the event carrier
            %at all (a dead event folds in unconditionally), which removed the
            %8 eventgoodids arm tests this census counted (a ninth loop used an
            %ALIASED spelling, "eg = ...eventgoodids; if iscell(eg) && ...",
            %which the regex below never matched -- so the census is a count of
            %the fully-qualified spelling, not of every event-carrier branch,
            %and an aliased reintroduction would pass at 0; pair any suspicion
            %with a grep for "eventgoodids;", which catches the alias
            %assignment); the three case-4 goodn assignments were rekeyed from
            %the event carrier to the group carrier (-3 eventgoodids, +3
            %goodids); and each of the nine loops gained a group-keyed
            %console-notice guard (+9 goodids). 11-8-3 = 0 and 11+3+9 = 23. If
            %either count moves again, re-derive it the same way -- never
            %adjust the pin to make an unreviewed edit pass.
            testCase.verifyEqual(ngood, 23, sprintf( ...
                ['Expected 23 emptiness-keyed goodids arm tests, found %d. ' ...
                'B29 converted 11 goodids sites and G39 added 12 more; a ' ...
                'different count means a site was deleted or missed.'], ngood));
            testCase.verifyEqual(nevent, 0, sprintf( ...
                ['Expected 0 emptiness-keyed eventgoodids arm tests, found ' ...
                '%d. The G39 strict fold removed every event-carrier arm ' ...
                'test; a new one appearing means a consumer is again ' ...
                'branching on a single event''s carrier.'], nevent));
        end

    end

    methods (Access = private)

        function psyrat_data = fixtureWithFirstIdRenamed(~, newid)
            %makeICRelDataMultiGroupEvent is the committed fixture RC-40
            %records as driving these arms live (plain IC, 2 groups x 2
            %events, relsummary 'sing' case 4). Rename its FIRST participant
            %in every per-event table; ids are a cellstr column there, which
            %is what psyrat_normalize_label_columns guarantees on the real
            %loader path too.
            psyrat_data = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
            for e = 1:numel(psyrat_data.rel.data)
                tbl = psyrat_data.rel.data{e};
                tbl.id(strcmp(tbl.id, 's1')) = {newid};
                psyrat_data.rel.data{e} = tbl;
            end
        end

        function [out, relerr] = runSummary(~, psyrat_data)
            %depcutoff 0.01 is the reachable cutoff the committed case-4 test
            %uses, so a cutoff IS found and participants ARE retained -- the
            %only state in which B29's wrong arm is observable.
            [out, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis',    'sing', ...
                'depcutoff',   0.01, ...
                'meascutoff',  2, ...
                'depcentmeas', 1, ...
                'gcoeff',      1, ...
                'CI',          0.95);
        end

        function p = relsummarySource(~)
            %absolute path derived from THIS file's location, never from pwd
            %or the MATLAB path (worktree copies have shadowed the toolbox
            %before)
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here), 'subroutines', 'calculation', ...
                'psyrat_relsummary.m');
        end

    end
end
