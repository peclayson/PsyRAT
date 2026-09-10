classdef TestStartprocNoneSentinel < PsyRATTestBase
    %Source pins for finding B25 -- every facet-absent ("none") guard in
    %psyrat_startproc.m must key on the sentinel INDEX (the last entry of
    %psyrat_data.proc.collist), never on the selected column's NAME.
    %
    %WHY SOURCE SCANS. All sixteen guard sites live in LOCAL functions of
    %psyrat_startproc.m (GUI callbacks and status helpers), so they cannot be
    %called from a test (the standing local-functions limitation), and
    %exercising them for real means clicking through the live GUI. These pins
    %hold the convention in place until that click-through. Pins were watched
    %to FAIL against the pre-fix source (2026-08-17).
    %
    %WHAT THE FIX IS FOR. The role popups append a 'none' sentinel to the
    %column list, so a data column literally named 'none' appears twice and is
    %selectable at its own index. The header resolution and the noneind idiom
    %already keyed on the INDEX (selection == length(collist)); the sixteen
    %design guards keyed on the NAME (strcmpi(collist{ind},'none')). With a
    %real 'none' column selected, the two conventions split: the guards
    %conclude the facet is absent while the loader loads the column, so the
    %engine self-dispatches a WITH-facet analysis every GUI gate believes it
    %prevented -- a silent estimand substitution. B25 in
    %PSYRAT_AUDIT_FINDINGS.md.

    methods (Test)

        function testNoStringKeyedNoneGuardRemains(testCase)
            %the sweeping pin: zero name-keyed sentinel compares anywhere in
            %the file (pre-fix source carried sixteen)
            src = fileread(testCase.startprocSource());
            pat = 'strcmpi?\s*\(\s*psyrat_data\.proc\.collist\s*\{[^}]+\}\s*,\s*''none''\s*\)';
            hits = regexp(src, pat, 'match');
            testCase.verifyEmpty(hits, sprintf( ...
                ['%d name-keyed ''none'' guard(s) remain in ' ...
                'psyrat_startproc.m; the sentinel must be detected by INDEX ' ...
                '(selection == length(collist)) or a real column named ' ...
                '''none'' re-splits the guards from the loader (B25)'], ...
                numel(hits)));
        end

        function testAutodetectCandidateListsNeverNameTheSentinel(testCase)
            %The precondition the ncol+1 anchor comment depends on. The
            %participant and measurement auto-detect loops iterate
            %1:length(proc.collist) -- the sentinel-extended list, ncol+1
            %entries -- but inp.id and inp.meas are consumed against
            %raw.colnames, which has only ncol (the Participant ID popup is
            %built with 'String',psyrat_data.raw.colnames). The extra
            %iteration is inert only while no candidate list can match the
            %sentinel. If one ever named it, the role index would be set to
            %ncol+1 and the popup Value would exceed its String.
            %
            %Scanned rather than called: the auto-detect block is inline in
            %psyrat_startproc, not a callable local (the standing
            %local-functions limitation this class exists to work around).
            src = fileread(testCase.startprocSource());

            lists = regexp(src, 'poss\s*=\s*\{[^}]*\}', 'match');
            testCase.verifyNotEmpty(lists, ...
                'The candidate lists must be findable, or this pin is vacuous.');

            offenders = lists(~cellfun(@isempty, ...
                regexpi(lists, '''\s*\(?\s*none\s*\)?\s*''', 'once')));
            testCase.verifyEmpty(offenders, sprintf( ...
                ['%d auto-detect candidate list(s) name the ''none'' ' ...
                'sentinel. A match there sets inp.id/inp.meas to ncol+1, ' ...
                'which is out of range for raw.colnames and for the ' ...
                'Participant ID popup''s String.'], numel(offenders)));
        end

        function testSubGuiGuardsUseTheIndexConvention(testCase)
            %the four sub-GUI entry guards (Select Groups / Select Events /
            %DoD contrast / Select Occasions) each carry exactly one
            %index-form facet-absent check
            expected = { ...
                'selectgroups_call', ...
                'psyrat_prefs.proc.inp.group == length(psyrat_data.proc.collist)'; ...
                'selectevents_call', ...
                'psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)'; ...
                'dodcontrast_call', ...
                'psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)'; ...
                'selecttime_call', ...
                'psyrat_prefs.proc.inp.time == length(psyrat_data.proc.collist)'};
            for k = 1:size(expected, 1)
                body = testCase.functionBody(expected{k, 1});
                testCase.verifyEqual( ...
                    numel(strfind(body, expected{k, 2})), 1, sprintf( ...
                    ['%s must guard facet absence with exactly one index ' ...
                    'comparison (%s); the name-keyed guard restores B25'], ...
                    expected{k, 1}, expected{k, 2}));
            end
        end

        function testExecFacetChecksUseTheIndexConvention(testCase)
            %psyrat_exec (the Analyze path) carries the header resolution
            %(one ==/~= pair per role, pre-existing) PLUS the nine converted
            %guards: the three which-X population conditions, the three
            %which-X handoffs, and the three design gates. Exact counts pin
            %both the conversion and its completeness.
            body = testCase.functionBody('psyrat_exec');
            expected = { ...
                'psyrat_prefs.proc.inp.group == length(psyrat_data.proc.collist)', 1; ...
                'psyrat_prefs.proc.inp.group ~= length(psyrat_data.proc.collist)', 3; ...
                'psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)', 2; ...
                'psyrat_prefs.proc.inp.event ~= length(psyrat_data.proc.collist)', 3; ...
                'psyrat_prefs.proc.inp.time == length(psyrat_data.proc.collist)', 1; ...
                'psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist)', 5};
            for k = 1:size(expected, 1)
                testCase.verifyEqual( ...
                    numel(strfind(body, expected{k, 1})), expected{k, 2}, ...
                    sprintf(['psyrat_exec must contain exactly %d of "%s" ' ...
                    '(header resolution plus the B25-converted guards); a ' ...
                    'different count means a guard was missed, dropped, or ' ...
                    'reworded away from the index convention'], ...
                    expected{k, 2}, expected{k, 1}));
            end
        end

        function testStatusHelpersUseTheIndexConvention(testCase)
            %psyrat_subsetstatus folds its stale-index defence and the old
            %name-keyed check into one index guard (>= covers both "the
            %sentinel" and "beyond the list")
            body = testCase.functionBody('psyrat_subsetstatus');
            testCase.verifyEqual(numel(strfind(body, ...
                'colind >= length(psyrat_data.proc.collist)')), 1, ...
                ['psyrat_subsetstatus must treat the facet as unused only ' ...
                'when colind >= length(collist): the sentinel index or a ' ...
                'stale out-of-range preference, never a column named ''none''']);
            testCase.verifyEmpty(strfind(body, ...
                'colind > length(psyrat_data.proc.collist)'), ...
                ['the old strict-greater stale-index check must fold into ' ...
                'the >= guard, not survive alongside it']);

            %the two DoD helpers key their event-column check on the index
            for name = {'psyrat_get_dod_events', 'psyrat_dodcontraststatus'}
                body = testCase.functionBody(name{1});
                testCase.verifyEqual(numel(strfind(body, ...
                    'psyrat_prefs.proc.inp.event == length(psyrat_data.proc.collist)')), 1, ...
                    sprintf(['%s must detect a missing event column by the ' ...
                    'sentinel index, not the name ''none'' (B25)'], name{1}));
            end
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

        function body = functionBody(testCase, name)
            %extract a local function's body: from its declaration line to
            %the next top-level function declaration
            src = fileread(testCase.startprocSource());
            a = regexp(src, ['\nfunction [^\n]*' name '\('], 'once');
            testCase.assertNotEmpty(a, sprintf( ...
                '%s not found in psyrat_startproc.m; the anchor has drifted', ...
                name));
            rest = src(a + 1:end);
            nextfun = regexp(rest(2:end), '\nfunction ', 'once');
            if isempty(nextfun)
                body = rest;
            else
                body = rest(1:nextfun);
            end
        end

    end
end
