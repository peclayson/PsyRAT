classdef TestStartprocPreflightAgreement < PsyRATTestBase
    % RC-24 regression: the GUI's early sserrvar+time guard must agree with
    % psyrat_preflight_validate. Since G34 (2026-08-29) the class also mirrors
    % the SECOND early guard, the diffest+time check, whose lack of a mirror is
    % how the same drift class recurred there.
    %
    % psyrat_startproc runs its own copy of the "subject-level error variance
    % with an occasion facet" check so the user gets a dialog before the data
    % file is loaded. psyrat_preflight_validate runs the authoritative version
    % later, just before CmdStan launches. Because the GUI copy returns first,
    % ANY combination the GUI rejects is unreachable, whatever preflight says.
    %
    % That is precisely how RC-24 happened: commit 8ae949f added the analysis-25
    % engine (trt_sserrvar) and relaxed preflight to exempt it, but never touched
    % psyrat_startproc, whose guard had been inherited from the ERA Toolbox port
    % ~4.5 months earlier. The two conditions silently disagreed and the stale
    % one won, leaving a fully implemented design unreachable from the GUI.
    %
    % psyrat_startproc is a long GUI callback function and cannot be invoked
    % headlessly, so this test pins the PREDICATE rather than the dialog: it
    % re-states the guard's condition as a local function that must be kept
    % character-identical to the source, and sweeps it against preflight over
    % the (diffest x dynrel x sserrvar x hasTime) grid. A future edit to either
    % side that breaks the agreement fails here.
    %
    % Scope: the sweep skips dynrel==2 with diffest==3. That combination is
    % governed by a DIFFERENT, earlier startproc gate -- DoD under dynamic
    % reliability is permitted ONLY as the person-specific nonconcurrent design
    % (analyses 28/29: subject-level reliability on, family gamma) and rejected
    % otherwise -- and preflight examines it under its own dod-specific checks,
    % not the sserrvar+time predicate under test here. Including the cell would
    % therefore compare two unrelated predicates. Note the GUI gate is stricter
    % than preflight in this cell: preflight enforces only sserrvar==2 for
    % DoD+dynrel, and only when a time facet is present (the shared
    % preflight:timeWithDiffUnsupported check exempts dynrel==2 && sserrvar==2),
    % and it never checks the family, whereas the GUI gate requires ssubjrel==2
    % unconditionally (the family requirement this note once described was
    % dropped when the Gaussian 28/29 arm landed; the gate's own comment says
    % the family is no longer part of it). That asymmetry is recorded under B17
    % in PSYRAT_AUDIT_FINDINGS.md, and this note previously said the
    % combination was rejected outright, which has been false since the gate
    % learned the 28/29 design.
    %
    % If this test fails, do NOT relax it. Reconcile the two conditions -- the
    % disagreement itself is the defect.

    methods (Test)

        function testGuardAgreesWithPreflightOnTheFullGrid(testCase)
            errId = 'preflight:timeWithSubjectReliabilityUnsupported';

            for diffest = [1 2 3]
                for dynrel = [1 2]
                    for sserrvar = [1 2]
                        for hasTime = [false true]

                            % Skip the DoD-under-dynrel cell: an EARLIER
                            % startproc gate governs it -- permitted only as
                            % the person-specific nonconcurrent design
                            % (analyses 28/29: subject-level reliability on,
                            % family gamma), rejected otherwise with its own
                            % dialog (grep the gate's 'difference-of-
                            % differences (DoD, diffest == 3)' comment in
                            % psyrat_startproc.m rather than a line number).
                            % That gate is a different concern from the
                            % sserrvar+time predicate under test here, and
                            % preflight examines the combination under its own
                            % dod/diff error ids (timeWithDiffUnsupported,
                            % dodWithSubjectReliabilityUnsupported), never the
                            % id this sweep asserts on.
                            if dynrel == 2 && diffest == 3
                                continue;
                            end

                            guardBlocks = localStartprocGuard( ...
                                hasTime, sserrvar, dynrel, diffest);

                            tbl = localMakeTable(hasTime);
                            report = psyrat_preflight_validate( ...
                                'datatable', tbl, ...
                                'sserrvar', sserrvar, ...
                                'diffest', diffest, ...
                                'dynrel', dynrel);
                            preflightBlocks = any(strcmp( ...
                                {report.errors.id}, errId));

                            testCase.verifyEqual(guardBlocks, preflightBlocks, ...
                                sprintf(['startproc guard and preflight disagree at ' ...
                                'diffest=%d dynrel=%d sserrvar=%d hasTime=%d ' ...
                                '(guard blocks=%d, preflight blocks=%d)'], ...
                                diffest, dynrel, sserrvar, hasTime, ...
                                guardBlocks, preflightBlocks));
                        end
                    end
                end
            end
        end

        function testAnalysis25IsAllowedByBothGates(testCase)
            % The specific combination RC-24 was about: subject-level error
            % variance + an occasion facet, non-difference and non-dynamic.
            % Neither gate may block it.
            testCase.verifyFalse(localStartprocGuard(true, 2, 1, 1), ...
                'the startproc guard still blocks trt_sserrvar (analysis 25)');

            report = psyrat_preflight_validate( ...
                'datatable', localMakeTable(true), ...
                'sserrvar', 2, 'diffest', 1, 'dynrel', 1);
            testCase.verifyFalse(any(strcmp({report.errors.id}, ...
                'preflight:timeWithSubjectReliabilityUnsupported')), ...
                'preflight still blocks trt_sserrvar (analysis 25)');
        end

        function testDifferenceDesignsWithTimeAreStillBlocked(testCase)
            % The guard must not become a no-op: non-dynamic DIFFERENCE designs
            % with a time facet stay rejected by both gates.
            for diffest = [2 3]
                testCase.verifyTrue(localStartprocGuard(true, 2, 1, diffest), ...
                    sprintf('guard no longer blocks diffest=%d with a time facet', diffest));
            end
        end

        function testDiffGuardAgreesWithPreflightOnTheFullGrid(testCase)
            % The SECOND startproc gate (the diffest + occasion check) had no
            % mirror here, and that gap is exactly how G34 happened: preflight
            % learned the two-facet difference pipeline (analysis 10) and this
            % gate never did, so a shipped design was unreachable from the GUI
            % for ~7 months with every suite green. Same sweep, same skip, the
            % diff-specific preflight id.
            errId = 'preflight:timeWithDiffUnsupported';

            for diffest = [1 2 3]
                for dynrel = [1 2]
                    for ssubjrel = [1 2]
                        for hasTime = [false true]

                            % Same scope cut as the first sweep: DoD under
                            % dynamic reliability is governed by an earlier,
                            % stricter startproc gate (analyses 28/29 only),
                            % so the two predicates are not comparable there.
                            if dynrel == 2 && diffest == 3
                                continue;
                            end

                            guardBlocks = localStartprocDiffGuard( ...
                                hasTime, dynrel, diffest, ssubjrel);

                            tbl = localMakeTable(hasTime);
                            report = psyrat_preflight_validate( ...
                                'datatable', tbl, ...
                                'sserrvar', ssubjrel, ...
                                'diffest', diffest, ...
                                'dynrel', dynrel);
                            preflightBlocks = any(strcmp( ...
                                {report.errors.id}, errId));

                            testCase.verifyEqual(guardBlocks, preflightBlocks, ...
                                sprintf(['startproc diff guard and preflight disagree at ' ...
                                'diffest=%d dynrel=%d ssubjrel=%d hasTime=%d ' ...
                                '(guard blocks=%d, preflight blocks=%d)'], ...
                                diffest, dynrel, ssubjrel, hasTime, ...
                                guardBlocks, preflightBlocks));
                        end
                    end
                end
            end
        end

        function testAnalysis10IsAllowedByBothGates(testCase)
            % The specific combination G34 was about: the static two-facet
            % difference score (diffest=2, subject-level off, occasion facet).
            % Preflight calls it "the supported two-facet (test-retest)
            % difference pipeline"; neither gate may block it.
            testCase.verifyFalse(localStartprocDiffGuard(true, 1, 2, 1), ...
                'the startproc diff guard still blocks analysis 10');
            testCase.verifyFalse(localStartprocGuard(true, 1, 1, 2), ...
                'the startproc sserrvar guard blocks analysis 10');

            report = psyrat_preflight_validate( ...
                'datatable', localMakeTable(true), ...
                'sserrvar', 1, 'diffest', 2, 'dynrel', 1);
            testCase.verifyFalse(any(strcmp({report.errors.id}, ...
                'preflight:timeWithDiffUnsupported')), ...
                'preflight blocks analysis 10');
        end

        function testDoDWithTimeStaysBlockedByTheDiffGuard(testCase)
            % The relaxed gate must not become a no-op: non-dynamic DoD with a
            % time facet stays rejected whatever the subject-level setting.
            for ssubjrel = [1 2]
                testCase.verifyTrue(localStartprocDiffGuard(true, 1, 3, ssubjrel), ...
                    sprintf(['diff guard no longer blocks DoD with a time ' ...
                    'facet at ssubjrel=%d'], ssubjrel));
            end
        end

        function testMirrorPredicatesAppearInTheSource(testCase)
            % THE WELD. The two sweeps above compare test-file MIRRORS against
            % preflight and never execute psyrat_startproc, so on their own
            % they would stay green if the real gate were edited or reverted
            % while the mirrors were not -- the exact drift class this file
            % exists to catch (G34's adversarial review demonstrated it). This
            % census ties each mirror's load-bearing predicate text to the
            % source: if either gate's condition changes, the count moves and
            % this fails, forcing the mirror and the sweep to be re-derived
            % together. A text weld, not an execution check -- the dialog halves
            % remain click-through territory.
            src = fileread(fullfile(fileparts(mfilename('fullpath')), ...
                '..', 'subroutines', 'guis', 'psyrat_startproc.m'));

            % gate 1 (ssubjrel + time), the middle line of its condition
            testCase.verifyEqual(numel(strfind(src, ...
                ['psyrat_prefs.proc.ssubjrel == 2 && ' ...
                'psyrat_prefs.proc.dynrel ~= 2 && ...'])), 1, ...
                ['the ssubjrel gate''s condition text moved or multiplied; ' ...
                're-derive localStartprocGuard against the source']);

            % gate 2 (diffest + time), the G34 exemption term
            testCase.verifyEqual(numel(strfind(src, ...
                '~(psyrat_prefs.proc.diffest == 2 && psyrat_prefs.proc.ssubjrel == 1)')), 1, ...
                ['the diffest gate''s analysis-10 exemption moved or ' ...
                'multiplied; re-derive localStartprocDiffGuard against the ' ...
                'source']);
        end

    end
end

function blocks = localStartprocGuard(hasTime, sserrvar, dynrel, diffest)
%Mirror of the guard at subroutines/guis/psyrat_startproc.m (the sserrvar +
%occasion check). KEEP THIS CONDITION CHARACTER-IDENTICAL TO THE SOURCE.
%In the source, hasTime is keyed on the sentinel INDEX, not on its label:
%  psyrat_prefs.proc.inp.time ~= length(psyrat_data.proc.collist)
%This comment previously quoted a name-keyed strcmpi(...,'none') form that the
%source had already dropped under B25. The index form is what to mirror, and
%it is also why relabelling the sentinel to '(none)' changed nothing here.
blocks = hasTime && sserrvar == 2 && dynrel ~= 2 && any(diffest == [2 3]);
end

function blocks = localStartprocDiffGuard(hasTime, dynrel, diffest, ssubjrel)
%Mirror of the SECOND guard in psyrat_startproc (the diffest + occasion check
%inside the `if any(psyrat_prefs.proc.diffest == [2 3])` block). KEEP THIS
%CONDITION CHARACTER-IDENTICAL TO THE SOURCE. hasTime is keyed on the sentinel
%index there, exactly as in localStartprocGuard.
blocks = any(diffest == [2 3]) && hasTime && dynrel ~= 2 && ...
    ~(diffest == 2 && ssubjrel == 1);
end

function tbl = localMakeTable(hasTime)
%Minimal valid long-format table, with or without a time/occasion column.
ids  = {'s1'; 's1'; 's2'; 's2'};
meas = [1.0; 1.2; 0.9; 1.1];
if hasTime
    tbl = table(ids, meas, {'t1';'t2';'t1';'t2'}, ...
        'VariableNames', {'id','meas','time'});
else
    tbl = table(ids, meas, 'VariableNames', {'id','meas'});
end
end
