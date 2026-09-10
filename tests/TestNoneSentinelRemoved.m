classdef TestNoneSentinelRemoved < PsyRATTestBase
    %Pins finding B34 -- psyrat_relsummary must NEVER store the CHAR 'none' in
    %the good/bad id carrier. Producers now assign cell(0,1).
    %
    %WHY THE CHAR HAD TO GO. The two Good IDs exporters
    %(psyrat_depoverallt / psyrat_trt_reloverallt) splice the carrier with
    %"gids = [gids; group(j).goodids(:)]" and do NO type test, so the char
    %escaped into user-facing output in three different ways depending on how
    %many groups had nothing retained:
    %  - SOME groups affected: MATLAB wraps the 4x1 char as ONE cell element,
    %    so the export silently lists a phantom participant named "none".
    %  - ALL groups affected: the splice yields a CHAR, so the XLSX branch
    %    (char(gids(i))) writes four junk one-letter rows per group and the CSV
    %    branch (gids{i}) throws MATLAB:cellRefFromNonCell.
    %It also crashed psyrat_relsummary itself: table('none') errors with
    %MATLAB:table:parseArgs:WrongNumberArgs because MATLAB reads a char row
    %vector as variable NAMES.
    %
    %The B29 arm tests are NOT what fixed this. Those choose which local
    %variable to build; the STORED field kept the char regardless, and the
    %exporters read the stored field.

    methods (Test)

        function testSubjectLevelAllExcludedStoresEmptyCellNotChar(testCase)
            %Route A from the B34 ledger entry: sing_sserr where the cutoff
            %excludes everyone. relerr.nogooddata is never set in this branch,
            %so nothing aborts and the carrier reaches the exporters.
            out = testCase.runSserr(0.9999);
            grp = out.relsummary.group;
            testCase.assertNotEmpty(grp, 'Fixture produced no groups.');
            for g = 1:numel(grp)
                testCase.verifyFalse(ischar(grp(g).goodids), sprintf( ...
                    'group(%d).goodids is a CHAR; the sentinel is back.', g));
                testCase.verifyTrue(iscell(grp(g).goodids), sprintf( ...
                    'group(%d).goodids must be a cell array.', g));
                testCase.verifyEmpty(grp(g).goodids, sprintf( ...
                    'group(%d).goodids must be EMPTY when nothing cleared.', g));
            end
        end

        function testSubjectLevelReachableCutoffStillRetainsIds(testCase)
            %POSITIVE CONTROL. Without it the test above would pass on a
            %fixture that never reached the site. At a reachable cutoff the
            %same fixture must retain real ids and report a non-zero goodn.
            out = testCase.runSserr(0.01);
            grp = out.relsummary.group;
            for g = 1:numel(grp)
                testCase.verifyTrue(iscell(grp(g).goodids), 'carrier type');
                testCase.verifyNotEmpty(grp(g).goodids, sprintf( ...
                    'group(%d) should retain ids at a reachable cutoff.', g));
                testCase.verifyGreaterThan(grp(g).event(1).goodn, 0, ...
                    'goodn must be positive when ids were retained.');
            end
        end

        function testExportSpliceStaysACellAndAddsNoPhantomRows(testCase)
            %The consequence that actually reached users. Replicates the
            %exporters' splice expression verbatim (psyrat_depoverallt.m and
            %psyrat_trt_reloverallt.m) against the REAL carrier, and asserts
            %it stays a cell and contributes NOTHING when nothing was retained.
            out = testCase.runSserr(0.9999);
            grp = out.relsummary.group;
            gids = [];
            for j = 1:numel(grp)
                gids = [gids; grp(j).goodids(:)]; %#ok<AGROW>
            end
            testCase.verifyFalse(ischar(gids), ...
                'Splice produced a CHAR: the XLSX export would emit junk rows.');
            testCase.verifyTrue(iscell(gids) || isempty(gids), ...
                'Splice must stay a cell array.');
            testCase.verifyEmpty(gids, ...
                'No participant cleared the cutoff, so no Good IDs rows.');
        end

        function testMixedGroupsExportNoParticipantNamedNone(testCase)
            %THE SILENT MODE, and the one that actually reached users. When
            %only SOME groups have nothing retained, the old char did NOT
            %produce a char splice: MATLAB wrapped the 4x1 char as ONE cell
            %element, so the export listed a phantom participant rendering as
            %the word "none" -- in BOTH the XLSX and CSV branches, with no
            %error and nothing visibly wrong.
            %
            %depcutoff 0.95 is a GENUINE mixed state on this fixture, measured
            %rather than constructed: group 1 retains nothing, group 2 retains
            %two ids. Building the mixed shape by hand instead would make this
            %assertion tautological -- it would pass on the unfixed source too.
            out = testCase.runSserr(0.95);
            grp = out.relsummary.group;
            nretained = arrayfun(@(x) numel(x.goodids), grp);
            testCase.assertTrue(any(nretained == 0) && any(nretained > 0), ...
                sprintf(['Fixture no longer produces a MIXED state at 0.95 ' ...
                '(retained per group = %s); this test proves nothing until ' ...
                'a mixed cutoff is restored.'], mat2str(nretained)));

            gids = [];
            for j = 1:numel(grp)
                gids = [gids; grp(j).goodids(:)]; %#ok<AGROW>
            end
            testCase.verifyTrue(iscell(gids), ...
                'Mixed splice must stay a cell array.');
            testCase.verifyNotEmpty(gids, ...
                'The retaining group must still contribute its ids.');
            for i = 1:numel(gids)
                testCase.verifyClass(gids{i}, 'char', ...
                    'each exported id must be a char row vector');
                testCase.verifyNotEqual(gids{i}, 'none', sprintf( ...
                    ['Exported id %d is the literal ''none'': the sentinel ' ...
                    'leaked back in as a phantom participant.'], i));
            end
        end

        function testNoProducerAssignsTheCharSentinel(testCase)
            %Source pin. Separates conversion from deletion: a whole-file zero
            %count on every spelling that assigned the char to the carrier.
            src = fileread(testCase.relsummarySource());
            pats = { ...
                'goodids\s*=\s*''none''', ...
                'tempids\{1\}\s*=\s*''none''', ...
                'eventgoodids\s*=\.\.\.\s*\n\s*''none''' };
            for k = 1:numel(pats)
                n = numel(regexp(src, pats{k}, 'match'));
                testCase.verifyEqual(n, 0, sprintf( ...
                    'Producer spelling %d is back in the source (%d hits).', k, n));
            end
            %and the replacement really is present, so this cannot pass by the
            %producers having been deleted outright
            nrep = numel(regexp(src, 'cell\(0,1\)', 'match'));
            testCase.verifyGreaterThanOrEqual(nrep, 33, sprintf( ...
                'Expected >= 33 cell(0,1) producers, found %d.', nrep));
        end

    end

    methods (Access = private)

        function out = runSserr(~, depcutoff)
            d = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
            d.proc.interactive = false;
            out = psyrat_relsummary('psyrat_data', d, 'analysis', 'sing_sserr', ...
                'depcutoff', depcutoff, 'meascutoff', 2, 'depcentmeas', 1, ...
                'gcoeff', 1, 'CI', 0.95);
        end

        function p = relsummarySource(~)
            %absolute path from THIS file's location, never from pwd or the
            %MATLAB path (worktree copies have shadowed the toolbox before)
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here), 'subroutines', 'calculation', ...
                'psyrat_relsummary.m');
        end

    end
end
