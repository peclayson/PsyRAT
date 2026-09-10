classdef TestStartprocColumnAutodetect < PsyRATTestBase
    % Column auto-detect preselection in psyrat_startproc: the Measurement
    % fallback must never land on a facet/role column.
    %
    % WHAT WENT WRONG. psyrat_startproc_gui seeds positional defaults
    % (id = 1, meas = 2, everything else = 'none') and then improves them by
    % matching generic header names. The Measurement candidate list is
    % {'Measurement','Meas'}. The then-shipped example file (the dev-only
    % fixture testdata.csv; the shipped stand-in with the same header row is
    % now legacy_headers_fixture.csv) has headers subjid,group,event,ern: its score
    % column 'ern' matches neither candidate, so Measurement stayed on the
    % positional default of column 2, which is the GROUP column. Two outcomes,
    % both bad:
    %   (a) 'group' also matched the Group candidate list, so a single column was
    %       preselected into two roles. The run then died at the Run button on
    %       the measurement/facet overlap dialog. The shipped example file was
    %       unrunnable exactly as the GUI preselected it.
    %   (b) Give the same file a numerically coded facet column that matches no
    %       literal (say 'grp' holding 1/2). Group stays 'none', so nothing
    %       overlaps, the numeric-measurement guard passes, and the run silently
    %       reports the reliability of a group code. That is the dangerous case:
    %       a plausible-looking coefficient for a variable that is not a score.
    %
    % A 2026-07-31 fix (CC8) corrected the candidate LITERALS only and recorded
    % that the matching logic was untouched. Literals cannot help a column named
    % 'ern'. The fallback is what changed here: when nothing matches, take the
    % last NUMERIC column that no other role already claimed.
    %
    % HOW THIS IS TESTED. The auto-detect block lives inside psyrat_startproc_gui,
    % a local function that builds a uicontrol figure immediately afterward, so it
    % cannot be called headlessly. Rather than hand-copy the logic into this file
    % (the mirror pattern used by TestStartprocPreflightAgreement, which can drift
    % from the source -- drift is exactly what caused RC-24), these tests EXTRACT
    % the block verbatim from psyrat_startproc.m at run time and execute it. The
    % code under test is therefore the shipped code, and a future edit to the
    % block is exercised here automatically.
    %
    % This pins the preselection only. It cannot verify what the user actually
    % sees in the Measurement popup; that still needs a live GUI click-through.

    methods (Test)

        function testShippedExampleDoesNotPreselectTheGroupColumn(testCase)
            % The canonical regression, run against a real shipped file
            % through the real reader so the column classes are whatever the
            % toolbox loads. legacy_headers_fixture.csv reproduces the header
            % row that originally triggered the bug (subjid,group,event,ern;
            % see the class comment) and, unlike the retired dev-only
            % testdata.csv, ships in public releases, so this regression stays
            % runnable from a release tree.
            root = testCase.projectRoot();
            tbl = psyrat_readtable('file', ...
                fullfile(root, 'test_data', 'legacy_headers_fixture.csv'));

            inp = localRunAutodetect(root, tbl);
            names = tbl.Properties.VariableNames;

            testCase.verifyEqual(names{inp.id}, 'subjid');
            testCase.verifyEqual(names{inp.group}, 'group');
            testCase.verifyEqual(names{inp.event}, 'event');

            % (a) Measurement is not the group column, and not any other role.
            testCase.verifyNotEqual(inp.meas, inp.group, ...
                'Measurement is still preselected to the group column');
            testCase.verifyFalse(any(inp.meas == [inp.id inp.group inp.event]), ...
                'Measurement is preselected to a column already used as a facet');

            % ...it is the score column.
            testCase.verifyEqual(names{inp.meas}, 'ern');
        end

        function testFallbackSkipsANumericallyCodedFacetThatMatchedNoLiteral(testCase)
            % Failure mode (b): the silent one. 'grp' holds 1/2 and matches no
            % candidate name, so Group stays 'none' and no overlap guard fires.
            % The old positional default would have picked it.
            root = testCase.projectRoot();
            tbl = table([1;1;2;2], [1;1;2;2], {'cor';'err';'cor';'err'}, ...
                [-2.8; -7.8; -3.1; -6.4], ...
                'VariableNames', {'subjid','grp','cond','ern'});

            inp = localRunAutodetect(root, tbl);

            testCase.verifyEqual(inp.group, 5, ...
                'precondition: ''grp'' must not match the Group candidates');
            testCase.verifyEqual(inp.meas, 4, ...
                'Measurement fell back onto the numerically coded facet column');
        end

        function testFallbackNeverPicksANonNumericColumn(testCase)
            % The last column is text, so the fallback has to skip it. Without
            % the numeric filter it would preselect a column the Run-button
            % numeric guard immediately rejects.
            root = testCase.projectRoot();
            tbl = table([1;1;2;2], [-2.8; -7.8; -3.1; -6.4], ...
                {'a';'b';'a';'b'}, ...
                'VariableNames', {'subjid','ern','label'});

            inp = localRunAutodetect(root, tbl);

            testCase.verifyEqual(inp.meas, 2, ...
                'Measurement fell back onto a non-numeric column');
        end

        function testDefaultIsLeftAloneWhenNoColumnQualifies(testCase)
            % Every numeric column is claimed by a role and the rest are text,
            % so there is no defensible measurement column. The fallback must
            % not invent one: leave the positional default so the Run-button
            % guards give the user a specific dialog.
            root = testCase.projectRoot();
            tbl = table([1;1;2;2], {'flk';'flk';'flk';'flk'}, ...
                {'cor';'err';'cor';'err'}, ...
                'VariableNames', {'subjid','group','event'});

            inp = localRunAutodetect(root, tbl);

            testCase.verifyEqual(inp.meas, 2, ...
                'the fallback picked a column when none qualified');
        end

        function testNameMatchingStillWinsOverTheFallback(testCase)
            % The fallback must not override a real header match. In each case
            % the matched column is NOT the last numeric column, so a fallback
            % that ran anyway would give a different answer.
            root = testCase.projectRoot();
            for nm = {'Measurement','Meas','meas','MEAS'}
                tbl = table([1;1;2;2], [-2.8; -7.8; -3.1; -6.4], ...
                    [410; 395; 502; 488], ...
                    'VariableNames', {'subjid', nm{1}, 'rt'});

                inp = localRunAutodetect(root, tbl);

                testCase.verifyEqual(inp.meas, 2, sprintf( ...
                    'header ''%s'' no longer auto-detects as Measurement', nm{1}));
            end
        end

        function testOccasionAndOtherRoleMatchingIsPreserved(testCase)
            % CC8 restored 'Session' (it used to carry a leading space against
            % an exact strcmpi). Re-pin the whole occasion list, plus the fact
            % that the new fallback then routes Measurement around it.
            root = testCase.projectRoot();
            for nm = {'Time','Occasion','Session','session'}
                tbl = table([1;1;2;2], [1;2;1;2], [-2.8; -7.8; -3.1; -6.4], ...
                    'VariableNames', {'subjid', nm{1}, 'ern'});

                inp = localRunAutodetect(root, tbl);

                testCase.verifyEqual(inp.time, 2, sprintf( ...
                    'header ''%s'' no longer auto-detects as Occasion', nm{1}));
                % ...and the occasion column, though numeric, is claimed, so
                % Measurement must skip it.
                testCase.verifyEqual(inp.meas, 3, sprintf( ...
                    'Measurement collided with the ''%s'' occasion column', nm{1}));
            end

            % Group and Event candidate literals.
            tbl = table([1;1;2;2], {'ctl';'ctl';'pat';'pat'}, ...
                {'cor';'err';'cor';'err'}, [-2.8; -7.8; -3.1; -6.4], ...
                'VariableNames', {'subjid','Group','Type','ern'});
            inp = localRunAutodetect(root, tbl);
            testCase.verifyEqual(inp.group, 2, '''Group'' no longer auto-detects');
            testCase.verifyEqual(inp.event, 3, '''Type'' no longer auto-detects');
            testCase.verifyEqual(inp.meas, 4);

            % Participant candidate literals still win column 1's role.
            for nm = {'Subject','ID','Participant','SubjID','Subj'}
                tbl = table([1;1;2;2], [-2.8; -7.8; -3.1; -6.4], ...
                    'VariableNames', {nm{1}, 'ern'});
                inp = localRunAutodetect(root, tbl);
                testCase.verifyEqual(inp.id, 1, sprintf( ...
                    'header ''%s'' no longer auto-detects as Participant ID', nm{1}));
                testCase.verifyEqual(inp.meas, 2);
            end
        end

        function testMeasurementIsNeverUnsetAcrossEveryCase(testCase)
            % Downstream sites assume inp.meas is a valid index into
            % raw.colnames: psyrat_measurepopupdisplay builds the popup Value
            % from it, and psyrat_coercemeasureindices maps an empty value to 1
            % (the ID column), so an unset default would silently preselect the
            % participant key. Sweep the fixtures and confirm the fallback always
            % leaves a usable scalar index.
            root = testCase.projectRoot();
            fixtures = { ...
                table([1;2], [3;4], 'VariableNames', {'subjid','ern'}), ...
                table([1;2], {'a';'b'}, 'VariableNames', {'subjid','lab'}), ...
                table({'a';'b'}, {'c';'d'}, 'VariableNames', {'x','y'}), ...
                table([1;2], 'VariableNames', {'onlycol'})};

            for k = 1:numel(fixtures)
                inp = localRunAutodetect(root, fixtures{k});
                ncol = width(fixtures{k});
                testCase.verifyNotEmpty(inp.meas, sprintf( ...
                    'fixture %d left Measurement unset', k));
                testCase.verifyTrue(isscalar(inp.meas) && inp.meas >= 1, ...
                    sprintf('fixture %d left a non-scalar/invalid Measurement', k));
                % A one-column file cannot satisfy the seeded default of 2;
                % psyrat_coercemeasureindices clamps that downstream, which is
                % the pre-existing behavior and is left unchanged here.
                if ncol > 1
                    testCase.verifyLessThanOrEqual(inp.meas, ncol, sprintf( ...
                        'fixture %d left Measurement past the last column', k));
                end
            end
        end

    end
end

function inp = localRunAutodetect(root, tbl)
%Execute the REAL auto-detect block from psyrat_startproc.m against a table.
%
%The block is a local function's straight-line body, so it cannot be called.
%Instead it is lifted verbatim out of the source: it starts at the unique line
%  if ~isfield(psyrat_prefs.proc,'inp')
%and ends at the first following line that is exactly 'end' at zero indent (every
%line inside the block is indented). Executing the extracted text means this test
%exercises the shipped code rather than a hand-copied mirror that could drift.
%
%The block reads only psyrat_prefs.proc, psyrat_data.proc.collist and
%psyrat_data.raw, so those are all that has to be built here. It is written to a
%temporary function file rather than eval'd so the variables it creates land in a
%normal function workspace.

src = fullfile(root, 'subroutines', 'guis', 'psyrat_startproc.m');
%deblank (not strtrim) strips trailing whitespace, including a CR on a
%CRLF checkout, while keeping leading indentation significant: the block's
%closing 'end' is the first one at zero indent.
lines = deblank(strsplit(fileread(src), newline));

startix = find(strcmp(strtrim(lines), 'if ~isfield(psyrat_prefs.proc,''inp'')'));
assert(isscalar(startix), ...
    ['Could not uniquely locate the auto-detect block in %s. ' ...
    'If the block was renamed or moved, update this extractor.'], src);

endoffset = find(strcmp(lines(startix+1:end), 'end'), 1);
assert(~isempty(endoffset), ...
    'Could not find the closing end of the auto-detect block in %s.', src);
blocktext = strjoin(lines(startix:startix+endoffset), newline);

%build the two structures the block consumes
colnames = tbl.Properties.VariableNames;
psyrat_data = struct();
psyrat_data.raw.colnames = colnames;
psyrat_data.raw.data = tbl;
psyrat_data.proc.collist = [colnames, {'(none)'}];
psyrat_prefs = struct();
psyrat_prefs.proc = struct();

%wrap the block in a throwaway function and run it. The function name is unique
%per call so MATLAB's function cache can never serve a stale probe from an
%earlier test method.
tmpdir = tempname;
mkdir(tmpdir);
cleaner = onCleanup(@() localCleanup(tmpdir)); %#ok<NASGU>
[~, tmpstem] = fileparts(tmpdir);
fname = matlab.lang.makeValidName(['psyrat_autodetect_probe_' tmpstem]);
fid = fopen(fullfile(tmpdir, [fname '.m']), 'w');
fprintf(fid, 'function psyrat_prefs = %s(psyrat_prefs,psyrat_data)\n', fname);
fprintf(fid, '%s\n', blocktext);
fprintf(fid, 'end\n');
fclose(fid);

addpath(tmpdir);
rehash path;
psyrat_prefs = feval(fname, psyrat_prefs, psyrat_data);
inp = psyrat_prefs.proc.inp;
end

function localCleanup(tmpdir)
%remove the probe from the path and from disk
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(tmpdir);
warning('on', 'MATLAB:rmpath:DirNotFound');
rmdir(tmpdir, 's');
end
