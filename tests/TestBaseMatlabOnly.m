classdef TestBaseMatlabOnly < PsyRATTestBase
    % Guards the shipped code against functions that do NOT exist in base
    % MATLAB, where a base-MATLAB equivalent is available and should be used.
    %
    % WHY THIS EXISTS. A developer machine usually has the Statistics and
    % Machine Learning Toolbox installed AND other EEG/ERP toolboxes on the path
    % - FieldTrip, for one, ships shims in external/stats. So a call like
    % range() resolves happily there, possibly to FieldTrip's copy rather than
    % to the toolbox at all, and then dies on a user's base MATLAB with
    % MATLAB:UndefinedFunction. That is exactly what happened to
    % TestDynrelSubjNondiffReliability: it went green locally and red on every
    % base-MATLAB run until the calls were rewritten as max-min.
    %
    % A local pass is therefore NOT evidence that the code is base-MATLAB clean;
    % this test is. It reads the source rather than calling the functions, so it
    % reports the same answer no matter what is installed on the machine running
    % it.
    %
    % THERE IS NO CI. The GitHub Actions lane that once ran this on a base-MATLAB
    % runner was retired 2026-08-20 (commit 2fbd5d7) and Actions is disabled at
    % the repo level, so this source scan is now the ONLY base-MATLAB signal that
    % exists. Earlier revisions of this comment cited
    % .github/workflows/matlab-validation.yml as the reason the guarantee held;
    % that file is gone, which makes this test more load-bearing, not less.
    %
    % SCOPE - what belongs on the denylist. Only functions that (a) are absent
    % from base MATLAB and (b) have a plain base-MATLAB equivalent, so the fix is
    % always a rewrite and never a skip. Deliberately excluded:
    %
    %   - quantile: used ~200 times in the calculation code. Base MATLAB R2026a
    %     supplies it; the CI run that caught range() proves this, since every
    %     posterior-CI test would otherwise have failed rather than two.
    %   - normcdf, gaminv, fitlme, hmcSampler, dlarray: genuinely toolbox-gated
    %     with no base equivalent. Tests that need them skip via assumeTrue
    %     instead (tests/TestNativeEngine.m, tests/TestGammaDiffCopula*.m). Do
    %     not add them here - this test would then demand a rewrite that cannot
    %     be written.
    %
    % Adding a denylist entry is cheap; removing one should be argued.

    properties (Constant)
        % {name, base-MATLAB replacement to suggest in the failure message}
        Denylist = {
            'range',      'max(x) - min(x)'
            'nanmean',    'mean(x, ''omitnan'')'
            'nanstd',     'std(x, ''omitnan'')'
            'nanvar',     'var(x, ''omitnan'')'
            'nanmedian',  'median(x, ''omitnan'')'
            'nansum',     'sum(x, ''omitnan'')'
            'nanmax',     'max(x, [], ''omitnan'')'
            'nanmin',     'min(x, [], ''omitnan'')'
            'prctile',    'quantile(x, p/100)'
            'zscore',     '(x - mean(x)) ./ std(x)'
            'iqr',        'diff(quantile(x, [0.25 0.75]))'
            'skewness',   'mean(((x - mean(x)) ./ std(x, 1)).^3)'
            'kurtosis',   'mean(((x - mean(x)) ./ std(x, 1)).^4)'
            'randsample', 'randperm / randi'
            'datasample', 'x(randi(numel(x), n, 1))'
            'tabulate',   'unique + histcounts (or groupcounts)'
            'grpstats',   'accumarray (or groupsummary)'
            };
    end

    methods (Test)

        function testShippedCodeUsesOnlyBaseMatlabFunctions(testCase)
            % Scans every .m file PsyRAT ships or tests with, INCLUDING the
            % vendored code under bundled_dependents/ (B40).
            %
            % That code used to be excluded, on two stated grounds: it is not
            % ours to constrain, and MatlabStan/MatlabProcessManager are loaded
            % only on the CmdStan path. The first is a defensible ownership
            % boundary. The second is true but the inference from it fails -
            % CmdStan is the DEFAULT engine, so "only on the CmdStan path"
            % describes the main path, not a niche one, and is no reason to
            % exempt it from a base-MATLAB guarantee the documentation makes.
            %
            % The exclusion is what let B39 survive a green suite: nansum was
            % already on the denylist below, and the two offending calls sat in
            % bundled_dependents/MatlabStan-2.15.1.0/StanFit.m the whole time,
            % unscanned. Scanning vendored code does not mean owning its style -
            % only its base-MATLAB compatibility, which is the one property the
            % rest of the toolbox's documentation promises on its behalf.
            root = testCase.projectRoot();
            files = localCollectSourceFiles(root);

            % Sanity check on the collector itself, so an empty or mis-rooted
            % file list can never masquerade as a clean tree.
            testCase.assertGreaterThan(numel(files), 50, ...
                ['Collected implausibly few .m files. The scan is probably ', ...
                'rooted in the wrong directory, which would make this test ', ...
                'pass vacuously.']);

            % Positive control on the vendored half specifically. The count
            % assert above cannot catch a silent loss of bundled_dependents/:
            % localCollectSourceFiles walks each root only 'if isfolder(...)',
            % so renaming or moving that folder would drop it with no error and
            % the scan would go green while blind - which is precisely how B39
            % survived a green suite before B40 was fixed. A total-count
            % threshold cannot fail in that direction; this can.
            testCase.assertNotEmpty( ...
                files(contains(files, [filesep 'bundled_dependents' filesep])), ...
                ['No .m files were collected from bundled_dependents/. Either ', ...
                'the folder moved or the root list lost it; in both cases the ', ...
                'base-MATLAB guarantee silently stops covering vendored code ', ...
                '(see B40 in PSYRAT_AUDIT_FINDINGS.md).']);

            offenders = {};
            for k = 1:numel(files)
                offenders = [offenders; ...
                    localScanFile(files{k}, root, TestBaseMatlabOnly.Denylist)]; %#ok<AGROW>
            end

            testCase.verifyEmpty(offenders, localFormatReport(offenders));
        end

        function testScannerDetectsAKnownViolationAndIgnoresLookalikes(testCase)
            % Pins the matcher itself. Without this, a regex typo would silence
            % the test above and it would still report green - the failure mode
            % this whole file exists to prevent, one level up.
            lines = {
                'testCase.verifyGreaterThan(range(T.gen_pt), 1e-6);'  % a real violation
                '    y = nanmean( x );'                               % whitespace before (
                '%the grid range (per group) is centered'             % prose in a comment
                '    z = psyrat_range(x);'                            % name is a suffix
                '    w = obj.range(x);'                               % method call
                '    v = arange;'                                     % no call parens
                '    u = max(x) - min(x);'                            % the correct form
                };
            hits = localScanLines(lines, TestBaseMatlabOnly.Denylist);

            testCase.verifyEqual(numel(hits), 2, ...
                'Scanner must flag exactly the two real violations.');
            testCase.verifyEqual(hits{1}.line, 1);
            testCase.verifyEqual(hits{1}.name, 'range');
            testCase.verifyEqual(hits{2}.line, 2);
            testCase.verifyEqual(hits{2}.name, 'nanmean');
        end

        function testCommentStripperKeepsCodeAndDropsComments(testCase)
            % The stripper is the part most likely to go subtly wrong, and a
            % stripper that ate whole lines would make the scan silently
            % vacuous. Pinned separately.
            testCase.verifyEqual(localStripComment('a = range(x); % note'), ...
                'a = range(x); ');
            % A whole-line comment strips to nothing. verifyEmpty rather than
            % verifyEqual(...,''): the cut returns a 1x0 char, not a 0x0 one,
            % and the callers only ever strtrim the result.
            testCase.verifyEmpty(localStripComment('% range(x) in prose'));
            % A percent sign inside a quoted string is not a comment. This case
            % is why the stripper counts quotes instead of cutting at the first %.
            testCase.verifyEqual(localStripComment('fprintf(''%d rows\n'', n);'), ...
                'fprintf(''%d rows\n'', n);');
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function files = localCollectSourceFiles(root)
%LOCALCOLLECTSOURCEFILES .m files that PsyRAT ships or tests with.
%   subroutines/ and tests/ are walked recursively; the project root is taken
%   non-recursively so the entry points (psyrat_start, psyrat_run, the run*
%   suites) are covered without pulling in scratch or artifact folders.

files = {};

roots = {fullfile(root, 'subroutines'), fullfile(root, 'tests'), ...
    fullfile(root, 'bundled_dependents')};
for k = 1:numel(roots)
    if isfolder(roots{k})
        found = dir(fullfile(roots{k}, '**', '*.m'));
        files = [files; localFullPaths(found)]; %#ok<AGROW>
    end
end

topLevel = dir(fullfile(root, '*.m'));
files = [files; localFullPaths(topLevel)];

% Agent worktree copies are not ours to scan.
% The worktree exclusion mirrors the addpath hygiene the repo already relies
% on: .claude/worktrees/* holds full copies of the toolbox, and scanning them
% would report the same finding many times over. The match is made against
% the path BELOW the scan root, not the absolute path: when the suite itself
% runs from a checkout under .claude/worktrees/ (the usual agent setup),
% every absolute path contains '/.claude/', and an unscoped match would
% empty the list and trip the anti-vacuity assert in
% testShippedCodeUsesOnlyBaseMatlabFunctions. projectRoot() returns no
% trailing filesep, so erase() leaves the below-root portion with its
% leading filesep and the match pattern itself is unchanged.
%
% A bundled_dependents/ skip filter used to sit here. It was DEAD CODE, and
% removing it is not what changed the scan's coverage: no path this function
% collected could ever contain '/bundled_dependents/', because that folder was
% not among the roots walked above and the project root is taken
% non-recursively. Adding it to roots is the change; the filter's removal
% merely stops the next reader from believing it was ever the gate.
files = files(~contains(erase(files, root), [filesep '.claude' filesep]));

% This file excludes itself: its scanner unit tests carry fixture strings
% that deliberately contain range( and nanmean( so the matcher can be shown
% to catch them. Scanning it would report those fixtures as violations. The
% file is not left unchecked as a result - testScannerDetects... and
% testCommentStripper... are what check it, and more directly than a text
% scan could.
files = files(~endsWith(files, [filesep 'TestBaseMatlabOnly.m']));
end

function paths = localFullPaths(entries)
%LOCALFULLPATHS dir() struct array to a cellstr column of absolute paths.
paths = cell(numel(entries), 1);
for k = 1:numel(entries)
    paths{k} = fullfile(entries(k).folder, entries(k).name);
end
end

function offenders = localScanFile(file, root, denylist)
%LOCALSCANFILE Reports denylisted calls in one file as file:line records.
text = fileread(file);
lines = regexp(text, '\r\n|\n|\r', 'split');

hits = localScanLines(lines, denylist);

relative = erase(file, [root filesep]);
offenders = cell(numel(hits), 1);
for k = 1:numel(hits)
    hits{k}.file = relative;
    offenders{k} = hits{k};
end
end

function hits = localScanLines(lines, denylist)
%LOCALSCANLINES Core matcher, shared by the file scan and its own unit test.
%   A hit is the denylisted name used as a CALL: followed by '(', and not
%   preceded by a word character or a dot, so psyrat_range( and obj.range( are
%   left alone. Comments are stripped first, because the toolbox legitimately
%   uses several of these words as English (for example "the grid range (per
%   group)" appears four times in psyrat_computevarcomp.m).

hits = {};
for ln = 1:numel(lines)
    code = localStripComment(lines{ln});
    if isempty(strtrim(code))
        continue
    end
    for d = 1:size(denylist, 1)
        name = denylist{d, 1};
        pattern = ['(?<![\w.])' name '\s*\('];
        if ~isempty(regexp(code, pattern, 'once'))
            hits{end+1} = struct('line', ln, 'name', name, ...
                'replacement', denylist{d, 2}, 'text', strtrim(lines{ln}), ...
                'file', ''); %#ok<AGROW>
        end
    end
end
hits = hits(:);
end

function code = localStripComment(line)
%LOCALSTRIPCOMMENT Drops a trailing MATLAB line comment.
%   Cuts at the first '%' that has an EVEN number of single quotes before it,
%   so a '%' inside a quoted string (a sprintf format, most often) is kept.
%   The transpose operator also spends a quote, so a line carrying BOTH a
%   transpose and a comment can be under-stripped; that direction is safe -
%   it can only produce a visible, fixable false positive, never a silent
%   miss. The repository uses no %{ %} block comments, so line-at-a-time
%   handling is sufficient.

code = line;
idx = strfind(line, '%');
for k = 1:numel(idx)
    if mod(count(line(1:idx(k)-1), ''''), 2) == 0
        code = line(1:idx(k)-1);
        return
    end
end
end

function msg = localFormatReport(offenders)
%LOCALFORMATREPORT Actionable failure text: where, what, and the fix.
if isempty(offenders)
    msg = '';
    return
end

parts = cell(numel(offenders), 1);
for k = 1:numel(offenders)
    o = offenders{k};
    parts{k} = sprintf('  %s:%d  %s(...)  ->  use %s\n      %s', ...
        o.file, o.line, o.name, o.replacement, o.text);
end

msg = sprintf(['%d call(s) to functions that do not exist in base MATLAB.\n' ...
    'A user on base MATLAB hits MATLAB:UndefinedFunction on these, even ' ...
    'though they resolve on a developer machine with the Statistics Toolbox ' ...
    'or FieldTrip on the path. This scan is the only signal for that - there ' ...
    'is no CI (retired 2026-08-20).\n\n%s\n\nHow to fix: rewrite using the ' ...
    'base-MATLAB form shown. Do NOT add an assumeTrue skip - that would ' ...
    'delete the coverage in the one environment where it matters.'], ...
    numel(offenders), strjoin(parts, sprintf('\n')));
end
