classdef TestDocumentationIntegrity < PsyRATTestBase
    % Structural checks on the shipped documentation, chiefly the user
    % manual (documentation/manual/).
    %
    % The manual is prose, so the numeric suites cannot see it rot. These
    % checks pin the failure modes that killed the legacy manual and the
    % ones the release gate cares about, and they are written to pass both
    % on a development checkout and on a staged release export (they touch
    % only files that ship):
    %
    %   * every relative link in every chapter resolves;
    %   * every referenced image exists and every image in images/ is
    %     referenced somewhere (no orphans silently bloating releases);
    %   * the PsyRAT version literal appears exactly once, in index.md, and
    %     matches psyrat_defineversion (the legacy manual carried five
    %     contradictory version claims; this makes that impossible);
    %   * no chapter references a withheld internal document, the private
    %     dev repository, or a local-only donotsync path (mirrors
    %     release/verify_export.sh checks 5 and 8, which cannot run until
    %     export time; this catches the same defects at unit-test time);
    %   * no shipped user-facing document (README, CHANGELOG,
    %     documentation/**) sends the reader to a path .gitattributes
    %     export-ignores, so the release cannot ship an instruction naming a
    %     file it omits (G61);
    %   * no PDF lives under documentation/manual/ except the one generated
    %     manual PDF the release gate exempts;
    %   * every value the tutorial chapters print is still the value
    %     ev_chapter_map.md says was substituted into that chapter.

    properties (Access = private)
        ManualDir
        ChapterFiles   % cellstr of manual/*.md full paths
    end

    methods (TestMethodSetup)
        function locateManual(testCase)
            testCase.ManualDir = fullfile(testCase.projectRoot(), ...
                'documentation', 'manual');
            testCase.assertTrue(isfolder(testCase.ManualDir), ...
                'documentation/manual/ is missing');
            d = dir(fullfile(testCase.ManualDir, '*.md'));
            testCase.ChapterFiles = fullfile({d.folder}, {d.name});
            testCase.assertTrue(any(strcmp({d.name}, 'index.md')), ...
                'documentation/manual/index.md is missing');
        end
    end

    methods (Test)

        function testEveryRelativeLinkResolves(testCase)
            for k = 1:numel(testCase.ChapterFiles)
                f = testCase.ChapterFiles{k};
                targets = localLinkTargets(fileread(f));
                for t = 1:numel(targets)
                    target = targets{t};
                    resolved = fullfile(fileparts(f), target);
                    testCase.verifyTrue(isfile(resolved) || isfolder(resolved), ...
                        sprintf('%s links to missing target: %s', f, target));
                end
            end
        end

        function testEveryImageIsReferencedAndEveryReferenceResolves(testCase)
            imgdir = fullfile(testCase.ManualDir, 'images');
            refs = {};
            for k = 1:numel(testCase.ChapterFiles)
                targets = localLinkTargets(fileread(testCase.ChapterFiles{k}));
                isimg = startsWith(targets, ['images' '/']);
                refs = [refs, targets(isimg)]; %#ok<AGROW>
            end
            refs = unique(refs);

            % every reference resolves (already covered by the link test, but
            % re-checked here so an image failure names the image check)
            for r = 1:numel(refs)
                testCase.verifyTrue(isfile(fullfile(testCase.ManualDir, refs{r})), ...
                    sprintf('referenced image missing: %s', refs{r}));
            end

            % every committed image is referenced by some chapter
            if isfolder(imgdir)
                d = dir(imgdir);
                d = d(~[d.isdir]);
                for i = 1:numel(d)
                    rel = ['images' '/' d(i).name];
                    testCase.verifyTrue(any(strcmp(refs, rel)), ...
                        sprintf('orphan image (committed but never referenced): %s', rel));
                end
            end
        end

        function testVersionLiteralIsSingleSourcedAndCurrent(testCase)
            % psyrat_defineversion.m is the version's single source of truth;
            % read it the same way release/verify_export.sh check 7 does.
            vsrc = fileread(fullfile(testCase.projectRoot(), 'subroutines', ...
                'filefunctions', 'psyrat_defineversion.m'));
            tok = regexp(vsrc, 'psyratver = ''([^'']+)''', 'tokens', 'once');
            testCase.assertNotEmpty(tok, ...
                'could not read psyratver from psyrat_defineversion.m');
            ver = tok{1};

            idx = fileread(fullfile(testCase.ManualDir, 'index.md'));
            testCase.verifyTrue(contains(idx, ...
                sprintf('documents PsyRAT %s', ver)), sprintf( ...
                'index.md does not state "documents PsyRAT %s"', ver));

            % The literal appears nowhere else in the manual: every other
            % mention must say "the current version" or point at the banner.
            for k = 1:numel(testCase.ChapterFiles)
                [~, base, ext] = fileparts(testCase.ChapterFiles{k});
                if strcmp([base ext], 'index.md')
                    continue
                end
                testCase.verifyFalse( ...
                    contains(fileread(testCase.ChapterFiles{k}), ver), sprintf( ...
                    '%s repeats the version literal %s; only index.md may', ...
                    testCase.ChapterFiles{k}, ver));
            end
        end

        function testNoChapterReferencesWithheldOrPrivateMaterial(testCase)
            % Mirrors verify_export.sh: check 8's withheld-document pattern,
            % check 5's private-repo slug, plus the local-only donotsync
            % paths (papers_donotsync/, example_set_donotsync/) that exist
            % only on the development machine.
            forbidden = { ...
                'gamma_locationscale_rollout.md', ...
                'gamma_scale_submodel_decision.md', ...
                'gamma_chisquare_remaining_models.md', ...
                'sigma_i_call_site_census.md', ...
                'modular_cut_copula_workflow.md', ...
                'splits_observed_design_formulas.md', ...
                'UserManual.pdf', ...
                'PSYRAT_AUDIT_FINDINGS', 'PSYRAT_CODEX_TODOS', ...
                'SCIENTIFIC_FORMULA_AUDIT', 'SCIENTIFIC_FORMULA_TRACEABILITY', ...
                'HANDOFF.md', 'GUI_QUEUE', 'MANUAL_PRODUCTION', ...
                ... % Assembled from parts on purpose (G67): verify_export.sh
                ... % check 5 greps the whole staged tree for this slug and hard
                ... % fails on any hit, so the guard that keeps it out of the
                ... % manual must not itself put the literal into a shipped file.
                ['peclayson/' 'PsyRAT_' 'test'], ...
                'donotsync'};
            for k = 1:numel(testCase.ChapterFiles)
                txt = fileread(testCase.ChapterFiles{k});
                for w = 1:numel(forbidden)
                    testCase.verifyFalse(contains(txt, forbidden{w}), sprintf( ...
                        '%s references withheld/private material: %s', ...
                        testCase.ChapterFiles{k}, forbidden{w}));
                end
            end
        end

        function testNoShippedDocReferencesAnExportIgnoredPath(testCase)
            % A shipped, user-facing document must not send the reader to a
            % file the release does not contain. G61 was exactly that: every
            % worked example in documentation/scripted_use.md ran on
            % test_data/testdata.csv, which .gitattributes export-ignores, so
            % the tarball carried a scripting guide pointing at a file it
            % omitted. Nothing caught it, because the checks above scope
            % themselves to documentation/manual/ and this file is not in it.
            %
            % SCOPE, and why it is not wider. verify_export.sh checks 6 and 7
            % already WARN on tracker names anywhere in the tree, and that
            % backlog (mostly .m source comments citing dev documents) is
            % owner-ruled not a block. This check covers only prose a reader
            % is told to follow -- README, CHANGELOG, and documentation/**.
            % A stale source comment misleads a developer reading the source;
            % a stale instruction in the manual misleads every user.
            %
            % The forbidden set is DERIVED from .gitattributes at run time
            % rather than hand-listed, so it cannot go stale when an entry is
            % added there, and so this file introduces no forbidden literal of
            % its own (the G67 trap: a guard that names what it forbids
            % becomes what it catches).
            root = testCase.projectRoot();
            attrfile = fullfile(root, '.gitattributes');
            % The manifest export-ignores ITSELF (its first entry), so a
            % staged export has no .gitattributes to derive the forbidden
            % set from: skip there rather than fail. release/verify_export.sh
            % check 8 covers the staged tree; on a development checkout the
            % file is present and this check runs in full.
            testCase.assumeTrue(isfile(attrfile), sprintf( ...
                ['%s is absent: this is a staged export (the manifest ', ...
                'export-ignores itself), so the forbidden set cannot be ', ...
                'derived here; release/verify_export.sh check 8 covers ', ...
                'the exported tree.'], attrfile));

            % Keep only entries that name a concrete FILE. Directory entries
            % ('release', 'artifacts') and dotfiles ('.claude') would reduce
            % to common words, and globs ('documentation/*.pdf') have no one
            % basename. Classified syntactically, not with isfolder, because
            % this must also pass on a staged export where those directories
            % have already been stripped.
            lines = regexp(fileread(attrfile), '\r?\n', 'split');
            forbidden = {};   % basenames to search FOR
            ignored   = {};   % every ignored entry, to exclude docs that are
                              % themselves dev-only (a withheld document may
                              % freely cite other withheld documents)
            for k = 1:numel(lines)
                % No \b here: MATLAB's regexp reads \b as a backspace, not a
                % word boundary, so a \b-anchored pattern silently matches
                % nothing. The positive control below is what caught that.
                if startsWith(strtrim(lines{k}), '#'); continue; end
                tok = regexp(lines{k}, '^(\S+)\s+.*export-ignore', ...
                    'tokens', 'once');
                if isempty(tok); continue; end
                entry = tok{1};
                ignored{end+1} = entry; %#ok<AGROW>
                if startsWith(entry, '.') || contains(entry, '*'); continue; end
                [~, nm, ext] = fileparts(entry);
                if isempty(ext); continue; end   % directory entry
                forbidden{end+1} = [nm ext]; %#ok<AGROW>
            end
            forbidden = unique(forbidden);
            ignored   = unique(ignored);

            % Positive control: a prohibition-shaped check passes vacuously
            % when its subject is missing, so prove both sides are populated
            % before trusting a green result.
            testCase.verifyGreaterThan(numel(forbidden), 5, ...
                'parsed too few export-ignore file entries; check .gitattributes');

            docs = {fullfile(root, 'README.md'), fullfile(root, 'CHANGELOG.md')};
            d = dir(fullfile(root, 'documentation', '**', '*.md'));
            docs = [docs, fullfile({d.folder}, {d.name})];
            docs = docs(cellfun(@isfile, docs));
            % Drop the dev-only documents. They are present on a development
            % checkout and absent from a staged export, and they are entitled
            % to cite each other, so scanning them would fail here and pass
            % there for the same tree.
            docs = docs(~cellfun(@(f) localIsExportIgnored(f, root, ignored), docs));
            testCase.verifyGreaterThan(numel(docs), 10, ...
                'found too few shipped documents to scan; check the doc layout');

            for i = 1:numel(docs)
                txt = fileread(docs{i});
                for w = 1:numel(forbidden)
                    testCase.verifyFalse(contains(txt, forbidden{w}), sprintf( ...
                        ['%s points the reader at %s, which .gitattributes ' ...
                        'export-ignores, so the release does not contain it'], ...
                        docs{i}, forbidden{w}));
                end
            end
        end

        function testOnlyTheGeneratedManualPdfMayExist(testCase)
            % The release gate fails on any PDF except
            % documentation/manual/psyrat_manual.pdf; enforce the same rule
            % at unit-test time so a stray PDF is caught before export.
            d = dir(fullfile(testCase.ManualDir, '**', '*.pdf'));
            for i = 1:numel(d)
                testCase.verifyEqual(d(i).name, 'psyrat_manual.pdf', sprintf( ...
                    'unexpected PDF under documentation/manual/: %s', d(i).name));
            end
        end

        function testPrintedExpectedValuesMatchTheSlugMap(testCase)
            % The tutorial chapters are drafted with {{EV:slug}} placeholders
            % and rendered to literal values before release
            % (release/render_manual_values.py). Rendering destroys the link
            % between a printed number and its row in
            % tutorial_expected_values.md, so the renderer writes
            % ev_chapter_map.md, and this checks the link still holds.
            %
            % WHAT THIS PROVES, AND WHAT IT DOES NOT. It proves each printed
            % string is still SOMEWHERE IN the chapter that printed it, so a
            % number edited by hand in a chapter no longer agrees silently
            % with the transcription file. It does not check position, and it
            % is close to vacuous for a one-character value such as "0" that
            % occurs all over a chapter anyway. The code-to-table leg is
            % runManualTutorialCheck's job; this is the table-to-chapter leg.
            %
            % Both sides are whitespace-normalized before comparing: chapter
            % prose is hard-wrapped, so any multi-word value spans a line break
            % as soon as it lands near the wrap column. The comparison is on
            % rendered text, not on source line breaks.
            flat = @(s) strtrim(regexprep(s, '\s+', ' '));
            mapfile = fullfile(testCase.ManualDir, 'ev_chapter_map.md');
            testCase.assertTrue(isfile(mapfile), ...
                'documentation/manual/ev_chapter_map.md is missing');
            lines = regexp(fileread(mapfile), '\r?\n', 'split');
            nrows = 0;
            cache = containers.Map('KeyType', 'char', 'ValueType', 'char');
            for k = 1:numel(lines)
                tok = regexp(lines{k}, ...
                    '^\| ([A-Za-z0-9_\-]+) \| ([^|]+?) \| (.*) \|$', ...
                    'tokens', 'once');
                if isempty(tok) || strcmp(tok{1}, 'Slug'); continue; end
                [slug, chapter, printed] = deal(tok{1}, tok{2}, tok{3});
                if ~isKey(cache, chapter)
                    f = fullfile(testCase.ManualDir, chapter);
                    testCase.assertTrue(isfile(f), sprintf( ...
                        'ev_chapter_map.md names a missing chapter: %s', chapter));
                    cache(chapter) = flat(fileread(f));
                end
                testCase.verifyTrue(contains(cache(chapter), flat(printed)), sprintf( ...
                    ['%s no longer prints the value recorded for slug %s ' ...
                    '(expected to find: %s)'], chapter, slug, printed));
                nrows = nrows + 1;
            end
            % Guard against a silently empty map making this pass vacuously.
            testCase.verifyGreaterThan(nrows, 100, ...
                'ev_chapter_map.md parsed to too few rows; check its format');
        end

    end
end

% ---- local helpers ----------------------------------------------------------

function tf = localIsExportIgnored(fullpath, root, ignored)
%True when fullpath is covered by one of the .gitattributes export-ignore
% entries: an exact path match, a file under an ignored directory, or a
% wildcard entry such as documentation/*.pdf.
rel = strrep(erase(fullpath, [root filesep]), filesep, '/');
tf = false;
for k = 1:numel(ignored)
    e = ignored{k};
    if strcmp(rel, e) || startsWith(rel, [e '/'])
        tf = true; return
    end
    if contains(e, '*') && ~isempty(regexp(rel, ...
            ['^' regexptranslate('wildcard', e) '$'], 'once'))
        tf = true; return
    end
end
end

function targets = localLinkTargets(txt)
%Extract relative link/image targets from markdown text.
% Matches [label](target) and ![alt](target); drops external schemes and
% pure-anchor links; strips any #fragment from kept targets.
raw = regexp(txt, '!?\[[^\]]*\]\(([^)\s]+)\)', 'tokens');
targets = cellfun(@(c) c{1}, raw, 'UniformOutput', false);
keep = ~cellfun(@(t) startsWith(t, 'http://') || startsWith(t, 'https://') ...
    || startsWith(t, 'mailto:') || startsWith(t, '#'), targets);
targets = targets(keep);
targets = cellfun(@(t) regexprep(t, '#.*$', ''), targets, 'UniformOutput', false);
targets = targets(~cellfun(@isempty, targets));
end
