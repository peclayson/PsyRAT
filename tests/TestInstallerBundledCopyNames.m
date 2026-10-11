classdef TestInstallerBundledCopyNames < PsyRATTestBase
    %Sweep follow-up F4 (ruled 2026-09-16): the installer's two bundled-copy
    %helpers built the destination folder from the NAME output of fileparts,
    %which splits a dotted version at its last dot, so the copied dependents
    %lost their final version component (MatlabProcessManager-0.5.1 became
    %MatlabProcessManager-0.5, corroborated live in ~/Documents/MATLAB/
    %PsyRATDependents; MatlabStan-2.15.1.0 would become MatlabStan-2.15.1; a
    %future two-component pin such as 0.6 would become MatlabProcessManager-0).
    %Every lookup that decides a verdict is by content (the one name-based
    %lookup, psyrat_findmatlabstandir, matches a 10-character prefix that both
    %spellings share), so the defect was a wrong label, not a wrong verdict,
    %and the fix keeps the extension: [n e].
    %
    %The helpers are local functions of psyrat_installdependents.m and their
    %copy arm only runs on the legacy-upgrade path (a saved path that predates
    %the vendored MatlabProcessManager), so they cannot be executed here; this
    %class pins the SOURCE of both helpers and demonstrates the naming rule on
    %the pinned version strings.

    methods (Test)

        function testBothCopyHelpersKeepTheDottedVersion(testCase)
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'installation', 'psyrat_installdependents.m'));
            helpers = {'psyrat_copybundledmatlabstan', 'psyrat_copybundledmpm'};
            for k = 1:numel(helpers)
                %anchor on the DEFINITION line (the helpers are also called by
                %name earlier in the file, and a call-site segment would hold
                %none of the strings below and pass the forbid check vacuously)
                idef = regexp(src, ['\nfunction [^\n]*' helpers{k} '\(wrkdir\)'], 'once');
                testCase.assertNotEmpty(idef, sprintf( ...
                    '%s must still be defined under this name', helpers{k}));
                seg = src(idef + 1:end);
                nextf = strfind(seg, [newline 'function ']);
                if ~isempty(nextf)
                    seg = seg(1:nextf(1));
                end
                testCase.assertEqual(numel(strfind(seg, 'copyfile(bundledsrc,')), 1, sprintf( ...
                    '%s segment must be the helper body (positive control)', helpers{k}));
                testCase.verifyEqual(numel(strfind(seg, '[~,foldername] = fileparts(bundledsrc);')), 0, ...
                    sprintf(['%s must not drop the last version component: the NAME ' ...
                    'output of fileparts splits a dotted folder name at its last dot'], helpers{k}));
                testCase.verifyEqual(numel(strfind(seg, '[~,n,e] = fileparts(bundledsrc);')), 1, ...
                    sprintf('%s must take the name AND extension of the bundled folder', helpers{k}));
                testCase.verifyEqual(numel(strfind(seg, 'foldername = [n e];')), 1, ...
                    sprintf('%s must rebuild the destination name from both parts', helpers{k}));
            end
        end

        function testTheNamingRuleReconstructsThePinnedFolders(testCase)
            %The demonstration behind the pin: on the version strings the
            %installer actually pins, the old expression truncates and the new
            %one reconstructs the bundled folder name exactly.
            depvers = psyrat_dependentsversions;
            names = {['MatlabStan-' depvers.matlabstan], ...
                ['MatlabProcessManager-' depvers.matlabprocessmanager]};
            for k = 1:numel(names)
                bundledsrc = fullfile(tempdir, 'bundled_dependents', names{k});
                [~, old] = fileparts(bundledsrc);
                [~, n, e] = fileparts(bundledsrc);
                testCase.verifyEqual([n e], names{k}, ...
                    'name plus extension must give back the bundled folder name');
                testCase.verifyTrue(contains(names{k}, '.'), ...
                    'the pinned versions are dotted, which is what makes the old expression truncate');
                testCase.verifyNotEqual(old, names{k}, ...
                    'the old NAME-only expression must be shown to truncate this pin');
            end
            %the bundled folders themselves carry the pinned names
            root = testCase.projectRoot();
            for k = 1:numel(names)
                testCase.verifyEqual(exist(fullfile(root, 'bundled_dependents', names{k}), 'dir'), 7, ...
                    sprintf('bundled_dependents/%s must exist', names{k}));
            end
        end

        function testTheRmolddepsCommentNoLongerClaimsIdenticalNames(testCase)
            %The comment above rmolddeps' targets said the installed copies
            %"carry identical folder names"; with the truncation that was false
            %for MatlabProcessManager. After the fix the copies do carry the
            %bundled names, and the comment says so in those words.
            src = fileread(fullfile(testCase.projectRoot(), ...
                'subroutines', 'installation', 'psyrat_installdependents.m'));
            testCase.verifyEqual(numel(strfind(src, 'also carry the bundled folder names')), 1, ...
                'the rmolddeps comment must describe the copied names as the bundled names');
        end

    end
end
