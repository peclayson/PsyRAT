classdef TestUpdaterUsesReleasesApi < matlab.unittest.TestCase
    %Source pins for psyrat_updatepsyrat (audit finding G87, owner ruling
    %2026-09-09): the toolbox updater resolves the latest stable release
    %through the GitHub releases API (psyrat_checkgithubreleases) and
    %downloads the PsyRAT.zip asset that release/PUBLISHING.md attaches to
    %every release. Before the ruling it scraped the releases web page for a
    %heading element, parsed a version out of a title the runbook never
    %writes, and appended a guessed asset name to a download URL no release
    %carried; the modal update prompt had already closed every PsyRAT window
    %when that failed.
    %
    %The function cannot run headlessly (uiputfile, websave, rmpath and
    %savepath), so these are source pins. Every forbid pin is paired with a
    %positive pin on the same file, because a prohibition alone is satisfied
    %by an empty file.

    methods (Test)
        function testUpdaterResolvesTheReleaseThroughTheApi(testCase)
            code = localCode();
            testCase.verifyEqual(numel(strfind(code, 'psyrat_checkgithubreleases(')), 1, ...
                'the updater must resolve the release through psyrat_checkgithubreleases');
            testCase.verifyGreaterThanOrEqual(numel(strfind(code, 'latestReleaseAssetURL')), 2, ...
                ['the updater must read the asset URL the release check exposes, ' ...
                'once at the gate and once for the download']);
            testCase.verifyEqual(numel(strfind(code, 'latestReleaseVersion')), 1, ...
                'the install folder name must come from the API version');
        end

        function testUpdaterNoLongerScrapesTheReleasesPage(testCase)
            src = localSource();
            %positive control for the forbid pins below: an empty file fails here
            testCase.verifyEqual(numel(strfind(src, 'websave(')), 1, ...
                'the download itself must still be there');
            testCase.verifyEqual(numel(strfind(src, 'webread(')), 0, ...
                'the releases web page is not fetched');
            testCase.verifyEqual(numel(strfind(src, 'release-title')), 0, ...
                'no dependence on the page markup');
            testCase.verifyEqual(numel(strfind(src, '''/PsyRAT.zip''')), 0, ...
                'the download URL is the asset''s own URL, never an appended guess');
            testCase.verifyEqual(numel(strfind(src, '''Version')), 0, ...
                'no title parse: no ''Version'' string literal to match a heading against');
        end

        function testUpdaterRefusesWithARelaunchWhenNothingIsDownloadable(testCase)
            src = localSource();
            pat = 'psyrat_start;\s+errordlg\(\{[^}]*\}, ''Update unavailable'', ''modal''\);';
            testCase.verifyEqual(numel(regexp(src, pat, 'start')), 1, ...
                ['when the API reports no stable release or no PsyRAT.zip asset, ' ...
                'the updater must relaunch psyrat_start FIRST and raise one modal ' ...
                'Update unavailable dialog LAST (B20 creation order)']);
        end
    end
end

function src = localSource()
%absolute path from THIS file's location, never from pwd or the MATLAB path
here = fileparts(mfilename('fullpath'));
src = fileread(fullfile(fileparts(here), 'subroutines', 'installation', ...
    'psyrat_updatepsyrat.m'));
end

function code = localCode()
%the same file with whole-line comments removed
lines = strsplit(localSource(), newline);
keep = ~cellfun(@(l) startsWith(strtrim(l), '%'), lines);
code = strjoin(lines(keep), newline);
end
