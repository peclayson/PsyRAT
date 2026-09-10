function [allPassed, results] = runAllUnitTests(varargin)
%RUNALLUNITTESTS Run PsyRAT unit tests under ./tests (excluding integration).
%
%   [allPassed, results] = runAllUnitTests()
%   [allPassed, results] = runAllUnitTests('Verbose', true, ...
%       'ReportJUnit', true, 'JUnitFile', 'test-results.xml', ...
%       'Select', 'TestCalculationFunctions')
%   [allPassed, results] = runAllUnitTests('Select', 'TestCalculationAccuracyOracle')
%
% Integration tests under ./tests/integration are excluded by default.
% Use runFeatureValidationSuite to execute integration/feature validation.
% Use runAccuracyValidationSuite to execute deterministic oracle-backed
% calculation accuracy checks.

% Copyright (C) 2016-2025 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see
%     <http://www.gnu.org/licenses/>.

parser = inputParser;
parser.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('ReportJUnit', false, @(x) islogical(x) || isnumeric(x));
parser.addParameter('JUnitFile', 'test-results.xml', @(x) ischar(x) || isstring(x));
%'Select' accepts a single name/folder, or a string array / cellstr of several, so a
%caller can assemble a lane from more than one test class (see runAccuracyValidationSuite).
parser.addParameter('Select', '', @(x) ischar(x) || isstring(x) || iscellstr(x));
parser.parse(varargin{:});
opts = parser.Results;

projectRoot = fileparts(mfilename('fullpath'));
testsRoot = fullfile(projectRoot, 'tests');

if ~isfolder(testsRoot)
    error('runAllUnitTests:MissingTestsFolder', ...
        ['Unable to find the tests folder at: %s. How to fix: run this ',...
        'function from the PsyRAT project root, or restore the ./tests ',...
        'directory.'], testsRoot);
end

% Ensure test helper classes/functions are available during discovery.
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'tests', 'helpers')));
addpath(genpath(fullfile(projectRoot, 'subroutines')));

suite = localBuildSuite(testsRoot, opts.Select);
if isempty(suite)
    error('runAllUnitTests:NoTestsSelected', ...
        ['No tests matched the requested selection: %s. How to fix: use ',...
        'a class name substring (for example ''TestCalculationFunctions'') ',...
        'or leave ''Select'' empty to run all non-integration tests.'], ...
        char(string(opts.Select)));
end

if opts.Verbose
    outputDetail = matlab.unittest.Verbosity.Detailed;
else
    outputDetail = matlab.unittest.Verbosity.Concise;
end

runner = matlab.unittest.TestRunner.withTextOutput('OutputDetail', outputDetail);
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin());

if opts.ReportJUnit
    junitFile = char(opts.JUnitFile);
    if ~isfile(junitFile)
        [junitDir, ~, ~] = fileparts(junitFile);
        if ~isempty(junitDir) && ~isfolder(junitDir)
            mkdir(junitDir);
        end
    end
    xmlPlugin = matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(junitFile);
    runner.addPlugin(xmlPlugin);
    fprintf('JUnit report enabled: %s\n', junitFile);
end

fprintf('Running %d test(s) from %s\n', numel(suite), testsRoot);
results = runner.run(suite);

nPassed = nnz([results.Passed]);
nFailed = nnz([results.Failed]);
nIncomplete = nnz([results.Incomplete]);

% Distinguish *intended* skips from *genuine* incompletes. MATLAB reports both
% as Incomplete, but they mean different things:
%   - Assumption-filtered skips: a test called assumeTrue/assumeFalse with a
%     false condition and was skipped cleanly. The native-engine tests do this
%     to skip when an optional MATLAB toolbox is absent (Statistics_Toolbox for
%     fitlme/hmcSampler; Deep Learning Toolbox for dlarray/AD). On a base-MATLAB
%     CI runner these legitimately skip and must NOT fail the run.
%   - Everything else (for example an error thrown during fixture setup) is an
%     Incomplete we do want to surface as a failure.
% This mirrors the workflow's seeded-smoke lane, which already gates on Failed
% rather than on Incomplete.
nAssumptionFiltered = localCountAssumptionFiltered(results);
nUnexpectedIncomplete = nIncomplete - nAssumptionFiltered;

% The run is green when there are no failures and no *unexpected* incompletes.
% This is the value the CI gate keys off of (the first output, 'ok').
allPassed = (nFailed == 0) && (nUnexpectedIncomplete == 0);

fprintf('\nTest Summary\n');
fprintf('Passed: %d\n', nPassed);
fprintf('Failed: %d\n', nFailed);
fprintf('Incomplete: %d (assumption-filtered skips: %d, unexpected: %d)\n', ...
    nIncomplete, nAssumptionFiltered, nUnexpectedIncomplete);
fprintf('Run OK (no failures, no unexpected incompletes): %d\n', allPassed);
fprintf('Coverage report: not generated (no external coverage toolbox required).\n');
end

function n = localCountAssumptionFiltered(results)
%LOCALCOUNTASSUMPTIONFILTERED Count tests skipped via a failed assumption.
%
% A test that fails an assumption (assumeTrue/assumeFalse) is reported as
% Incomplete and records a diagnostic whose Event is 'AssumptionFailed'. This
% counts those so the caller can treat them as expected skips rather than
% failures. Requires the DiagnosticsRecordingPlugin (added above) to populate
% Details.DiagnosticRecord; an Incomplete test without that marker (for example
% one that errored in fixture setup) is deliberately NOT counted, so it still
% registers as an unexpected incomplete and fails the run.
n = 0;
for k = 1:numel(results)
    if ~results(k).Incomplete
        continue;
    end
    details = results(k).Details;
    if ~isfield(details, 'DiagnosticRecord')
        continue;
    end
    recs = details.DiagnosticRecord;
    for r = 1:numel(recs)
        if strcmp(recs(r).Event, 'AssumptionFailed')
            n = n + 1;
            break;
        end
    end
end
end

function suite = localBuildSuite(testsRoot, selectArg)
selectArg = string(selectArg);
selectArg = selectArg(strlength(selectArg) > 0);   %drop empties so {} and '' behave alike
allSuite = matlab.unittest.TestSuite.fromFolder(testsRoot, 'IncludingSubfolders', true);

if isempty(selectArg)
    integrationToken = [filesep 'tests' filesep 'integration'];
    isIntegration = contains({allSuite.BaseFolder}, integrationToken, ...
        'IgnoreCase', true);
    suite = allSuite(~isIntegration);
    return;
end

%A single selector naming a folder keeps its original meaning: run that folder.
if isscalar(selectArg) && isfolder(char(selectArg))
    suite = matlab.unittest.TestSuite.fromFolder(char(selectArg), 'IncludingSubfolders', true);
    return;
end

%Otherwise match names/folders against every selector and take the union, so a lane can
%be assembled from several test classes.
mask = false(1, numel(allSuite));
for k = 1:numel(selectArg)
    selectValue = char(selectArg(k));
    mask = mask | ...
        contains({allSuite.Name}, selectValue, 'IgnoreCase', true) | ...
        contains({allSuite.BaseFolder}, selectValue, 'IgnoreCase', true);
end
suite = allSuite(mask);
end
