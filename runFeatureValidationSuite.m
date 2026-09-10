function [allPassed, results] = runFeatureValidationSuite(varargin)
%RUNFEATUREVALIDATIONSUITE Run end-to-end PsyRAT feature validation tests.
%
% This wrapper runs the integration tests that validate the workflow from
% variance component estimation in CmdStan through reliability summaries.
%
% Optional name-value pairs:
%   'Verbose'    : true/false (default: true)
%   'ReportJUnit': true/false (default: false)
%   'JUnitFile'  : JUnit XML output path (default: feature-validation.xml)
%
% Dependency path overrides can be set with environment variables:
%   PSYRAT_TEST_MPM_PATH
%   PSYRAT_TEST_MATLABSTAN_PATH
%   PSYRAT_TEST_CMDSTAN_PATH

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
parser.addParameter('JUnitFile', 'feature-validation.xml', @(x) ischar(x) || isstring(x));
parser.parse(varargin{:});
opts = parser.Results;

[allPassed, results] = runAllUnitTests( ...
    'Verbose', opts.Verbose, ...
    'ReportJUnit', opts.ReportJUnit, ...
    'JUnitFile', opts.JUnitFile, ...
    'Select', 'tests/integration');

% A CmdStan-less machine skips every integration test via assumptions, and the
% runner correctly does not fail intended skips, so without this guard the
% suite can return green having executed nothing. The release gate needs
% "green AND something ran".
nExecuted = numel(results) - nnz([results.Incomplete]);
fprintf('Feature suite executed (non-skipped) tests: %d\n', nExecuted);
if nExecuted == 0
    warning('PsyRAT:featureSuiteAllSkipped', ...
        ['No integration test actually executed (CmdStan/MatlabStan likely ' ...
        'unavailable); reporting failure rather than a vacuous pass.']);
    allPassed = false;
end
end
