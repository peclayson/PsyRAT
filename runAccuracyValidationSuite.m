function [allPassed, results] = runAccuracyValidationSuite(varargin)
%RUNACCURACYVALIDATIONSUITE Run deterministic PsyRAT calculation accuracy tests.
%
%   [allPassed, results] = runAccuracyValidationSuite()
%   [allPassed, results] = runAccuracyValidationSuite('Verbose', true, ...
%       'ReportJUnit', true, 'JUnitFile', 'accuracy-validation.xml')
%
% This suite is intentionally focused on deterministic calculation correctness
% tests -- no CmdStan, no sampling noise.
%
% WHAT EACH MEMBER OF THIS LANE CAN AND CANNOT PROVE. The distinction matters,
% because a green run here is easy to over-read:
%
%   TestCalculationAccuracyOracle    REGRESSION ONLY. PsyRATAccuracyOracle
%       duplicates the production formulas, so this proves the formulas are
%       UNCHANGED, not that they are correct -- both sides would have to be wrong
%       in the same way for it to notice.
%
%   TestIndependentOracleAgreement   EXTERNAL CORROBORATION. PsyRATIndependentOracle
%       states the published formulas of Rocha et al. (2026) in the paper's own
%       symbols. Its expressions have been checked term by term against Tables 2,
%       3 and 6 as re-extracted from the source PDF in documentation/, and they
%       agree. That verified agreement with the printed tables -- not a claim
%       about who typed them in -- is what makes agreement here evidence the
%       formulas are RIGHT, not merely unchanged. The class header records that
%       its Table 6 block has weaker provenance than Tables 2 and 3.
%
%   TestExternalBenchmark            EXTERNAL NUMERIC ANCHOR. Validates against
%       trial-count tables produced by a separate analysis (Heindorf et al.).
%       Covers one-facet, one-facet difference, and difference-of-differences only.
%
% Do not cite a green run of this suite as evidence of accuracy for a design
% family that none of the three covers.

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
parser.addParameter('JUnitFile', 'accuracy-validation.xml', @(x) ischar(x) || isstring(x));
parser.parse(varargin{:});
opts = parser.Results;

[allPassed, results] = runAllUnitTests( ...
    'Verbose', opts.Verbose, ...
    'ReportJUnit', opts.ReportJUnit, ...
    'JUnitFile', opts.JUnitFile, ...
    'Select', {'TestCalculationAccuracyOracle', ...
               'TestIndependentOracleAgreement', ...
               'TestExternalBenchmark'});
end
