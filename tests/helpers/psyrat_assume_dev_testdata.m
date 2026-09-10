function psyrat_assume_dev_testdata(testCase)
%PSYRAT_ASSUME_DEV_TESTDATA Skip (not fail) without the dev-only source fixture.
%
% test_data/testdata.csv is a development-only fixture: it is export-ignored
% and does not ship in public releases (owner ruling 2026-08-29, "don't ship
% any of our test data"). PsyRATCmdStanIntegrationFactory derives every one
% of its scenario tables from that file, so any test that goes through the
% factory must ASSUME the file is present, mirroring the established pattern
% for the real-EP fixture in TestImportWorkflow: on a development checkout
% the tests run in full; on a public release tree they skip with this
% message, and runAllUnitTests counts assumption-filtered skips as expected.
%
% Call this at the top of any Test method (or testCase-carrying local
% function) that reaches PsyRATCmdStanIntegrationFactory.makeScenarioTable.

fixture = fullfile(testCase.projectRoot(), 'test_data', 'testdata.csv');
testCase.assumeTrue(isfile(fixture), sprintf( ...
    ['Development-only source fixture not present: %s\n', ...
    'It is not distributed in public releases, so the scenario tables ', ...
    'derived from it cannot be built here. On a development checkout this ', ...
    'file is tracked; if it is missing there, restore it from git.'], fixture));
end
