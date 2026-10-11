classdef TestStansummaryCsvParse < matlab.unittest.TestCase
    %TESTSTANSUMMARYCSVPARSE Named, full-precision parsing of the stansummary
    %CSV (work queue Q4; post-beta review PBA-01).
    %
    % The text parser inside psyrat_storeconv read the third-from-last token of
    % each printed row as n_eff and the last as R-hat, at the two significant
    % figures print(fit) requests. That token is N_Eff on CmdStan < 2.33 and
    % ESS_tail on >= 2.33 (the tested 2.38 layout ends
    % "ESS_bulk ESS_tail ESS_bulk/s R_hat"), so the stored n_eff silently
    % changed meaning with the CmdStan version and bulk ESS was lost.
    % psyrat_parse_stansummary_csv reads the CSV that
    % "stansummary --csv_filename" writes, maps columns by header name, stores
    % min(ESS_bulk, ESS_tail) as the gate's n_eff and keeps bulk ESS, tail ESS
    % and MCSE beside it. The modern fixture below is copied from a real
    % stansummary 2.39.0 run on tests/integration/cmdstan_artifacts/output-1.csv
    % (--sig_figs=8); the legacy fixture is the pre-2.33 column set.

    methods (Test)

        function testModernLayoutStoresMinOfBulkAndTail(testCase)
            f = localWrite(testCase, localModernCsv());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'lp__', 'b'}));
            testCase.verifyEqual(rows(:, 1), {'lp__'; 'b[1]'; 'b[2]'; 'b_sigma[1]'; 'b_sigma[2]'});
            testCase.verifySize(rows, [5 6]);
            % b[1]: bulk 2.7181084 and tail 7, so the gate's n_eff is the bulk value
            testCase.verifyEqual(rows{2, 2}, 2.7181084, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 3}, 2.1048832, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 4}, 2.7181084, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 5}, 7, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 6}, 0.0013412905, 'AbsTol', 1e-12);
            % lp__: bulk 2.7130949, tail 7 -> n_eff is again the bulk value; the
            % text parser stored the tail value 7 for both rows
            testCase.verifyEqual(rows{1, 2}, 2.7130949, 'AbsTol', 1e-9);
        end

        function testMinGoesEitherWay(testCase)
            % bulk 140 / tail 276 stores 140; bulk 300 / tail 120 stores 120
            f = localWrite(testCase, localModernCsvTwoRows());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'mu', 'sig_u'}));
            testCase.verifySize(rows, [2 6]);
            testCase.verifyEqual(rows{1, 2}, 140, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 2}, 120, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{1, 4}, 140, 'AbsTol', 1e-9);   % ess_bulk kept
            testCase.verifyEqual(rows{1, 5}, 276, 'AbsTol', 1e-9);   % ess_tail kept
        end

        function testLegacyLayoutRetainsNEffWithoutInventingBulkTail(testCase)
            f = localWrite(testCase, localLegacyCsv());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'lp__', 'mu', 'sig_u'}));
            testCase.verifyEqual(rows(:, 1), {'lp__'; 'mu'; 'sig_u'});
            testCase.verifySize(rows, [3 6]);
            testCase.verifyEqual(rows{2, 2}, 1834.5, 'AbsTol', 1e-9);   % n_eff = N_Eff
            testCase.verifyTrue(isnan(rows{2, 4}));   % bulk ESS was not reported
            testCase.verifyTrue(isnan(rows{2, 5}));   % tail ESS was not reported
            testCase.verifyEqual(rows{2, 3}, 1.0012345, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{2, 6}, 0.0123, 'AbsTol', 1e-12);
        end

        function testPrefixFilterDropsUnmonitoredRows(testCase)
            f = localWrite(testCase, localModernCsv());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'b_sigma'}));
            testCase.verifyEqual(rows(:, 1), {'b_sigma[1]'; 'b_sigma[2]'});
            none = psyrat_parse_stansummary_csv(f, localKeep({'zzz'}));
            testCase.verifyTrue(isempty(none) || size(none, 1) == 0);
        end

        function testQuotedNamesWithCommasSurvive(testCase)
            % matrix elements print as "cor_id[1,2]": the comma sits inside the
            % quotes and must not split the row
            f = localWrite(testCase, localModernCsvMatrix());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'cor_id', 'sd_id'}));
            testCase.verifyEqual(rows(:, 1), {'sd_id[1]'; 'cor_id[1,1]'; 'cor_id[1,2]'});
            testCase.verifyEqual(rows{3, 3}, 1.0007, 'AbsTol', 1e-9);
            testCase.verifyEqual(rows{3, 2}, 812, 'AbsTol', 1e-9);   % min(1234, 812)
            % the constant diagonal keeps nan in every diagnostic
            testCase.verifyTrue(isnan(rows{2, 2}) && isnan(rows{2, 3}));
        end

        function testNanAndQuotesAreHandled(testCase)
            % sampler diagnostics rows carry literal nan; names are double-quoted
            f = localWrite(testCase, localModernCsv());
            rows = psyrat_parse_stansummary_csv(f, localKeep({'stepsize__', 'lp__'}));
            % rows come back in FILE order, not in the order of the prefixes
            testCase.verifyEqual(rows(:, 1), {'lp__'; 'stepsize__'});
            testCase.verifyTrue(isnan(rows{2, 2}) && isnan(rows{2, 3}) ...
                && isnan(rows{2, 4}) && isnan(rows{2, 5}) && isnan(rows{2, 6}), ...
                'a constant row keeps nan in every diagnostic, including n_eff');
            testCase.verifyFalse(isnan(rows{1, 2}));
            testCase.verifyFalse(any(contains(rows(:, 1), '"')));
        end

        function testMissingFileErrors(testCase)
            testCase.verifyError(@() psyrat_parse_stansummary_csv( ...
                fullfile(tempdir, 'psyrat_no_such_summary_file.csv'), localKeep({'lp__'})), ...
                'PsyRAT:stansummaryCsvMissing');
        end

    end
end

function keep = localKeep(prefixes)
% the same prefix rule psyrat_storeconv applies through psyrat_hasprefix
keep = @(name) any(cellfun(@(p) strncmp(name, p, length(p)), prefixes));
end

function f = localWrite(testCase, txt)
f = [tempname '.csv'];
fid = fopen(f, 'w');
testCase.assertGreaterThan(fid, 0, 'could not open a temp file');
fprintf(fid, '%s', txt);
fclose(fid);
testCase.addTeardown(@() delete(f));
end

function txt = localModernCsv()
txt = sprintf(['name,Mean,MCSE,StdDev,MAD,5%%,50%%,95%%,ESS_bulk,ESS_tail,ESS_bulk/s,R_hat\n' ...
    '"lp__",-221749.27,64000.625,142597.75,86525.767,-502583.86,-186331.28,-120808.68,2.7130949,7,387.58499,2.1048832\n' ...
    '"accept_stat__",1,nan,0,0,1,1,1,nan,nan,nan,nan\n' ...
    '"stepsize__",2.5009742e-05,nan,0,0,2.5009742e-05,2.5009742e-05,2.5009742e-05,nan,nan,nan,nan\n' ...
    '"divergent__",0,nan,0,0,0,0,0,nan,nan,nan,nan\n' ...
    '"b[1]",0.001404698,0.0013412905,0.0026135436,0.0033194634,-0.0010051405,0.0015478678,0.0060018418,2.7181084,7,388.3012,2.1048832\n' ...
    '"b[2]",1.863838,4.1914841e-05,8.5716788e-05,5.9304e-05,1.8637212,1.8638242,1.8639864,2.9346157,7,419.23082,2.1048832\n' ...
    '"b_sigma[1]",0.77441914,0.0044341892,0.0086672387,0.010531664,0.7591659,0.77411359,0.78249556,2.7130949,7,387.58499,2.1048832\n' ...
    '"b_sigma[2]",-0.37532421,0.00030108951,0.00056129196,0.00088367408,-0.37621836,-0.37544539,-0.3747614,2.7339434,7,390.56334,2.1048832\n' ...
    '# For each parameter, ESS_bulk and ESS_tail measure the effective sample size for the entire sample (bulk) and for the .\n' ...
    '# and R_hat measures the potential scale reduction on split chains. At convergence R_hat will be very close to 1.00.\n']);
end

function txt = localModernCsvTwoRows()
txt = sprintf(['name,Mean,MCSE,StdDev,MAD,5%%,50%%,95%%,ESS_bulk,ESS_tail,ESS_bulk/s,R_hat\n' ...
    '"mu",5.7,0.054,0.64,0.67,4.6,5.7,6.7,140,276,8.1,1.0023\n' ...
    '"sig_u",2.1,0.02,0.3,0.31,1.6,2.1,2.7,300,120,17,1.0011\n']);
end

function txt = localModernCsvMatrix()
txt = sprintf(['name,Mean,MCSE,StdDev,MAD,5%%,50%%,95%%,ESS_bulk,ESS_tail,ESS_bulk/s,R_hat\n' ...
    '"sd_id[1]",0.62,0.004,0.08,0.08,0.5,0.62,0.75,900,650,40,1.0021\n' ...
    '"cor_id[1,1]",1,nan,0,0,1,1,1,nan,nan,nan,nan\n' ...
    '"cor_id[1,2]",0.31,0.006,0.12,0.12,0.1,0.31,0.5,1234,812,55,1.0007\n']);
end

function txt = localLegacyCsv()
txt = sprintf(['name,Mean,MCSE,StdDev,5%%,50%%,95%%,N_Eff,N_Eff/s,R_hat\n' ...
    '"lp__",-1200.5,0.5,12.1,-1220,-1200,-1180,600.25,30.1,1.0034567\n' ...
    '"mu",0.5,0.0123,0.25,0.1,0.5,0.9,1834.5,90.2,1.0012345\n' ...
    '"sig_u",2.1,0.02,0.3,1.6,2.1,2.7,900,45,1.0009\n']);
end
