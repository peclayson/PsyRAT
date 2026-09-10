classdef TestReadtableDelimiterOption < PsyRATTestBase
    %Audit 2026-09-05: psyrat_readtable documented a 'delimiter' option but
    %its parser looked for the key 'idcol', so the option was unreachable and
    %an unusual delimiter could be read only if MATLAB's own auto-detection
    %or the built-in cycling list (comma, space, tab, semicolon, pipe)
    %happened to cover it.
    %
    %RED against the pre-fix code: an '@'-delimited file is outside the
    %cycling list AND outside readtable's auto-detection, so without the
    %option the reader errors; with the option honored it reads two columns.
    %The delimiter is deliberately '@': the first version used '#', which
    %readtable on R2025b (Update 5) auto-detects (as it does nearly every
    %single character; '@', '$' and '+' were the exceptions when probed on
    %2026-09-07), so both tests passed with or without the fix. The choice is
    %MATLAB-version dependent; if a future release auto-detects '@', the
    %positive control below fails first and says so.

    methods (Test)

        function testDelimiterOptionIsHonored(testCase)
            tdir = tempname; mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
            f = fullfile(tdir, 'at_delimited.txt');
            fid = fopen(f, 'w');
            fprintf(fid, 'id@meas\ns1@1.1\ns2@2.2\n');
            fclose(fid);

            out = psyrat_readtable('file', f, 'delimiter', '@');
            testCase.verifyEqual(width(out), 2);
            testCase.verifyEqual(height(out), 2);
            testCase.verifyTrue(any(strcmpi(out.Properties.VariableNames, 'meas')));
        end

        function testWithoutTheOptionAnUnlistedDelimiterStillFails(testCase)
            %positive control for the premise: neither readtable's auto-detection
            %nor the cycling list splits on '@', so the option is what made the
            %read above succeed
            tdir = tempname; mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
            f = fullfile(tdir, 'at_delimited.txt');
            fid = fopen(f, 'w');
            fprintf(fid, 'id@meas\ns1@1.1\ns2@2.2\n');
            fclose(fid);
            testCase.verifyError(@() psyrat_readtable('file', f), 'loadingfile:delimiter');
        end

    end
end
