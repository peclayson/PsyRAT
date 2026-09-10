classdef TestCheckconvFailsClosed < PsyRATTestBase
    %Audit 2026-09-05: psyrat_checkconv certified convergence in two edge
    %states. A NaN R-hat or n_eff (what stansummary prints for a parameter
    %with zero variance across draws, parsed by str2double as NaN) was never
    %flagged, because NaN >= 1.1 and NaN < floor are both false; and a
    %convergence table with no parameter rows (no summary row matched the
    %design's prefix list) passed because find([] ...) is empty. Both now
    %read as not converged. The positive controls pin that finite values are
    %judged exactly as before, so ordinary runs are unchanged.
    %
    %Refined 2026-09-07 (owner ruling, at the tutorial refit): a row with NaN in
    %BOTH columns is a parameter that is constant by construction (correlation-
    %matrix diagonals fixed at 1; disabled covariance and residual-correlation
    %placeholders fixed at 0), which stansummary reports exactly that way, and is
    %excluded from the gate; the first form of this fix flagged those rows and
    %chapter 15's subject-level difference fit could never converge. A NaN in
    %only one column, an Inf R-hat, and an all-constant table still fail.
    %
    %RED against the pre-fix code: the four fail-closed tests and the Inf and
    %all-constant arms. The finite-values control and the constant-rows test
    %pass on both builds (the old gate ignored NaN altogether), so they are
    %compatibility pins, not differentials.

    methods (Test)

        function testNanRhatIsNotConverged(testCase)
            rel = localRel({{ ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', 120, NaN; ...
                'sig', 140, 1.01}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'A NaN R-hat must fail the convergence gate.');
        end

        function testNanNeffIsNotConverged(testCase)
            rel = localRel({{ ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', NaN, 1.00; ...
                'sig', 140, 1.01}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'A NaN n_eff must fail the convergence gate.');
        end

        function testHeaderOnlyTableIsNotConverged(testCase)
            %nested-cell layout (one fit) and the flat layout both
            rel = localRel({{'name', 'n_eff', 'Rhat'}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'A table with no monitored parameter must not read as converged (nested).');
            rel.out.conv.data = {'name', 'n_eff', 'Rhat'};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'A table with no monitored parameter must not read as converged (flat).');
        end

        function testNoTablesIsNotConvergedAndDoesNotError(testCase)
            rel = localRel(cell(1, 0));
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'A result with no convergence tables must fail closed, not error.');
        end

        function testFiniteValuesJudgedAsBefore(testCase)
            %positive controls: a clean table converges, a bad R-hat and a low
            %n_eff still fail, in both layouts
            rel = localRel({{ ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', 120, 1.00; ...
                'sig', 140, 1.01}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1);

            rel.out.conv.data = {{'name', 'n_eff', 'Rhat'; 'mu', 120, 1.10}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, 'R-hat at the 1.1 threshold fails');

            rel.out.conv.data = {{'name', 'n_eff', 'Rhat'; 'mu', 39, 1.00}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, 'n_eff below 10 x chains fails');

            rel.out.conv.data = {'name', 'n_eff', 'Rhat'; 'mu', 40, 1.09};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1, 'flat layout at the boundary passes');
        end

        function testMultiFitNanInSecondTableFails(testCase)
            %the nested layout loops over fits; a NaN in any table fails the run
            rel = localRel({ ...
                {'name', 'n_eff', 'Rhat'; 'mu', 120, 1.00}, ...
                {'name', 'n_eff', 'Rhat'; 'mu', 120, NaN}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0);
        end

        function testConstantRowsAreExcludedFromTheGate(testCase)
            %rows with NaN in BOTH columns are constants by construction (what
            %stansummary prints for cor_id[1,1] = 1 or a disabled rescor = 0)
            %and must neither fail a converged fit nor mask a failing one, in
            %both layouts
            rel = localRel({{ ...
                'name', 'n_eff', 'Rhat'; ...
                'sig_id', 120, 1.00; ...
                'cor_id[1,1]', NaN, NaN; ...
                'rescor[1]', NaN, NaN}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1, ...
                'Constant rows must not fail a converged fit (nested).');

            rel.out.conv.data = {'name', 'n_eff', 'Rhat'; ...
                'sig_id', 120, 1.00; 'cor_id[1,1]', NaN, NaN};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1, ...
                'Constant rows must not fail a converged fit (flat).');

            rel.out.conv.data = {{'name', 'n_eff', 'Rhat'; ...
                'sig_id', 120, 1.20; 'cor_id[1,1]', NaN, NaN}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'Constant rows must not mask a bad R-hat on a varying row.');

            rel.out.conv.data = {{'name', 'n_eff', 'Rhat'; ...
                'sig_id', 12, 1.00; 'cor_id[1,1]', NaN, NaN}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'Constant rows must not mask a low n_eff on a varying row.');
        end

        function testAllConstantAndInfRhatStillFail(testCase)
            %a table whose every row is a constant showed nothing converging
            %(same reading as the header-only table), and chains stuck at
            %different constants print an Inf R-hat, which still fails
            rel = localRel({{'name', 'n_eff', 'Rhat'; 'cor_id[1,1]', NaN, NaN}});
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'An all-constant table must not read as converged (nested).');

            rel.out.conv.data = {'name', 'n_eff', 'Rhat'; 'cor_id[1,1]', NaN, NaN};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'An all-constant table must not read as converged (flat).');

            rel.out.conv.data = {{'name', 'n_eff', 'Rhat'; 'mu', 120, Inf}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0, ...
                'An Inf R-hat must fail the convergence gate.');
        end

    end
end

function rel = localRel(convdata)
rel = struct();
rel.nchains = 4;
rel.out.conv.data = convdata;
end
