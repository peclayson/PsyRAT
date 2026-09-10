classdef TestInterpretationSafety < PsyRATTestBase
    %Unit coverage for the interpretation-safety fixes (council review items
    %1-3): the viewer convergence/divergence notice (psyrat_convergence_notice),
    %the divergence totaller (psyrat_total_divergences), and the headless
    %gating of psyrat_relsummary's threshold-failure dialog.
    %
    %These exercise the new pure helpers headlessly. The end-to-end surfacing
    %(modal in the live GUI; divergent__ extraction from a real CmdStan fit) is
    %verified separately with a live run and a GUI click-through.

    methods (Test)

        % ---- Item 1: psyrat_convergence_notice ---------------------------

        function testNoticeEmptyForConvergedNoDivergences(testCase)
            rel = local_rel(1, 0);
            testCase.verifyEmpty(psyrat_convergence_notice(rel));
        end

        function testNoticeFlagsNonConvergence(testCase)
            rel = local_rel(0, 0);
            msg = psyrat_convergence_notice(rel);
            testCase.verifyNotEmpty(msg);
            testCase.verifyTrue(contains(lower(msg), 'did not converge'));
        end

        function testNoticeFlagsDivergences(testCase)
            rel = local_rel(1, 7);
            msg = psyrat_convergence_notice(rel);
            testCase.verifyNotEmpty(msg);
            testCase.verifyTrue(contains(lower(msg), 'divergent'));
            testCase.verifyTrue(contains(msg, '7'));
        end

        function testNoticeFlagsBothWhenBothApply(testCase)
            rel = local_rel(0, 3);
            msg = psyrat_convergence_notice(rel);
            testCase.verifyTrue(contains(lower(msg), 'did not converge'));
            testCase.verifyTrue(contains(lower(msg), 'divergent'));
        end

        function testNoticeNaNDivergencesSilent(testCase)
            %NaN means "not assessed" (native engines / unavailable): no notice
            %about divergences, and a converged run stays silent overall.
            rel = local_rel(1, NaN);
            testCase.verifyEmpty(psyrat_convergence_notice(rel));
        end

        function testNoticeLegacyStructNoError(testCase)
            %A results file predating the diagnostic fields must produce '' and
            %never error.
            testCase.verifyEmpty(psyrat_convergence_notice(struct()));
            rel = struct('out', struct('foo', 1));   % out present, no conv
            testCase.verifyEmpty(psyrat_convergence_notice(rel));
        end

        % ---- Item 2: psyrat_total_divergences ----------------------------

        function testTotalSingleTableHeaderCount(testCase)
            %Single convergence table (the a~=1 shape) with the count stamped in
            %the header cell (column 4).
            tbl = local_conv_table(5);
            testCase.verifyEqual(psyrat_total_divergences(tbl), 5);
        end

        function testTotalCellOfTablesSums(testCase)
            %Cell-of-tables (the a==1 append shape): sum across fits.
            data = {local_conv_table(2), local_conv_table(4), local_conv_table(0)};
            testCase.verifyEqual(psyrat_total_divergences(data), 6);
        end

        function testTotalThreeColumnTableIsNaN(testCase)
            %Legacy/native tables without a 4th column -> not assessed.
            tbl = {'name','n_eff','r_hat'; 'mu',120,1.00};
            testCase.verifyTrue(isnan(psyrat_total_divergences(tbl)));
        end

        function testTotalAllNaNIsNaN(testCase)
            data = {local_conv_table(NaN), local_conv_table(NaN)};
            testCase.verifyTrue(isnan(psyrat_total_divergences(data)));
        end

        function testTotalMixedRealAndNaNSumsReal(testCase)
            data = {local_conv_table(3), local_conv_table(NaN)};
            testCase.verifyEqual(psyrat_total_divergences(data), 3);
        end

        function testTotalNonCellIsNaN(testCase)
            testCase.verifyTrue(isnan(psyrat_total_divergences([])));
            testCase.verifyTrue(isnan(psyrat_total_divergences(42)));
        end

        % ---- Item 3: psyrat_relsummary headless gating -------------------

        function testRelSummaryAcceptsNonInteractiveFlag(testCase)
            %The new interactive-flag parsing must run cleanly and leave the
            %numeric result unchanged on good data (the flag only gates the
            %threshold-failure dialog, not the computation).
            base = PsyRATTestDataFactory.makeRawSingRelData();

            args = {'analysis','sing','depcutoff',0.70, ...
                'meascutoff',2,'depcentmeas',1,'CI',0.95};

            [outDefault, errDefault] = psyrat_relsummary('psyrat_data', base, args{:});

            headless = base;
            headless.proc.interactive = false;
            [outHeadless, errHeadless] = psyrat_relsummary('psyrat_data', headless, args{:});

            testCase.verifyEqual(errHeadless.nogooddata, errDefault.nogooddata);
            %on good data (no threshold failure) the computed dependability must
            %be identical: the interactive flag gates only the dialog, not math.
            if errDefault.nogooddata == 0
                testCase.verifyEqual( ...
                    outHeadless.relsummary.group(1).event(1).dep.m, ...
                    outDefault.relsummary.group(1).event(1).dep.m, ...
                    'AbsTol', 1e-12);
            end
        end

    end
end

% ---- local helpers ---------------------------------------------------------

function rel = local_rel(converged, ndivergent)
%Minimal rel struct carrying just the convergence fields the notice reads.
rel = struct();
rel.out.conv.converged = converged;
rel.out.conv.ndivergent = ndivergent;
end

function tbl = local_conv_table(ndiv)
%A convergence table in the {name,n_eff,r_hat} layout with the divergence count
%stamped in the header cell (column 4), exactly as psyrat_storeconv produces.
tbl = {'name','n_eff','r_hat'; 'mu',120,1.00; 'sig',140,1.01};
tbl{1,4} = ndiv;
end
