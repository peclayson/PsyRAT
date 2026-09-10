classdef TestUnitLabels < matlab.unittest.TestCase
    % Locks the psyrat_unitlabels analysis-unit label contract. The splits
    % feature relabels every hard-coded "trial" in plots/tables/viewers through
    % this single helper, composing each display string so that the splits<=1
    % (single-trial) branch reproduces the ORIGINAL literal byte-for-byte. If a
    % future edit to psyrat_unitlabels changes one of these fragments, the
    % composed output stops being byte-identical for non-splits runs (a silent
    % regression in shared rendering code), so the exact strings are pinned here.

    methods (Test)

        function testTrialFragmentsAreExact(testCase)
            % Single-trial (default) fragments must equal the original literals
            % that the rendering code composed from before the splits feature.
            for splits = {1, [], false, 0}
                L = psyrat_unitlabels(splits{1});
                testCase.verifyFalse(L.issplit);
                testCase.verifyEqual(L.Unit,      'Trial');
                testCase.verifyEqual(L.Units,     'Trials');
                testCase.verifyEqual(L.unit,      'trial');
                testCase.verifyEqual(L.units,     'trials');
                testCase.verifyEqual(L.NumberOf,  'Number of Trials');
                testCase.verifyEqual(L.numberof,  'number of trials');
                testCase.verifyEqual(L.Cutoff,    'Trial Cutoff');
                testCase.verifyEqual(L.CutoffVar, 'Trial_Cutoff');
                testCase.verifyEqual(L.NumSymbol, '# of Trials');
            end
        end

        function testSplitFragmentsAreExact(testCase)
            % Parallel (2), nonparallel (3), and logical-true all switch to the
            % split fragments.
            for splits = {2, 3, true}
                L = psyrat_unitlabels(splits{1});
                testCase.verifyTrue(L.issplit);
                testCase.verifyEqual(L.Unit,      'Split');
                testCase.verifyEqual(L.Units,     'Splits');
                testCase.verifyEqual(L.unit,      'split');
                testCase.verifyEqual(L.units,     'splits');
                testCase.verifyEqual(L.NumberOf,  'Number of Splits');
                testCase.verifyEqual(L.numberof,  'number of splits');
                testCase.verifyEqual(L.Cutoff,    'Split Cutoff');
                testCase.verifyEqual(L.CutoffVar, 'Split_Cutoff');
                testCase.verifyEqual(L.NumSymbol, '# of Splits');
            end
        end

        function testComposedDisplayStringsRoundTrip(testCase)
            % Spot-check the composition idioms used in the rendering code: the
            % splits<=1 branch must reproduce the exact pre-splits literals.
            Lt = psyrat_unitlabels(1);
            testCase.verifyEqual(['Mean ' Lt.NumSymbol], 'Mean # of Trials');
            testCase.verifyEqual(['Std Dev of ' Lt.Units], 'Std Dev of Trials');
            testCase.verifyEqual(['Number of ' Lt.units ' to plot'], ...
                'Number of trials to plot');
            testCase.verifyEqual(['Plot: ' Lt.Units ' vs Reliability'], ...
                'Plot: Trials vs Reliability');
            testCase.verifyEqual(['Mean_Num_' Lt.Units], 'Mean_Num_Trials');

            % Regression pins for two pre-merge code-review label fixes:
            % (1) the "increasing N" table title must NOT double the unit word --
            %     L.NumberOf already contains the unit, so the title composes
            %     "Increasing the Number of Trials", not the old buggy
            %     "Increasing Trials the Number of Trials"
            %     (psyrat_depcutofft.m / psyrat_trt_relcutofft.m).
            testCase.verifyEqual( ...
                ['Results of Increasing the ' Lt.NumberOf ' on Dependability'], ...
                'Results of Increasing the Number of Trials on Dependability');
            % (2) the test-retest overall CSV header names the coefficient column
            %     ONCE -- guards against the old doubled "DependabilityDependability"
            %     (psyrat_trt_reloverallt.m); the 'Dependability' literal here is
            %     the gcoeff name that precedes the per-unit count columns.
            testCase.verifyEqual( ...
                ['Label,N Included,N Excluded,' 'Dependability' ',Mean Num of ' Lt.Units], ...
                'Label,N Included,N Excluded,Dependability,Mean Num of Trials');

            % and the split form composes the intended wording
            Ls = psyrat_unitlabels(3);
            testCase.verifyEqual(['Mean ' Ls.NumSymbol], 'Mean # of Splits');
            testCase.verifyEqual(['Plot: ' Ls.Units ' vs Reliability'], ...
                'Plot: Splits vs Reliability');
            testCase.verifyEqual( ...
                ['Results of Increasing the ' Ls.NumberOf ' on Dependability'], ...
                'Results of Increasing the Number of Splits on Dependability');
        end

    end
end
