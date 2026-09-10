classdef TestExternalBenchmark < PsyRATTestBase
    % External validation against published (under-review) trial-count tables
    % from the effort-doors task (Heindorf, Rocha, Chen, Vispoel, & Clayson,
    % under review). See PsyRATHeindorfFixture for the source (co)variances,
    % published counts, and the caveats that motivate the tolerances below.
    %
    % Unlike PsyRATAccuracyOracle (which re-implements the formulas and so
    % tests internal consistency), these tests compare against numbers from a
    % separate analysis, so they validate correctness of the dependability,
    % difference-score, and difference-of-differences math plus the trial
    % cutoff logic.
    %
    % Tolerances are empirically grounded:
    %   - Single conditions reproduce the published counts to within 1 trial
    %     (30 of 32 exactly; the two that differ are at the 0.95 threshold and
    %     reflect two-decimal rounding of the published inputs).
    %   - Difference scores and the DoD contrast reproduce the published counts
    %     to a few percent. The point-estimate (single-draw) reliability used
    %     here sits slightly above the full-posterior mean reliability, so the
    %     computed counts run modestly below the published values. A relative
    %     tolerance covers this documented gap while still catching real
    %     formula errors, which would move counts far more.

    properties (Constant)
        CI = 0.95;
        SingleCountTol = 1;        % absolute trials
        DiffRelTol = 0.05;         % relative
        DodRelTol = 0.06;          % relative
    end

    methods (Test)
        function testHeindorfSingleConditionTrialCounts(testCase)
            comps = PsyRATHeindorfFixture.components();
            for k = 1:numel(comps)
                c = comps(k);
                for i = 1:numel(c.condition_order)
                    for j = 1:numel(c.thresholds)
                        n = PsyRATHeindorfFixture.singleConditionCount( ...
                            c.bp(i,i), c.bt(i,i), c.resid(i), ...
                            c.thresholds(j), testCase.CI);
                        pub = c.single(i,j);
                        testCase.verifyEqual(isnan(n), false, ...
                            sprintf('%s %s @%.2f: no crossing found', ...
                            c.name, c.condition_order{i}, c.thresholds(j)));
                        testCase.verifyLessThanOrEqual(abs(n - pub), testCase.SingleCountTol, ...
                            sprintf('%s %s @%.2f: computed %d vs published %d', ...
                            c.name, c.condition_order{i}, c.thresholds(j), n, pub));
                    end
                end
            end
        end

        function testHeindorfDifferenceScoreTrialCounts(testCase)
            comps = PsyRATHeindorfFixture.components();
            for k = 1:numel(comps)
                c = comps(k);
                pairs = { ...
                    [2 1], c.diff_LEgain_LEloss, 'LE-Gain - LE-Loss'; ...
                    [4 3], c.diff_HEgain_HEloss, 'HE-Gain - HE-Loss'};
                for p = 1:size(pairs, 1)
                    ix = pairs{p, 1};
                    pubRow = pairs{p, 2};
                    label = pairs{p, 3};
                    [bp2, bt2, erv2] = localTwoCondition(c, ix);
                    for j = 1:numel(c.thresholds)
                        n = PsyRATHeindorfFixture.diffScoreCount( ...
                            bp2, bt2, erv2, c.thresholds(j), testCase.CI);
                        pub = pubRow(j);
                        testCase.verifyEqual(isnan(n), false, ...
                            sprintf('%s %s @%.2f: no crossing found', ...
                            c.name, label, c.thresholds(j)));
                        testCase.verifyLessThanOrEqual(abs(n - pub) / pub, testCase.DiffRelTol, ...
                            sprintf('%s %s @%.2f: computed %d vs published %d (rel %.3f)', ...
                            c.name, label, c.thresholds(j), n, pub, abs(n - pub) / pub));
                    end
                end
            end
        end

        function testHeindorfDoDTrialCounts(testCase)
            comps = PsyRATHeindorfFixture.components();
            for k = 1:numel(comps)
                c = comps(k);
                bp4 = reshape(c.bp, [1 4 4]);
                bt4 = reshape(c.bt, [1 4 4]);
                erv4 = 0.5 * log(c.resid);
                for j = 1:numel(c.thresholds)
                    n = PsyRATHeindorfFixture.dodCount( ...
                        bp4, bt4, erv4, c.thresholds(j), testCase.CI);
                    pub = c.dod(j);
                    testCase.verifyEqual(isnan(n), false, ...
                        sprintf('%s DoD @%.2f: no crossing found', ...
                        c.name, c.thresholds(j)));
                    testCase.verifyLessThanOrEqual(abs(n - pub) / pub, testCase.DodRelTol, ...
                        sprintf('%s DoD @%.2f: computed %d vs published %d (rel %.3f)', ...
                        c.name, c.thresholds(j), n, pub, abs(n - pub) / pub));
                end
            end
        end
    end
end

function [bp2, bt2, erv2] = localTwoCondition(c, ix)
% Build [1 x 2 x 2] (co)variance draws and [1 x 2] log residual SD for the
% two constituent conditions of a difference score (indices ix = [a b]).
a = ix(1);
b = ix(2);
bp2 = reshape([c.bp(a,a) c.bp(a,b); c.bp(b,a) c.bp(b,b)], [1 2 2]);
bt2 = reshape([c.bt(a,a) c.bt(a,b); c.bt(b,a) c.bt(b,b)], [1 2 2]);
erv2 = 0.5 * log([c.resid(a) c.resid(b)]);
end
