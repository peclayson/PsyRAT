classdef TestDiffTrlCutoff < PsyRATTestBase
    % Tests for psyrat_diffrel_trlcutoff, which derives a SINGLE common trial
    % count at which a difference score reaches a reliability threshold.
    %
    % The difference-score cutoff row used to be evaluated at the two
    % INDEPENDENTLY derived per-event cutoffs, each chosen so that its own
    % event separately cleared the threshold, and the trial count was then
    % suppressed in the table. That answered a question about the two events
    % rather than about the difference score.
    %
    % With both counts held equal the coefficient has a closed form:
    % HM(n,n) = n, so the error variance is the single-observation variance of
    % the difference divided by n, and the coefficient is
    % uni / (uni + denom/n). The expected trial counts below are therefore
    % solved by hand from that identity rather than read back from the
    % production code -- these are literal anchors, not mirror tests.

    methods (Test)
        function testGeneralizabilityCutoffMatchesClosedForm(testCase)
            [bp, bt, er_var, wp_cov] = localFixture();

            % uni = 3, denom_rel = 16 + 1 - 2*2.8 = 11.4
            % 3/(3 + 11.4/n) >= 0.70  =>  0.9n >= 7.98  =>  n >= 8.867  =>  9
            [n, ds] = psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'gen', 'CI', 0.95, 'depcutoff', 0.70, ...
                'meascutoff', 2, 'ntrials', 1000);

            testCase.verifyEqual(n, 9);
            testCase.verifyEqual(ds.pt, 3/(3 + 11.4/9), 'AbsTol', 1e-12);
        end

        function testDependabilityCutoffMatchesClosedForm(testCase)
            [bp, bt, er_var, wp_cov] = localFixture();

            % denom_abs = 11.4 + (1 + 1 - 2*0.9) = 11.6
            % 3/(3 + 11.6/n) >= 0.70  =>  0.9n >= 8.12  =>  n >= 9.022  =>  10
            [n, ds] = psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'dep', 'CI', 0.95, 'depcutoff', 0.70, ...
                'meascutoff', 2, 'ntrials', 1000);

            testCase.verifyEqual(n, 10);
            testCase.verifyEqual(ds.pt, 3/(3 + 11.6/10), 'AbsTol', 1e-12);
        end

        function testReturnedCountIsTheSmallestThatClearsTheThreshold(testCase)
            % The defining contract, checked against production rather than
            % against the closed form: the coefficient reaches the threshold
            % at the returned count and does not reach it one trial earlier.
            [bp, bt, er_var, wp_cov] = localFixture();
            cutoff = 0.65;

            [n, ds] = psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'gen', 'CI', 0.95, 'depcutoff', cutoff, ...
                'meascutoff', 2, 'ntrials', 1000);

            testCase.verifyGreaterThan(n, 1);
            testCase.verifyGreaterThanOrEqual(ds.pt, cutoff);

            below = psyrat_diffrel('bp', bp, 'bt', bt, 'er_var', er_var, ...
                'wp_cov', wp_cov, 'obs', [n-1 n-1], 'est', 'gen', 'CI', 0.95);
            testCase.verifyLessThan(below.pt, cutoff);
        end

        function testReturnedScoreIsEvaluatedAtEqualTrialCounts(testCase)
            % The point of the change: the reported score must come from the
            % difference score at a COMMON n, not from a pair of per-event
            % counts. Equal counts also mean the guarded error variance never
            % fires, so the call is warning free.
            [bp, bt, er_var, wp_cov] = localFixture();

            [n, ds] = psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'dep', 'CI', 0.95, 'depcutoff', 0.60, ...
                'meascutoff', 2, 'ntrials', 1000);

            direct = testCase.verifyWarningFree(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', [n n], 'est', 'dep', 'CI', 0.95));

            testCase.verifyEqual(ds.pt, direct.pt, 'AbsTol', 0);
            testCase.verifyEqual(ds.ll, direct.ll, 'AbsTol', 0);
            testCase.verifyEqual(ds.ul, direct.ul, 'AbsTol', 0);
        end

        function testCredibleIntervalLimitNeedsMoreTrialsThanPointEstimate(testCase)
            % meascutoff selects which part of the credible interval has to
            % clear the threshold. With draws that actually vary, the lower
            % limit is the most demanding and the upper limit the least, so
            % the three cutoffs must be ordered. Constant draws could not
            % detect a meascutoff that was ignored.
            [bp, bt, er_var, wp_cov] = localVaryingFixture();

            args = {'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'gen', 'CI', 0.95, 'depcutoff', 0.60, 'ntrials', 5000};

            nLower = psyrat_diffrel_trlcutoff(args{:}, 'meascutoff', 1);
            nPoint = psyrat_diffrel_trlcutoff(args{:}, 'meascutoff', 2);
            nUpper = psyrat_diffrel_trlcutoff(args{:}, 'meascutoff', 3);

            testCase.verifyGreaterThan(nLower, nPoint);
            testCase.verifyGreaterThan(nPoint, nUpper);
        end

        function testUnreachableThresholdReturnsSentinel(testCase)
            % -1 is the same "not reached" sentinel the per-event trial
            % cutoffs use, and the caller is expected to fall back rather than
            % report a score, so no score is returned.
            [bp, bt, er_var, wp_cov] = localFixture();

            [n, ds] = psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'gen', 'CI', 0.95, 'depcutoff', 0.999999, ...
                'meascutoff', 2, 'ntrials', 20);

            testCase.verifyEqual(n, -1);
            testCase.verifyEmpty(ds);
        end

        function testInputValidation(testCase)
            [bp, bt, er_var, wp_cov] = localFixture();
            good = {'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'est', 'gen', 'CI', 0.95, 'depcutoff', 0.70, ...
                'meascutoff', 2, 'ntrials', 100};

            % odd number of inputs
            testCase.verifyError(@() psyrat_diffrel_trlcutoff(good{1:end-1}), ...
                'varargin:incomplete');

            % a required input missing entirely (depcutoff here)
            testCase.verifyError(@() psyrat_diffrel_trlcutoff( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'est', 'gen', ...
                'CI', 0.95, 'meascutoff', 2, 'ntrials', 100), ...
                'varargin:missing');

            badMeas = localReplace(good, 'meascutoff', 4);
            testCase.verifyError(@() psyrat_diffrel_trlcutoff(badMeas{:}), ...
                'varargin:meascutoff');

            badN = localReplace(good, 'ntrials', 0);
            testCase.verifyError(@() psyrat_diffrel_trlcutoff(badN{:}), ...
                'varargin:ntrials');
        end
    end
end

function [bp, bt, er_var, wp_cov] = localFixture()
%Constant-across-draws components, so every quantity has a closed form.
%uni = 2 + 2 - 2*0.5 = 3; residual variances 16 and 1 with a residual
%covariance of 2.8 (rho = 0.70); between-trial variances 1 and 1 with a
%covariance of 0.9.
nd = 6;
one = ones(nd,1);

bp = zeros(nd,2,2);
bp(:,1,1) = 2; bp(:,2,2) = 2; bp(:,1,2) = 0.5; bp(:,2,1) = 0.5;

bt = zeros(nd,2,2);
bt(:,1,1) = 1; bt(:,2,2) = 1; bt(:,1,2) = 0.9; bt(:,2,1) = 0.9;

er_var = log([4*one 1*one]);
wp_cov = 2.8*one;
end

function args = localReplace(args, name, value)
%Swap the value of one name/value pair, leaving the rest of the argument list
%alone, so each validation case differs from the valid call in exactly one way.
ind = find(strcmpi(name, args), 1);
args{ind+1} = value;
end

function [bp, bt, er_var, wp_cov] = localVaryingFixture()
%Draws that genuinely vary, so the lower limit, point estimate and upper
%limit of the credible interval separate.
bpDraws = [1.2; 1.6; 2.0; 2.4; 2.8; 3.2; 3.6; 4.0];
nd = numel(bpDraws);

bp = zeros(nd,2,2);
bp(:,1,1) = bpDraws;
bp(:,2,2) = bpDraws;
bp(:,1,2) = 0.4*bpDraws;
bp(:,2,1) = 0.4*bpDraws;

bt = zeros(nd,2,2);
bt(:,1,1) = 0.30; bt(:,2,2) = 0.30;
bt(:,1,2) = 0.05; bt(:,2,1) = 0.05;

er_var = log(repmat([2.0 2.0], nd, 1));
wp_cov = 0.5*ones(nd,1);
end
