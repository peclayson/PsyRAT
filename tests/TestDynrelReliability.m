classdef TestDynrelReliability < PsyRATTestBase
    % Deterministic unit tests for the dynamic/conditional reliability
    % (Rast & Clayson, analysis 11) subject-level dimension standardizer
    % (psyrat_zscore_subjectlevel) and reliability calculator
    % (psyrat_rel_dynrel). No CmdStan/MatlabStan required: the calculator is
    % checked against the closed-form one-facet coefficients
    %   G(z) = sigma_p^2 / (sigma_p^2 + sigma_e(z)^2 / n')
    %   D(z) = sigma_p^2 / (sigma_p^2 + (sigma_e(z)^2 + sigma_i^2) / n')
    % with single-draw inputs so the credible interval collapses to the point
    % estimate. The live known-truth CmdStan recovery is a separate slow test.

    methods (Test)

        % ---- standardizer ----------------------------------------------------

        function testStandardizerMeanZeroSdOne(testCase)
            T = localSubjTable([10 12 14 16 18], 4);
            z = psyrat_zscore_subjectlevel(T, 'dim1');
            zsub = z(1:4:end);  % first row of each of the 5 subjects
            testCase.verifyEqual(mean(zsub), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(std(zsub), 1, 'AbsTol', 1e-12);
        end

        function testStandardizerBroadcastsWithinSubject(testCase)
            T = localSubjTable([10 12 14 16 18], 4);
            z = psyrat_zscore_subjectlevel(T, 'dim1');
            % all four rows of subject 1 share the same standardized value
            testCase.verifyEqual(z(1:4), repmat(z(1), 4, 1), 'AbsTol', 0);
        end

        function testStandardizerReStandardizationIsIdentity(testCase)
            % Pre-standardized values (mean 0, sd 1) should be returned
            % unchanged, so simulating with standardized z then passing the same
            % values to the model yields the same z (used by the recovery test).
            raw = (-2:1)';                 % 4 distinct subject values
            raw = (raw - mean(raw)) / std(raw);
            T = localSubjTableVec(raw, 3);
            z = psyrat_zscore_subjectlevel(T, 'dim1');
            testCase.verifyEqual(z(1:3:end), raw, 'AbsTol', 1e-12);
        end

        function testStandardizerErrorsOnWithinIdVariation(testCase)
            T = localSubjTable([10 12 14 16 18], 4);
            T.dim1(1) = 99;  % subject 1 no longer constant
            testCase.verifyError(@() psyrat_zscore_subjectlevel(T, 'dim1'), ...
                'dimprep:notconstant');
        end

        function testStandardizerErrorsOnNoVariance(testCase)
            T = localSubjTable([7 7 7 7 7], 4);  % all subjects identical
            testCase.verifyError(@() psyrat_zscore_subjectlevel(T, 'dim1'), ...
                'dimprep:novariance');
        end

        % ---- calculator: one dimension --------------------------------------

        function testCalcOneDimMatchesClosedForm(testCase)
            sig_p = 4; sig_log0 = log(3); b_sigma = 0.25; sig_i = 2; nprime = 10;
            z1 = [-2 -1 0 1 2];
            out = psyrat_rel_dynrel('sig_p', sig_p, 'sig_log0', sig_log0, ...
                'b_sigma', b_sigma, 'sig_i', sig_i, 'ndim', 1, 'z1', z1, ...
                'obs', nprime, 'CI', .95);
            for a = 1:numel(z1)
                sige = exp(sig_log0 + b_sigma * z1(a));
                Gcf = sig_p^2 / (sig_p^2 + sige^2 / nprime);
                Dcf = sig_p^2 / (sig_p^2 + (sige^2 + sig_i^2) / nprime);
                iccg = sig_p^2 / (sig_p^2 + sige^2);
                testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.ICCg.pt(a), iccg, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.sige_pt(a), sige, 'AbsTol', 1e-10);
            end
            % single-draw inputs: CI bounds collapse to the point estimate
            testCase.verifyEqual(out.G.ll, out.G.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.G.ul, out.G.pt, 'AbsTol', 1e-12);
        end

        function testCalcDependabilityNotAboveGeneralizability(testCase)
            out = psyrat_rel_dynrel('sig_p', 4, 'sig_log0', log(3), ...
                'b_sigma', 0.25, 'sig_i', 2, 'ndim', 1, 'z1', -2:0.5:2, ...
                'obs', 10, 'CI', .95);
            testCase.verifyLessThanOrEqual(out.D.pt, out.G.pt + 1e-12);
        end

        function testCalcReliabilityDecreasesAsResidualGrows(testCase)
            % positive scale slope -> residual grows with z -> reliability falls
            out = psyrat_rel_dynrel('sig_p', 4, 'sig_log0', log(3), ...
                'b_sigma', 0.4, 'sig_i', 2, 'ndim', 1, 'z1', -2:1:2, ...
                'obs', 10, 'CI', .95);
            testCase.verifyTrue(all(diff(out.G.pt) < 0));
        end

        % ---- calculator: two dimensions -------------------------------------

        function testCalcTwoDimMatchesClosedFormAndShape(testCase)
            sig_p = 4; sig_log0 = log(3); sig_i = 2; nprime = 10;
            b_sigma = [0.2 -0.1 0.05];          % main1, main2, interaction
            z1 = [-1 0 1]; z2 = [-1 1];
            out = psyrat_rel_dynrel('sig_p', sig_p, 'sig_log0', sig_log0, ...
                'b_sigma', b_sigma, 'sig_i', sig_i, 'ndim', 2, 'z1', z1, ...
                'z2', z2, 'obs', nprime, 'CI', .95);
            testCase.verifySize(out.G.pt, [numel(z1) numel(z2)]);
            for a = 1:numel(z1)
                for c = 1:numel(z2)
                    xrow = [z1(a) z2(c) z1(a) * z2(c)];
                    sige = exp(sig_log0 + b_sigma * xrow');
                    Gcf = sig_p^2 / (sig_p^2 + sige^2 / nprime);
                    testCase.verifyEqual(out.G.pt(a, c), Gcf, 'AbsTol', 1e-10);
                end
            end
        end

        % ---- calculator: marginal (lognormal) option ------------------------

        function testCalcMarginalLognormal(testCase)
            sig_p = 4; sig_log0 = log(3); b_sigma = 0.25; sig_i = 2;
            nprime = 10; sig_dp = 0.5;
            out = psyrat_rel_dynrel('sig_p', sig_p, 'sig_log0', sig_log0, ...
                'b_sigma', b_sigma, 'sig_i', sig_i, 'ndim', 1, 'z1', 0, ...
                'obs', nprime, 'CI', .95, 'sig_dp', sig_dp, 'marginal', 1);
            sige_marg = exp(sig_log0) * exp(sig_dp^2);
            Gm = sig_p^2 / (sig_p^2 + sige_marg^2 / nprime);
            testCase.verifyEqual(out.G.pt(1), Gm, 'AbsTol', 1e-10);
        end

        % ---- calculator: input validation -----------------------------------

        function testCalcRejectsBadNdim(testCase)
            testCase.verifyError(@() psyrat_rel_dynrel('sig_p', 4, ...
                'sig_log0', 1, 'b_sigma', 0.2, 'sig_i', 2, 'ndim', 3, ...
                'z1', 0, 'obs', 10, 'CI', .95), 'varargin:ndim');
        end

        function testCalcRequiresZ2WhenTwoDim(testCase)
            testCase.verifyError(@() psyrat_rel_dynrel('sig_p', 4, ...
                'sig_log0', 1, 'b_sigma', [0.2 0.1 0], 'sig_i', 2, 'ndim', 2, ...
                'z1', 0, 'obs', 10, 'CI', .95), 'varargin:z2');
        end

        function testCalcRejectsWrongSlopeCount(testCase)
            % ndim=1 expects 1 slope column; passing 3 should error
            testCase.verifyError(@() psyrat_rel_dynrel('sig_p', 4, ...
                'sig_log0', 1, 'b_sigma', [0.2 0.1 0], 'sig_i', 2, 'ndim', 1, ...
                'z1', 0, 'obs', 10, 'CI', .95), 'varargin:bsigma');
        end

    end
end

function T = localSubjTable(subjvals, ntrl)
% Build a long table: numel(subjvals) subjects x ntrl trials, with a constant
% subject-level dim1 value per subject.
T = localSubjTableVec(subjvals(:), ntrl);
end

function T = localSubjTableVec(subjvals, ntrl)
ids = {}; dim1 = []; meas = [];
for s = 1:numel(subjvals)
    for t = 1:ntrl
        ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
        dim1(end+1, 1) = subjvals(s); %#ok<AGROW>
        meas(end+1, 1) = 0; %#ok<AGROW>
    end
end
T = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});
end
