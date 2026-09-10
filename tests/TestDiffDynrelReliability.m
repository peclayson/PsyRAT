classdef TestDiffDynrelReliability < PsyRATTestBase
    % Deterministic unit tests for the group-level dynamic difference-score
    % reliability calculator (psyrat_rel_diffdynrel, Rast & Clayson, analysis
    % 12). No CmdStan: the calculator is checked against the closed-form
    % gain-minus-loss coefficients
    %   sigma^2_p_delta = bp1 + bp2 - 2*bp_cov
    %   G_delta(z) = sigma^2_p_delta / (sigma^2_p_delta + (sigma^2_e,1(z)+sigma^2_e,2(z))/n')
    %   D_delta(z) = G_delta with the trial difference variance added
    % with the per-event residual log-SD a function of the standardized
    % dimensions (sigma_e,event(z) = exp(b_sigma[event] + b_sigma_dim[event].z)).
    % Single-draw inputs collapse the credible interval to the point estimate.

    methods (Test)

        function testOneDimMatchesClosedForm(testCase)
            [idvc, trvc, uni] = localVarcov();
            ls1 = log(3); ls2 = log(2.5); nprime = 20;
            bg = 0.2; bl = -0.1;                    % event-specific scale slopes
            st1 = trvc(1,1,1); st2 = trvc(1,2,2); tcov = trvc(1,2,1);
            z1 = [-2 -1 0 1 2];
            out = psyrat_rel_diffdynrel('id_varcov', idvc, 'trl_varcov', trvc, ...
                'b_sigma', [ls1 ls2], 'b_sigma_dim', [bg bl], 'ndim', 1, ...
                'z1', z1, 'obs', nprime, 'CI', .95);
            for a = 1:numel(z1)
                wp1 = exp(ls1 + bg*z1(a))^2; wp2 = exp(ls2 + bl*z1(a))^2;
                rel = (wp1 + wp2) / nprime;
                abserr = rel + (st1 + st2)/nprime - 2*tcov/nprime;
                Gcf = uni / (uni + rel);
                Dcf = uni / (uni + abserr);
                testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
            end
            % single-draw: CI bounds collapse to the point estimate
            testCase.verifyEqual(out.G.ll, out.G.pt, 'AbsTol', 1e-12);
        end

        function testDependabilityNotAboveGeneralizability(testCase)
            [idvc, trvc] = localVarcov();
            out = psyrat_rel_diffdynrel('id_varcov', idvc, 'trl_varcov', trvc, ...
                'b_sigma', [log(3) log(2.5)], 'b_sigma_dim', [0.2 -0.1], ...
                'ndim', 1, 'z1', -2:0.5:2, 'obs', 20, 'CI', .95);
            testCase.verifyLessThanOrEqual(out.D.pt, out.G.pt + 1e-12);
        end

        function testTwoDimMatchesClosedFormAndShape(testCase)
            [idvc, trvc, uni] = localVarcov();
            ls1 = log(3); ls2 = log(2.5); nprime = 20;
            bsd = [0.2 -0.1 0.15 0.05 0.03 -0.02]; % [g_z1 l_z1 g_z2 l_z2 g_z1z2 l_z1z2]
            z1 = [-1 0 1]; z2 = [-1 1];
            out = psyrat_rel_diffdynrel('id_varcov', idvc, 'trl_varcov', trvc, ...
                'b_sigma', [ls1 ls2], 'b_sigma_dim', bsd, 'ndim', 2, ...
                'z1', z1, 'z2', z2, 'obs', nprime, 'CI', .95);
            testCase.verifySize(out.G.pt, [numel(z1) numel(z2)]);
            for a = 1:numel(z1)
                for c = 1:numel(z2)
                    zt = [z1(a); z2(c); z1(a)*z2(c)];
                    er1 = ls1 + [bsd(1) bsd(3) bsd(5)] * zt;
                    er2 = ls2 + [bsd(2) bsd(4) bsd(6)] * zt;
                    wp1 = exp(er1)^2; wp2 = exp(er2)^2;
                    Gcf = uni / (uni + (wp1 + wp2)/nprime);
                    testCase.verifyEqual(out.G.pt(a, c), Gcf, 'AbsTol', 1e-10);
                end
            end
        end

        function testRejectsWrongSlopeCount(testCase)
            [idvc, trvc] = localVarcov();
            % ndim=1 expects 2 slope columns; passing 6 should error
            testCase.verifyError(@() psyrat_rel_diffdynrel('id_varcov', idvc, ...
                'trl_varcov', trvc, 'b_sigma', [1 1], ...
                'b_sigma_dim', [1 2 3 4 5 6], 'ndim', 1, 'z1', 0, ...
                'obs', 10, 'CI', .95), 'varargin:bsigmadim');
        end

        function testSummaryDispatchesDifferenceVariant(testCase)
            % psyrat_dynrel_summary must accept the ic_diff_dynrel layout and
            % produce the same uniform surface struct with a sigma^2_p_delta row.
            REL = localFakeDiffREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 12);
            testCase.verifyEqual(summ.strata(1).label, 'none');
            testCase.verifySize(summ.strata(1).G.pt, [1 12]);
            vc = summ.strata(1).varcomp;
            testCase.verifyTrue(any(strcmp(vc.symbol, 'sigma^2_p_delta')));
        end

    end
end

function [idvc, trvc, uni] = localVarcov()
% Single-draw person/trial event variance-covariances.
idvc = zeros(1,2,2); idvc(1,1,1) = 9; idvc(1,2,2) = 7; idvc(1,2,1) = 2; idvc(1,1,2) = 2;
trvc = zeros(1,2,2); trvc(1,1,1) = 1.5; trvc(1,2,2) = 1.0; trvc(1,2,1) = 0.3; trvc(1,1,2) = 0.3;
uni = idvc(1,1,1) + idvc(1,2,2) - 2*idvc(1,2,1); % = 12
end

function REL = localFakeDiffREL()
% Minimal ic_diff_dynrel result with a handful of posterior draws.
nd = 50;
REL.analysis = 'ic_diff_dynrel'; REL.ndim = 1; REL.dim_names = {'dim1'};
idvc = zeros(nd,2,2); trvc = zeros(nd,2,2);
for k = 1:nd
    idvc(k,1,1) = 9; idvc(k,2,2) = 7; idvc(k,1,2) = 2; idvc(k,2,1) = 2;
    trvc(k,1,1) = 1.5; trvc(k,2,2) = 1.0; trvc(k,1,2) = 0.3; trvc(k,2,1) = 0.3;
end
o.id_varcov = {num2cell(idvc)};
o.trl_varcov = {num2cell(trvc)};
o.b_sigma = {num2cell(repmat([log(3) log(2.5)], nd, 1))};
o.b_sigma_dim = {num2cell(repmat([0.2 -0.1], nd, 1))};
o.b_dim = {num2cell(repmat([1 0.5], nd, 1))};
o.dimz = {linspace(-1.5, 1.5, 20)'};
o.ntrials = 20;
o.labels = {'none'};
REL.out = o;
end
