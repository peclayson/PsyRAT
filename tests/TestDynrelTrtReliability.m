classdef TestDynrelTrtReliability < PsyRATTestBase
    % Deterministic unit tests for the trial + occasion two-facet
    % dynamic/conditional reliability (Rast & Clayson, analysis 14) calculator
    % (psyrat_rel_dynrel_trt) and the summary istrt branch
    % (psyrat_dynrel_summary). No CmdStan/MatlabStan required: the calculator is
    % checked against the closed-form two-facet coefficients (the same forms as
    % psyrat_rel_trt, with the residual sigma_pto,e(z) = exp(alpha + b_sigma.z))
    % for all three reltypes, using single-draw inputs so the credible interval
    % collapses to the point estimate. The live known-truth CmdStan recovery is a
    % separate slow test (TestDynrelTrtRecovery).

    methods (Test)

        % ---- calculator: one dimension, all three reltypes -------------------

        function testCalcReltype3OneDimMatchesClosedForm(testCase)
            p = localPars();
            z1 = [-2 -1 0 1 2];
            out = localCalc(p, 1, z1, [], 3);
            for a = 1:numel(z1)
                sige = exp(p.log0 + p.bsig * z1(a));
                Gcf = localTrtClosed(3, 2, p, sige);
                Dcf = localTrtClosed(3, 1, p, sige);
                testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.sige_pt(a), sige, 'AbsTol', 1e-10);
            end
            % single-draw inputs: CI bounds collapse to the point estimate
            testCase.verifyEqual(out.G.ll, out.G.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.G.ul, out.G.pt, 'AbsTol', 1e-12);
        end

        function testCalcReltype1And2MatchClosedForm(testCase)
            p = localPars();
            z1 = [-1 0 1];
            for rt = [1 2]
                out = localCalc(p, 1, z1, [], rt);
                for a = 1:numel(z1)
                    sige = exp(p.log0 + p.bsig * z1(a));
                    Gcf = localTrtClosed(rt, 2, p, sige);
                    Dcf = localTrtClosed(rt, 1, p, sige);
                    testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                    testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
                end
            end
        end

        function testCalcDependabilityNotAboveGeneralizability(testCase)
            p = localPars();
            for rt = [1 2 3]
                out = localCalc(p, 1, -2:0.5:2, [], rt);
                testCase.verifyLessThanOrEqual(out.D.pt, out.G.pt + 1e-12);
            end
        end

        function testCalcReliabilityDecreasesAsResidualGrows(testCase)
            % positive scale slope -> residual grows with z -> reliability falls
            p = localPars(); p.bsig = 0.4;
            out = localCalc(p, 1, -2:1:2, [], 3);
            testCase.verifyTrue(all(diff(out.G.pt) < 0));
        end

        % ---- occasion n' --------------------------------------------------

        function testMoreOccasionsRaisesReliability(testCase)
            % with both facets random, adding occasions shrinks the occasion and
            % residual error contributions, so reliability should not decrease.
            p = localPars();
            out1 = localCalc(p, 1, 0, [], 3); % nocc default 1
            p2 = p; p2.nocc = 4;
            out4 = localCalc(p2, 1, 0, [], 3);
            testCase.verifyGreaterThan(out4.G.pt, out1.G.pt);
        end

        function testNoccDefaultsToOne(testCase)
            p = localPars();
            out = psyrat_rel_dynrel_trt('sig_p', p.sig_p, 'sig_log0', p.log0, ...
                'b_sigma', p.bsig, 'sig_occ', p.sig_occ, 'sig_trl', p.sig_trl, ...
                'sig_trlxid', p.txp, 'sig_occxid', p.oxp, ...
                'sig_trlxocc', p.txo, 'ndim', 1, 'z1', 0, 'obs', p.nt, ...
                'reltype', 3, 'CI', .95); % nocc omitted
            testCase.verifyEqual(out.nocc, 1);
        end

        % ---- calculator: two dimensions -------------------------------------

        function testCalcTwoDimMatchesClosedFormAndShape(testCase)
            p = localPars(); p.bsig = [0.2 -0.1 0.05];
            z1 = [-1 0 1]; z2 = [-1 1];
            out = localCalc(p, 2, z1, z2, 3);
            testCase.verifySize(out.G.pt, [numel(z1) numel(z2)]);
            for a = 1:numel(z1)
                for c = 1:numel(z2)
                    xrow = [z1(a) z2(c) z1(a) * z2(c)];
                    sige = exp(p.log0 + p.bsig * xrow');
                    Gcf = localTrtClosed(3, 2, p, sige);
                    testCase.verifyEqual(out.G.pt(a, c), Gcf, 'AbsTol', 1e-10);
                end
            end
        end

        % ---- calculator: marginal (lognormal) option ------------------------

        function testCalcMarginalLognormal(testCase)
            p = localPars(); sig_dp = 0.5;
            out = psyrat_rel_dynrel_trt('sig_p', p.sig_p, 'sig_log0', p.log0, ...
                'b_sigma', p.bsig, 'sig_occ', p.sig_occ, 'sig_trl', p.sig_trl, ...
                'sig_trlxid', p.txp, 'sig_occxid', p.oxp, ...
                'sig_trlxocc', p.txo, 'ndim', 1, 'z1', 0, 'obs', p.nt, ...
                'nocc', p.nocc, 'reltype', 3, 'CI', .95, ...
                'sig_dp', sig_dp, 'marginal', 1);
            sige_marg = exp(p.log0) * exp(sig_dp^2);
            pm = p; % closed form with the marginal residual
            Gm = localTrtClosed(3, 2, pm, sige_marg);
            testCase.verifyEqual(out.G.pt(1), Gm, 'AbsTol', 1e-10);
        end

        % ---- calculator: input validation -----------------------------------

        function testCalcRejectsBadNdim(testCase)
            p = localPars();
            testCase.verifyError(@() localCalc(p, 3, 0, [], 3), 'varargin:ndim');
        end

        function testCalcRequiresZ2WhenTwoDim(testCase)
            p = localPars(); p.bsig = [0.2 0.1 0];
            testCase.verifyError(@() psyrat_rel_dynrel_trt('sig_p', p.sig_p, ...
                'sig_log0', p.log0, 'b_sigma', p.bsig, 'sig_occ', p.sig_occ, ...
                'sig_trl', p.sig_trl, 'sig_trlxid', p.txp, 'sig_occxid', p.oxp, ...
                'sig_trlxocc', p.txo, 'ndim', 2, 'z1', 0, 'obs', p.nt, ...
                'reltype', 3, 'CI', .95), 'varargin:z2');
        end

        function testCalcRejectsBadReltype(testCase)
            p = localPars();
            testCase.verifyError(@() localCalc(p, 1, 0, [], 4), 'varargin:reltype');
        end

        function testCalcRejectsBadNocc(testCase)
            p = localPars(); p.nocc = 0;
            testCase.verifyError(@() localCalc(p, 1, 0, [], 3), 'varargin:nocc');
        end

        % ---- summary istrt branch (end to end, fabricated REL) --------------

        function testSummaryTrtSurfaceMatchesClosedForm(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1); % one dimension
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3, ...
                'nocc', []); % default occasion n' = 1
            st = summ.strata(1);
            testCase.verifyEqual(st.obs, p.nt);   % per-(id,occasion) median
            testCase.verifyEqual(st.nocc, 1);     % owner default
            for a = 1:numel(st.z1)
                sige = exp(p.log0 + p.bsig * st.z1(a));
                pn = p; pn.nocc = 1;
                Gcf = localTrtClosed(3, 2, pn, sige);
                testCase.verifyEqual(st.G.pt(a), Gcf, 'AbsTol', 1e-8);
            end
        end

        function testSummaryObservedOccasions(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1);
            REL.out.nocc(1) = 5; % observed occasions in this stratum
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            testCase.verifyEqual(summ.strata(1).nocc, 5);
            % spot-check one grid point against the closed form at nocc = 5
            st = summ.strata(1);
            sige = exp(p.log0 + p.bsig * st.z1(1));
            pn = p; pn.nocc = 5;
            Gcf = localTrtClosed(3, 2, pn, sige);
            testCase.verifyEqual(st.G.pt(1), Gcf, 'AbsTol', 1e-8);
        end

        function testSummaryVarcompTable(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1);
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3);
            T = summ.strata(1).varcomp;
            % point estimates equal the (squared) fabricated SDs / values
            testCase.verifyEqual(localVc(T, 'sigma^2_p'), p.sig_p^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_o'), p.sig_occ^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_i'), p.sig_trl^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_pt'), p.txp^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_po'), p.oxp^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_to'), p.txo^2, 'AbsTol', 1e-8);
            testCase.verifyEqual(localVc(T, 'sigma^2_pto,e(0)'), exp(p.log0)^2, ...
                'AbsTol', 1e-8);
        end

    end
end

% ---- helpers ------------------------------------------------------------------

function p = localPars()
% A fixed set of two-facet variance components (SD units) + design counts.
p.sig_p   = 4;        % between-person SD
p.log0    = log(3);   % log residual SD at z = 0
p.bsig    = 0.25;     % scale dimension slope(s)
p.bmean   = 0.5;      % mean dimension slope(s)
p.sig_occ = 1.5;      % occasion main-effect SD
p.sig_trl = 2;        % trial main-effect SD
p.txp     = 1;        % trial x person SD
p.oxp     = 0.8;      % occasion x person SD
p.txo     = 0.5;      % trial x occasion SD
p.nt      = 8;        % trial n'
p.nocc    = 2;        % occasion n'
end

function out = localCalc(p, ndim, z1, z2, reltype)
args = {'sig_p', p.sig_p, 'sig_log0', p.log0, 'b_sigma', p.bsig, ...
    'sig_occ', p.sig_occ, 'sig_trl', p.sig_trl, 'sig_trlxid', p.txp, ...
    'sig_occxid', p.oxp, 'sig_trlxocc', p.txo, 'ndim', ndim, 'z1', z1, ...
    'obs', p.nt, 'nocc', p.nocc, 'reltype', reltype, 'CI', .95};
if ~isempty(z2); args = [args, {'z2', z2}]; end
out = psyrat_rel_dynrel_trt(args{:});
end

function rel = localTrtClosed(reltype, gcoeff, p, sige)
% Closed-form two-facet coefficient, mirroring psyrat_rel_trt exactly (variance
% components = squared SDs). gcoeff 1 = dependability (absolute), 2 =
% generalizability (relative). nt/no are the trial/occasion n'.
bp = p.sig_p^2; bo = p.sig_occ^2; bt = p.sig_trl^2;
txp = p.txp^2; oxp = p.oxp^2; txo = p.txo^2; err = sige^2;
nt = p.nt; no = p.nocc;
switch reltype
    case 1 % coefficient of equivalence (occasion fixed)
        uni = bp + oxp / no;
        if gcoeff == 1
            e = txp/nt + err/(nt*no) + bt/nt + txo/(nt*no);
        else
            e = txp/nt + err/(nt*no);
        end
    case 2 % coefficient of stability (trial fixed)
        uni = bp + txp / nt;
        if gcoeff == 1
            e = oxp/no + err/(nt*no) + bo/no + txo/(nt*no);
        else
            e = oxp/no + err/(nt*no);
        end
    case 3 % coefficient of trial equivalence and stability (both random)
        uni = bp;
        if gcoeff == 1
            e = txp/nt + oxp/no + err/(nt*no) + bt/nt + bo/no + txo/(nt*no);
        else
            e = txp/nt + oxp/no + err/(nt*no);
        end
end
rel = uni / (uni + e);
end

function REL = localFabricateRel(p, ndim)
% Fabricate a one-stratum ic_dynrel_trt REL with single-draw (constant) posterior
% values, matching the storage layout case 14 produces, so the summary istrt
% branch can be exercised without CmdStan. The credible interval collapses to the
% point estimate (one draw).
REL = struct();
REL.analysis = 'ic_dynrel_trt';
REL.ndim = ndim;
if ndim == 1
    REL.dim_names = {'dim1'};
    bsig_row = p.bsig; bmean_row = p.bmean;
else
    REL.dim_names = {'dim1', 'dim2'};
    bsig_row = [p.bsig p.bsig/2 0]; bmean_row = [p.bmean p.bmean/2 0];
end

% person location-scale cholesky factor with correlation rho = 0.3
rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

REL.out = struct();
REL.out.labels = {'none'};
REL.out.dimz = {(-1.5:1:1.5)'};        % 4 subject-level standardized values
REL.out.ntrials = p.nt;                % per-(id,occasion) median trial count
REL.out.nocc = p.nocc;                 % observed occasion count
REL.out.gro_sds = {num2cell([p.sig_p, 0.5])};   % (sigma_p, sigma_delta_p)
REL.out.pop_sdlog = p.log0;            % Intercept_sigma
REL.out.sig_occ = p.sig_occ;
REL.out.sig_trl = p.sig_trl;
REL.out.sig_trlxid = p.txp;
REL.out.sig_occxid = p.oxp;
REL.out.sig_trlxocc = p.txo;
REL.out.b = {num2cell(bmean_row)};
REL.out.b_sigma = {num2cell(bsig_row)};
REL.out.chol_corrmat = {num2cell(choldraw)};
end

function v = localVc(T, sym)
% Point estimate for a component symbol from a varcomp table.
v = T.estimate(strcmp(T.symbol, sym));
end
