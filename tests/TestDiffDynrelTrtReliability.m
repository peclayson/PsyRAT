classdef TestDiffDynrelTrtReliability < PsyRATTestBase
    % Deterministic unit tests for the group-level trial + occasion two-facet
    % dynamic difference-score reliability calculator (psyrat_rel_diffdynrel_trt,
    % Rast & Clayson, analysis 15) and the summary isdifftrt branch
    % (psyrat_dynrel_summary). No CmdStan: the calculator is checked against the
    % closed-form two-facet gain-minus-loss coefficients (the same forms as
    % psyrat_diffrel_trt, every component a gain-minus-loss contrast variance,
    % with the per-event residual sigma_e,event(z) = exp(b_sigma[event] +
    % b_sigma_dim[event].z)) for all three reltypes, using single-draw inputs so
    % the credible interval collapses to the point estimate. Events are balanced
    % (obs1=obs2, nocc1=nocc2), so the D-study harmonic means reduce to the per-
    % facet contrast variance divided by the per-facet n'. The live known-truth
    % recovery is a separate slow test (TestDiffDynrelTrtRecovery).

    methods (Test)

        % ---- calculator: all three reltypes vs closed form -------------------

        function testReltype3OneDimMatchesClosedForm(testCase)
            p = localPars();
            z1 = [-2 -1 0 1 2];
            out = localCalc(p, 1, z1, [], 3);
            for a = 1:numel(z1)
                [wp1, wp2] = localResid(p, z1(a), []);
                Gcf = localClosed(3, 'gen', p, wp1, wp2);
                Dcf = localClosed(3, 'dep', p, wp1, wp2);
                testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
            end
            % single-draw inputs: CI bounds collapse to the point estimate
            testCase.verifyEqual(out.G.ll, out.G.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.G.ul, out.G.pt, 'AbsTol', 1e-12);
        end

        function testReltype1And2MatchClosedForm(testCase)
            p = localPars();
            z1 = [-1 0 1];
            for rt = [1 2]
                out = localCalc(p, 1, z1, [], rt);
                for a = 1:numel(z1)
                    [wp1, wp2] = localResid(p, z1(a), []);
                    Gcf = localClosed(rt, 'gen', p, wp1, wp2);
                    Dcf = localClosed(rt, 'dep', p, wp1, wp2);
                    testCase.verifyEqual(out.G.pt(a), Gcf, 'AbsTol', 1e-10);
                    testCase.verifyEqual(out.D.pt(a), Dcf, 'AbsTol', 1e-10);
                end
            end
        end

        function testDependabilityNotAboveGeneralizability(testCase)
            p = localPars();
            for rt = [1 2 3]
                out = localCalc(p, 1, -2:0.5:2, [], rt);
                testCase.verifyLessThanOrEqual(out.D.pt, out.G.pt + 1e-12);
            end
        end

        function testReliabilityDecreasesAsResidualGrows(testCase)
            % positive event-1 and event-2 scale slopes -> residual grows with z
            p = localPars(); p.bg = 0.4; p.bl = 0.4;
            out = localCalc(p, 1, -2:1:2, [], 3);
            testCase.verifyTrue(all(diff(out.G.pt) < 0));
        end

        % ---- occasion n' --------------------------------------------------

        function testMoreOccasionsRaisesReliability(testCase)
            p = localPars();
            out1 = localCalc(p, 1, 0, [], 3);        % nocc default 1
            p2 = p; p2.nocc = 4;
            out4 = localCalc(p2, 1, 0, [], 3);
            testCase.verifyGreaterThan(out4.G.pt, out1.G.pt);
        end

        function testNoccDefaultsToOne(testCase)
            p = localPars();
            ba = localBlockArgs(p);
            out = psyrat_rel_diffdynrel_trt(ba{:}, ...
                'b_sigma', [p.ls1 p.ls2], 'b_sigma_dim', [p.bg p.bl], ...
                'ndim', 1, 'z1', 0, 'obs', p.nt, 'reltype', 3, 'CI', .95);
            testCase.verifyEqual(out.nocc, 1);
        end

        % ---- two dimensions -------------------------------------------------

        function testTwoDimMatchesClosedFormAndShape(testCase)
            p = localPars();
            p.bsd = [0.2 -0.1 0.15 0.05 0.03 -0.02]; % [g_z1 l_z1 g_z2 l_z2 g_z1z2 l_z1z2]
            z1 = [-1 0 1]; z2 = [-1 1];
            out = localCalc(p, 2, z1, z2, 3);
            testCase.verifySize(out.G.pt, [numel(z1) numel(z2)]);
            for a = 1:numel(z1)
                for c = 1:numel(z2)
                    [wp1, wp2] = localResid(p, z1(a), z2(c));
                    Gcf = localClosed(3, 'gen', p, wp1, wp2);
                    testCase.verifyEqual(out.G.pt(a, c), Gcf, 'AbsTol', 1e-10);
                end
            end
        end

        % ---- validation -----------------------------------------------------

        function testRejectsWrongSlopeCount(testCase)
            p = localPars();
            ba = localBlockArgs(p);
            testCase.verifyError(@() psyrat_rel_diffdynrel_trt(ba{:}, ...
                'b_sigma', [p.ls1 p.ls2], 'b_sigma_dim', [1 2 3 4 5 6], ...
                'ndim', 1, 'z1', 0, 'obs', p.nt, 'reltype', 3, 'CI', .95), ...
                'varargin:bsigmadim');
        end

        function testRejectsBadReltype(testCase)
            p = localPars();
            testCase.verifyError(@() localCalc(p, 1, 0, [], 4), 'varargin:reltype');
        end

        % ---- summary isdifftrt branch (end to end, fabricated REL) ----------

        function testSummaryTrtDiffSurfaceMatchesClosedForm(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1);
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3, 'nocc', []);
            st = summ.strata(1);
            testCase.verifyEqual(st.obs, p.nt);   % per-event per-cell median
            testCase.verifyEqual(st.nocc, 1);     % owner default
            for a = 1:numel(st.z1)
                [wp1, wp2] = localResid(p, st.z1(a), []);
                pn = p; pn.nocc = 1;
                Gcf = localClosed(3, 'gen', pn, wp1, wp2);
                testCase.verifyEqual(st.G.pt(a), Gcf, 'AbsTol', 1e-8);
            end
        end

        function testSummaryObservedOccasions(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1);
            REL.out.nocc(1) = 5;
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            testCase.verifyEqual(summ.strata(1).nocc, 5);
        end

        function testSummaryVarcompTable(testCase)
            p = localPars();
            REL = localFabricateRel(p, 1);
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'reltype', 3);
            vc = summ.strata(1).varcomp;
            % the difference universe-score contrast variance must be present and
            % equal the hand value
            testCase.verifyTrue(any(strcmp(vc.symbol, 'sigma^2_p_delta')));
            uni = localDelta(p.idvc);
            testCase.verifyEqual(vc.estimate(strcmp(vc.symbol,'sigma^2_p_delta')), ...
                uni, 'AbsTol', 1e-8);
            % occasion difference variance row exists (two-facet specific)
            testCase.verifyTrue(any(strcmp(vc.symbol, 'sigma^2_o_delta')));
        end

    end
end

% ---- helpers ------------------------------------------------------------------

function p = localPars()
% Single-draw event 2x2 covariance blocks + per-event residual log-SDs + design.
p.idvc  = localBlock(9, 7, 2);     % person (-> bp)
p.trvc  = localBlock(1.5, 1.0, 0.3); % trial (-> bt)
p.ocvc  = localBlock(1.2, 0.9, 0.2); % occasion (-> bo)
p.tidvc = localBlock(1.0, 0.8, 0.25);% person x trial (-> bpi)
p.oidvc = localBlock(0.9, 0.7, 0.2); % person x occasion (-> bpo)
p.tovc  = localBlock(0.5, 0.4, 0.1); % trial x occasion (-> boi)
p.ls1 = log(3); p.ls2 = log(2.5);   % per-event residual log-SD at z = 0
p.bg = 0.2; p.bl = -0.1;            % event-specific scale slopes (one dimension)
p.bsd = [];                        % two-dimension slope row (set in 2-dim test)
p.nt = 8; p.nocc = 2;              % per-event trial / occasion n'
end

function M = localBlock(v11, v22, c12)
M = zeros(1,2,2); M(1,1,1) = v11; M(1,2,2) = v22; M(1,1,2) = c12; M(1,2,1) = c12;
end

function d = localDelta(M)
% gain-minus-loss contrast variance from a single-draw 2x2 block
d = M(1,1,1) + M(1,2,2) - 2*M(1,2,1);
end

function [wp1, wp2] = localResid(p, z1, z2)
% per-event residual VARIANCE at (z1[,z2]) (exp of the log-SD, squared)
if isempty(p.bsd)
    er1 = p.ls1 + p.bg*z1;
    er2 = p.ls2 + p.bl*z1;
else
    zt = [z1; z2; z1*z2];
    er1 = p.ls1 + [p.bsd(1) p.bsd(3) p.bsd(5)] * zt;
    er2 = p.ls2 + [p.bsd(2) p.bsd(4) p.bsd(6)] * zt;
end
wp1 = exp(er1)^2; wp2 = exp(er2)^2;
end

function args = localBlockArgs(p)
args = {'id_varcov', p.idvc, 'trl_varcov', p.trvc, 'occ_varcov', p.ocvc, ...
    'tid_varcov', p.tidvc, 'oid_varcov', p.oidvc, 'to_varcov', p.tovc};
end

function out = localCalc(p, ndim, z1, z2, reltype)
if isempty(p.bsd)
    bsd = [p.bg p.bl];
else
    bsd = p.bsd;
end
args = [localBlockArgs(p), {'b_sigma', [p.ls1 p.ls2], 'b_sigma_dim', bsd, ...
    'ndim', ndim, 'z1', z1, 'obs', p.nt, 'nocc', p.nocc, 'reltype', reltype, ...
    'CI', .95}];
if ~isempty(z2); args = [args, {'z2', z2}]; end
out = psyrat_rel_diffdynrel_trt(args{:});
end

function rel = localClosed(reltype, est, p, wp1, wp2)
% Closed-form two-facet gain-minus-loss coefficient, mirroring psyrat_diffrel_trt
% for balanced events (obs1=obs2=nt, nocc1=nocc2=no -> harmonic means = nt/no).
% Every component is the gain-minus-loss contrast variance / per-facet n'.
bp  = localDelta(p.idvc);
bpi = localDelta(p.tidvc);   % person x trial
bpo = localDelta(p.oidvc);   % person x occasion
bt  = localDelta(p.trvc);    % trial
bo  = localDelta(p.ocvc);    % occasion
boi = localDelta(p.tovc);    % trial x occasion
bpoi = wp1 + wp2;            % residual difference variance (er_cov = 0)
nt = p.nt; no = p.nocc;
switch reltype
    case 1 % equivalence (occasion fixed)
        uni = bp + bpo/no;
        rel_err = bpi/nt + bpoi/(nt*no);
        abs_err = rel_err + bt/nt + boi/(nt*no);
    case 2 % stability (trial fixed)
        uni = bp + bpi/nt;
        rel_err = bpo/no + bpoi/(nt*no);
        abs_err = rel_err + bo/no + boi/(nt*no);
    case 3 % equivalence and stability (both random)
        uni = bp;
        rel_err = bpi/nt + bpo/no + bpoi/(nt*no);
        abs_err = rel_err + bt/nt + bo/no + boi/(nt*no);
end
switch est
    case 'gen'; rel = uni / (uni + rel_err);
    case 'dep'; rel = uni / (uni + abs_err);
end
end

function REL = localFabricateRel(p, ndim)
% Fabricate a one-stratum ic_diff_dynrel_trt REL with single-draw (constant)
% posterior values matching the case-15 storage layout.
REL = struct();
REL.analysis = 'ic_diff_dynrel_trt';
REL.ndim = ndim;
if ndim == 1
    REL.dim_names = {'dim1'};
    bsd = [p.bg p.bl]; bdim = [1 0.5];
else
    REL.dim_names = {'dim1', 'dim2'};
    bsd = [p.bg p.bl 0.1 0.05 0 0]; bdim = [1 0.5 0.3 0.2 0 0];
end
REL.out = struct();
REL.out.labels = {'none'};
REL.out.dimz = {(-1.5:1:1.5)'};
REL.out.ntrials = p.nt;
REL.out.nocc = p.nocc;
REL.out.id_varcov  = {num2cell(p.idvc)};
REL.out.trl_varcov = {num2cell(p.trvc)};
REL.out.occ_varcov = {num2cell(p.ocvc)};
REL.out.tid_varcov = {num2cell(p.tidvc)};
REL.out.oid_varcov = {num2cell(p.oidvc)};
REL.out.to_varcov  = {num2cell(p.tovc)};
REL.out.b_sigma = {num2cell([p.ls1 p.ls2])};
REL.out.b_sigma_dim = {num2cell(bsd)};
REL.out.b_dim = {num2cell(bdim)};
end
