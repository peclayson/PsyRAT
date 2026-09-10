classdef TestDiffDynrelRescorReliability < PsyRATTestBase
    % Deterministic unit tests for the CONCURRENT / residual-correlation (rescor)
    % GROUP-LEVEL dynamic difference-score reliability (Rast & Clayson, Stage 4c):
    % one-facet (analysis 17, ic_diff_dynrel_rescor) and trial + occasion two-facet
    % (analysis 18, ic_diff_dynrel_trt_rescor). No CmdStan: psyrat_dynrel_summary is
    % driven with a fabricated result whose posterior draws are constant (so credible
    % intervals collapse to the point estimate) and a KNOWN residual correlation rho.
    %
    % The reliability surface is checked against the validated coefficients
    % psyrat_diffrel / psyrat_diffrel_trt evaluated with the z-conditioned residual
    % covariance er_cov(z) = rho * sigma_1(z) * sigma_2(z), where
    % sigma_e(z) = exp(b_sigma[e] + b_sigma_dim[e].z). Setting rho = 0 must reproduce
    % the non-concurrent surface (cases 12/15). The residual-correlation varcomp row
    % is also verified.

    methods (Test)

        function testOneFacetConcurrentSurface(testCase)
            % Case 17: the surface at the grid endpoints equals psyrat_diffrel with
            % the z-conditioned residual covariance from the known rho.
            rho = 0.5;
            REL = localFakeOneFacet(rho);
            idvc = localIdvc(); trvc = localTrvc();
            bsig = [log(3) log(2.5)]; bsdim = [0.2 -0.1];
            ntr = REL.out.ntrials(1);
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12);
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                er = [bsig(1)+bsdim(1)*zv, bsig(2)+bsdim(2)*zv];
                ercov = rho * exp(er(1)) * exp(er(2));
                g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',er,'obs',[ntr ntr],...
                    'est','gen','wp_cov',ercov,'CI',.95);
                d = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',er,'obs',[ntr ntr],...
                    'est','dep','wp_cov',ercov,'CI',.95);
                testCase.verifyEqual(st.G.pt(idx), g.pt, 'AbsTol', 1e-9);
                testCase.verifyEqual(st.D.pt(idx), d.pt, 'AbsTol', 1e-9);
            end
        end

        function testTwoFacetConcurrentSurfaceAllReltypes(testCase)
            % Case 18: surface vs psyrat_diffrel_trt with the z-conditioned er_cov,
            % for all three reltypes; default occasion n' = 1.
            rho = 0.4;
            REL = localFakeTwoFacet(rho);
            blocks = localTrtBlocks();
            bsig = [log(3) log(2.5)]; bsdim = [0.2 -0.1];
            ntr = REL.out.ntrials(1);
            for rt = [1 2 3]
                summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',rt);
                st = summ.strata(1);
                for idx = [1 numel(st.z1)]
                    zv = st.z1(idx);
                    er = [bsig(1)+bsdim(1)*zv, bsig(2)+bsdim(2)*zv];
                    ercov = rho * exp(er(1)) * exp(er(2));
                    g = psyrat_diffrel_trt(blocks{:},'er_var',er,'obs',[ntr ntr],...
                        'nocc',[1 1],'reltype',rt,'est','gen','er_cov',ercov,'CI',.95);
                    d = psyrat_diffrel_trt(blocks{:},'er_var',er,'obs',[ntr ntr],...
                        'nocc',[1 1],'reltype',rt,'est','dep','er_cov',ercov,'CI',.95);
                    testCase.verifyEqual(st.G.pt(idx), g.pt, 'AbsTol', 1e-9);
                    testCase.verifyEqual(st.D.pt(idx), d.pt, 'AbsTol', 1e-9);
                end
            end
        end

        function testRescorZeroReducesToNonConcurrent(testCase)
            % With rho = 0 the concurrent surface equals the non-concurrent surface
            % (the rescor term drops out): G/D match psyrat_diffrel with wp_cov = 0.
            REL = localFakeOneFacet(0);
            idvc = localIdvc(); trvc = localTrvc();
            bsig = [log(3) log(2.5)]; bsdim = [0.2 -0.1];
            ntr = REL.out.ntrials(1);
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',8);
            st = summ.strata(1);
            for idx = 1:numel(st.z1)
                zv = st.z1(idx);
                er = [bsig(1)+bsdim(1)*zv, bsig(2)+bsdim(2)*zv];
                g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',er,'obs',[ntr ntr],...
                    'est','gen','wp_cov',0,'CI',.95);
                testCase.verifyEqual(st.G.pt(idx), g.pt, 'AbsTol', 1e-9);
            end
        end

        function testRescorVarcompRow(testCase)
            % Both variants append the residual-correlation rows: rho_e,12 equals
            % the known rho and sigma_e,12(0) equals err_varcov(2,1).
            rho = 0.5;
            for REL = {localFakeOneFacet(rho), localFakeTwoFacet(rho)}
                summ = psyrat_dynrel_summary(REL{1},'CI',.95,'ngrid',8);
                vc = summ.strata(1).varcomp;
                testCase.verifyEqual(localVal(vc,'rho_e,12'), rho, 'AbsTol', 1e-10);
                cov0true = rho * 3 * 2.5; % exp(log 3) * exp(log 2.5)
                testCase.verifyEqual(localVal(vc,'sigma_e,12(0)'), cov0true, 'AbsTol', 1e-9);
            end
        end

        function testFigureBuildsGroupLevel(testCase)
            % Group-level: the figure builds with no per-participant overlay.
            REL = localFakeTwoFacet(0.4);
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',8,'reltype',3);
            testCase.verifyEmpty(summ.strata(1).ssrel_table);
            fh = psyrat_dynrelplot(summ,'gcoeff',2);
            closer = onCleanup(@() close(fh)); %#ok<NASGU>
            testCase.verifyClass(fh, 'matlab.ui.Figure');
        end

    end
end

% --- fixtures -----------------------------------------------------------------

function M = localBlk(v1,v2,cov)
M = zeros(1,2,2); M(1,1,1) = v1; M(1,2,2) = v2; M(1,1,2) = cov; M(1,2,1) = cov;
end

function idvc = localIdvc(); idvc = localBlk(9,7,2); end
function trvc = localTrvc(); trvc = localBlk(1.5,1.0,0.3); end

function blocks = localTrtBlocks()
% bp/bpi/bpo/bt/bo/boi = id/tid/oid/trl/occ/to (single-draw [1x2x2]).
blocks = {'bp',localBlk(9,7,2),'bpi',localBlk(1.2,0.9,0.2),...
    'bpo',localBlk(0.7,0.5,0.15),'bt',localBlk(1.5,1.0,0.3),...
    'bo',localBlk(0.8,0.6,0.1),'boi',localBlk(0.4,0.3,0.05)};
end

function z = localZ(); z = [-1; 0; 1]; end

function errvc = localErrvc(rho,nd)
% residual covariance at z = 0 = quad_form_diag([1 rho; rho 1], exp(b_sigma)),
% with b_sigma = [log 3, log 2.5] -> SDs 3 and 2.5.
s1 = 3; s2 = 2.5;
M = zeros(1,2,2);
M(1,1,1) = s1^2; M(1,2,2) = s2^2;
M(1,1,2) = rho*s1*s2; M(1,2,1) = rho*s1*s2;
errvc = repmat(M,nd,1);
end

function REL = localFakeOneFacet(rho)
% Minimal ic_diff_dynrel_rescor result with constant draws and a fixed rho.
nd = 30;
REL.analysis = 'ic_diff_dynrel_rescor';
REL.ndim = 1; REL.dim_names = {'dim1'};
o = struct();
o.id_varcov  = {num2cell(repmat(localIdvc(),nd,1))};
o.trl_varcov = {num2cell(repmat(localTrvc(),nd,1))};
o.b_sigma     = {num2cell(repmat([log(3) log(2.5)],nd,1))};
o.b_sigma_dim = {num2cell(repmat([0.2 -0.1],nd,1))};
o.b_dim       = {num2cell(repmat([1 0.5],nd,1))};
o.rescor      = {num2cell(repmat(rho,nd,1))};
o.err_varcov  = {num2cell(localErrvc(rho,nd))};
o.dimz = {localZ()};
o.ntrials = 20;
o.labels = {'none'};
REL.out = o;
end

function REL = localFakeTwoFacet(rho)
% Minimal ic_diff_dynrel_trt_rescor result with constant draws and a fixed rho.
nd = 30;
REL.analysis = 'ic_diff_dynrel_trt_rescor';
REL.ndim = 1; REL.dim_names = {'dim1'};
bl = localTrtBlocks();
o = struct();
o.id_varcov  = {num2cell(repmat(bl{2},  nd, 1))};
o.tid_varcov = {num2cell(repmat(bl{4},  nd, 1))};
o.oid_varcov = {num2cell(repmat(bl{6},  nd, 1))};
o.trl_varcov = {num2cell(repmat(bl{8},  nd, 1))};
o.occ_varcov = {num2cell(repmat(bl{10}, nd, 1))};
o.to_varcov  = {num2cell(repmat(bl{12}, nd, 1))};
o.b_sigma     = {num2cell(repmat([log(3) log(2.5)],nd,1))};
o.b_sigma_dim = {num2cell(repmat([0.2 -0.1],nd,1))};
o.b_dim       = {num2cell(repmat([1 0.5],nd,1))};
o.rescor      = {num2cell(repmat(rho,nd,1))};
o.err_varcov  = {num2cell(localErrvc(rho,nd))};
o.dimz = {localZ()};
o.ntrials = 20;
o.nocc = 2;
o.labels = {'none'};
REL.out = o;
end

function v = localVal(vc, sym)
v = vc.estimate(strcmp(vc.symbol, sym));
end
