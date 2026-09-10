classdef TestSubjDiffDynrelTrtReliability < PsyRATTestBase
    % Deterministic unit tests for the SUBJECT-LEVEL trial + occasion two-facet
    % dynamic difference-score reliability (Rast & Clayson, analysis 16 / Stage 4b).
    % No CmdStan: psyrat_dynrel_summary is driven with a fabricated
    % ic_diff_dynrel_sserrvar_trt result whose posterior draws are constant (so
    % credible intervals collapse to the point estimate). The per-participant table
    % and the typical-person reference surface are checked against the validated
    % two-facet coefficient psyrat_diffrel_trt evaluated with the SAME inputs the
    % production path feeds it:
    %   - per subject s, at that subject's own z_s and own per-event residual
    %     er_e,s(z_s) = er_var_ss_e[s] + b_sigma_dim[e].z_s;
    %   - the reference surface, at the typical person (scale RE = 0), residual
    %     er_e(z) = b_sigma[e] + b_sigma_dim[e].z.
    % All three reltypes (CE/CS/CES) are exercised. This mirrors
    % TestSubjDiffDynrelReliability (case 13) and TestDiffDynrelTrtReliability
    % (case 15), composing the two.

    methods (Test)

        function testPerSubjectMatchesDiffrelTrt(testCase)
            % For every subject and every reltype, the per-participant gen/dep/ICC
            % equals psyrat_diffrel_trt at that subject's own z and residual, with
            % the default occasion n' = 1.
            REL = localFakeREL();
            blocks = localBlocks();
            bsdim = [0.2 -0.1];                 % event-specific scale slopes
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();

            for rt = [1 2 3]
                summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',rt);
                tab = summ.strata(1).ssrel_table;
                for s = 1:numel(zsub)
                    er = [erss1(s)+bsdim(1)*zsub(s), erss2(s)+bsdim(2)*zsub(s)];
                    g = psyrat_diffrel_trt(blocks{:},'er_var',er,...
                        'obs',[n1(s) n2(s)],'nocc',[1 1],'reltype',rt,...
                        'est','gen','er_cov',0,'CI',.95);
                    d = psyrat_diffrel_trt(blocks{:},'er_var',er,...
                        'obs',[n1(s) n2(s)],'nocc',[1 1],'reltype',rt,...
                        'est','dep','er_cov',0,'CI',.95);
                    row = find(tab.id2 == s,1);
                    testCase.verifyEqual(tab.gen_pt(row), g.pt, 'AbsTol', 1e-9);
                    testCase.verifyEqual(tab.dep_pt(row), d.pt, 'AbsTol', 1e-9);
                    testCase.verifyEqual(tab.iccg_pt(row), g.icc_pt, 'AbsTol', 1e-9);
                    testCase.verifyEqual(tab.iccd_pt(row), d.icc_pt, 'AbsTol', 1e-9);
                    % constant draws: the CI collapses to the point estimate
                    testCase.verifyEqual(tab.gen_ll(row), tab.gen_pt(row), 'AbsTol', 1e-8);
                end
            end
        end

        function testObservedOccasionsChangesNprime(testCase)
            % The 'observed' occasion option feeds each subject's stratum occasion
            % count (2 here) instead of the default 1; the per-subject coefficient
            % then matches psyrat_diffrel_trt at nocc = 2.
            REL = localFakeREL();
            blocks = localBlocks();
            bsdim = [0.2 -0.1];
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',3,...
                'nocc','observed');
            tab = summ.strata(1).ssrel_table;
            for s = 1:numel(zsub)
                er = [erss1(s)+bsdim(1)*zsub(s), erss2(s)+bsdim(2)*zsub(s)];
                d = psyrat_diffrel_trt(blocks{:},'er_var',er,...
                    'obs',[n1(s) n2(s)],'nocc',[2 2],'reltype',3,...
                    'est','dep','er_cov',0,'CI',.95);
                row = find(tab.id2 == s,1);
                testCase.verifyEqual(tab.dep_pt(row), d.pt, 'AbsTol', 1e-9);
            end
        end

        function testReferenceSurfaceMatchesTypicalPerson(testCase)
            % The reference surface (scale RE = 0) at the grid endpoints equals
            % psyrat_diffrel_trt at the typical-person residual b_sigma + slope.z,
            % using the stratum trial n' and the default occasion n' = 1.
            REL = localFakeREL();
            blocks = localBlocks();
            bsig = [log(3) log(2.5)]; bsdim = [0.2 -0.1];
            ntr = REL.out.ntrials(1);
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',3);
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                er = [bsig(1)+bsdim(1)*zv, bsig(2)+bsdim(2)*zv];
                g = psyrat_diffrel_trt(blocks{:},'er_var',er,'obs',[ntr ntr],...
                    'nocc',[1 1],'reltype',3,'est','gen','er_cov',0,'CI',.95);
                d = psyrat_diffrel_trt(blocks{:},'er_var',er,'obs',[ntr ntr],...
                    'nocc',[1 1],'reltype',3,'est','dep','er_cov',0,'CI',.95);
                testCase.verifyEqual(st.G.pt(idx), g.pt, 'AbsTol', 1e-9);
                testCase.verifyEqual(st.D.pt(idx), d.pt, 'AbsTol', 1e-9);
            end
        end

        function testDependabilityNotAboveGeneralizability(testCase)
            REL = localFakeREL();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',3);
            tab = summ.strata(1).ssrel_table;
            testCase.verifyLessThanOrEqual(tab.dep_pt, tab.gen_pt + 1e-12);
        end

        function testScaleVarcompRows(testCase)
            % The joint 4x4 person block reports the per-event scale SDs (sd_id[3:4])
            % and, with an identity person correlation, zero couplings; the
            % difference-component table rows are still present.
            REL = localFakeREL();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,'reltype',3);
            vc = summ.strata(1).varcomp;
            sd = localSdId();
            testCase.verifyEqual(localVal(vc,'sigma_delta_p,1'), sd(3), 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc,'sigma_delta_p,2'), sd(4), 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc,'rho_p1,sigma1'), 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc,'rho_sigma1,sigma2'), 0, 'AbsTol', 1e-10);
            testCase.verifyTrue(any(strcmp(vc.symbol,'sigma^2_p_delta')));
            % the occasion difference-component row from the two-facet table
            testCase.verifyTrue(any(strcmp(vc.symbol,'sigma^2_o_delta')));
        end

        function testSummaryProducesSubjectTableAndSurface(testCase)
            REL = localFakeREL();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',10,'reltype',3);
            testCase.verifyEqual(summ.strata(1).label, 'none');
            testCase.verifySize(summ.strata(1).G.pt, [1 10]);
            testCase.verifyEqual(summ.strata(1).nocc, 1); % default occasion n'
            tab = summ.strata(1).ssrel_table;
            testCase.verifyClass(tab, 'table');
            testCase.verifyEqual(height(tab), 3);
            for col = {'gen_pt','dep_pt','iccg_pt','iccd_pt','sem_gen_pt','sem_dep_pt'}
                testCase.verifyTrue(any(strcmp(tab.Properties.VariableNames, col{1})));
            end
        end

        function testFigureBuildsWithOverlay(testCase)
            REL = localFakeREL();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',10,'reltype',3);
            fh = psyrat_dynrelplot(summ,'gcoeff',2);
            closer = onCleanup(@() close(fh)); %#ok<NASGU>
            testCase.verifyClass(fh, 'matlab.ui.Figure');
        end

    end
end

% --- fixtures -----------------------------------------------------------------

function blocks = localBlocks()
% The six single-draw [1x2x2] event covariance blocks in psyrat_diffrel_trt order
% (bp/bpi/bpo/bt/bo/boi = id/tid/oid/trl/occ/to).
blocks = {'bp',localBlk(9,7,2),'bpi',localBlk(1.2,0.9,0.2),...
    'bpo',localBlk(0.7,0.5,0.15),'bt',localBlk(1.5,1.0,0.3),...
    'bo',localBlk(0.8,0.6,0.1),'boi',localBlk(0.4,0.3,0.05)};
end

function M = localBlk(v1,v2,cov)
M = zeros(1,2,2); M(1,1,1) = v1; M(1,2,2) = v2; M(1,1,2) = cov; M(1,2,1) = cov;
end

function z = localZ();      z = [-1; 0; 1];                 end
function e = localErss1();  e = [log(3); log(3.2); log(2.8)]; end
function e = localErss2();  e = [log(2.5); log(2.6); log(2.4)]; end
function t = localTrls1();  t = [20; 18; 22];               end
function t = localTrls2();  t = [20; 16; 24];               end
function s = localSdId();   s = [2.5 2.0 0.4 0.3];          end

function REL = localFakeREL()
% Minimal ic_diff_dynrel_sserrvar_trt result with constant posterior draws so the
% point estimates equal the closed form. 3 subjects, 1 dimension, 2 occasions,
% identity person correlation (zero location-scale couplings).
nd = 30; nsub = 3;
REL.analysis = 'ic_diff_dynrel_sserrvar_trt';
REL.ndim = 1; REL.dim_names = {'dim1'};

bl = localBlocks();
o = struct();
% expand each single-draw block to nd identical draws and store as cells
o.id_varcov  = {num2cell(repmat(bl{2},  nd, 1))};
o.tid_varcov = {num2cell(repmat(bl{4},  nd, 1))};
o.oid_varcov = {num2cell(repmat(bl{6},  nd, 1))};
o.trl_varcov = {num2cell(repmat(bl{8},  nd, 1))};
o.occ_varcov = {num2cell(repmat(bl{10}, nd, 1))};
o.to_varcov  = {num2cell(repmat(bl{12}, nd, 1))};

o.b_sigma     = {num2cell(repmat([log(3) log(2.5)], nd, 1))};
o.b_sigma_dim = {num2cell(repmat([0.2 -0.1], nd, 1))};
o.b_dim       = {num2cell(repmat([1 0.5], nd, 1))};

% per-subject z = 0 log-residuals (draws x NSUB), constant across draws
o.er_var_ss1 = {num2cell(repmat(localErss1()', nd, 1))};
o.er_var_ss2 = {num2cell(repmat(localErss2()', nd, 1))};

% joint person block: constant SDs and identity correlation cholesky
sd = localSdId();
o.sd_id = {num2cell(repmat(sd, nd, 1))};
Lid = zeros(nd, 4, 4);
for k = 1:nd
    Lid(k, :, :) = eye(4);
end
o.L_id = {num2cell(Lid)};

% per-subject lookup aligned to id2 (with per-event occasion counts)
id = (101:100+nsub)';
occs = [2; 2; 2];
ssinfo = table(id, (1:nsub)', localZ(), localTrls1(), localTrls2(), ...
    round((localTrls1()+localTrls2())/2), occs, occs, ...
    'VariableNames', {'id','id2','z1','trls1','trls2','trls','occs1','occs2'});
o.ssinfo = {ssinfo};

o.dimz = {localZ()};
o.ntrials = 20;
o.nocc = 2;            % observed occasions per stratum
o.labels = {'none'};
REL.out = o;
end

function v = localVal(vc, sym)
v = vc.estimate(strcmp(vc.symbol, sym));
end
