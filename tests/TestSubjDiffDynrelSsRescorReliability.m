classdef TestSubjDiffDynrelSsRescorReliability < PsyRATTestBase
    % Deterministic unit tests for the SUBJECT-LEVEL CONCURRENT dynamic
    % difference-score reliability with a PER-SUBJECT residual correlation (Rast &
    % Clayson, Stage 4e): one-facet (analysis 21, ic_diff_dynrel_sserrvar_rho) and
    % trial + occasion two-facet (analysis 22, ic_diff_dynrel_sserrvar_trt_rho). No
    % CmdStan: psyrat_dynrel_summary is driven with a fabricated result whose draws
    % are constant, a per-person scale random effect (er_var_ss), and a KNOWN
    % per-subject residual correlation vector rho_ss (one rho per participant).
    %
    % These models = the population-rescor models (cases 19/20) with the single
    % shared residual correlation replaced by a per-person rho. The per-participant
    % phi at each subject's own z is therefore checked against psyrat_diffrel /
    % psyrat_diffrel_trt with that subject's OWN z-conditioned per-event SDs AND
    % that subject's OWN residual covariance rho[s] * sigma_1(z,s) * sigma_2(z,s).
    % Setting every rho[s] = 0 must reproduce the non-concurrent case-13/16 table.
    % The typical-person reference surface uses the population-AVERAGE correlation
    % mean(rho).

    methods (Test)

        function testOneFacetPerSubjectRho(testCase)
            rho = localRhoVec();              % distinct per-subject correlations
            REL = localFakeOneFacet(rho);
            idvc = localIdvc(); trvc = localTrvc();
            bsdim = [0.2 -0.1];
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12);
            tab = summ.strata(1).ssrel_table;
            for s = 1:numel(zsub)
                er1 = erss1(s)+bsdim(1)*zsub(s);
                er2 = erss2(s)+bsdim(2)*zsub(s);
                wpcov = rho(s)*exp(er1)*exp(er2);   % each subject's OWN rho
                g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',[er1 er2],...
                    'obs',[n1(s) n2(s)],'est','gen','wp_cov',wpcov,'CI',.95);
                d = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',[er1 er2],...
                    'obs',[n1(s) n2(s)],'est','dep','wp_cov',wpcov,'CI',.95);
                row = find(tab.id2 == s,1);
                testCase.verifyEqual(tab.gen_pt(row), g.pt, 'AbsTol', 1e-9);
                testCase.verifyEqual(tab.dep_pt(row), d.pt, 'AbsTol', 1e-9);
            end
        end

        function testTwoFacetPerSubjectRhoAllReltypes(testCase)
            rho = localRhoVec();
            REL = localFakeTwoFacet(rho);
            blocks = localTrtBlocks();
            bsdim = [0.2 -0.1];
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();
            for rt = [1 2 3]
                summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12,...
                    'reltype',rt,'nocc','observed');
                tab = summ.strata(1).ssrel_table;
                for s = 1:numel(zsub)
                    er1 = erss1(s)+bsdim(1)*zsub(s);
                    er2 = erss2(s)+bsdim(2)*zsub(s);
                    ercov = rho(s)*exp(er1)*exp(er2);   % each subject's OWN rho
                    g = psyrat_diffrel_trt(blocks{:},'er_var',[er1 er2],...
                        'obs',[n1(s) n2(s)],'nocc',[2 2],'reltype',rt,...
                        'est','gen','er_cov',ercov,'CI',.95);
                    d = psyrat_diffrel_trt(blocks{:},'er_var',[er1 er2],...
                        'obs',[n1(s) n2(s)],'nocc',[2 2],'reltype',rt,...
                        'est','dep','er_cov',ercov,'CI',.95);
                    row = find(tab.id2 == s,1);
                    testCase.verifyEqual(tab.gen_pt(row), g.pt, 'AbsTol', 1e-9);
                    testCase.verifyEqual(tab.dep_pt(row), d.pt, 'AbsTol', 1e-9);
                end
            end
        end

        function testPerSubjectRhoZeroReducesToNonConcurrent(testCase)
            % every rho[s] = 0 -> the per-subject covariance drops out;
            % per-participant phi matches psyrat_diffrel with wp_cov = 0 (case 13).
            REL = localFakeOneFacet([0 0 0]);
            idvc = localIdvc(); trvc = localTrvc();
            bsdim = [0.2 -0.1];
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',8);
            tab = summ.strata(1).ssrel_table;
            for s = 1:numel(zsub)
                er1 = erss1(s)+bsdim(1)*zsub(s);
                er2 = erss2(s)+bsdim(2)*zsub(s);
                g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',[er1 er2],...
                    'obs',[n1(s) n2(s)],'est','gen','wp_cov',0,'CI',.95);
                row = find(tab.id2 == s,1);
                testCase.verifyEqual(tab.gen_pt(row), g.pt, 'AbsTol', 1e-9);
            end
        end

        function testReferenceSurfaceUsesMeanRho(testCase)
            % The typical-person reference surface (scale RE = 0) uses the
            % population-AVERAGE correlation mean(rho) and b_sigma + slope.z: G at
            % the grid endpoints matches psyrat_diffrel with that residual covariance.
            rho = localRhoVec();
            REL = localFakeOneFacet(rho);
            idvc = localIdvc(); trvc = localTrvc();
            bsig = [log(3) log(2.5)]; bsdim = [0.2 -0.1];
            ntr = REL.out.ntrials(1);
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',12);
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                er1 = bsig(1)+bsdim(1)*zv; er2 = bsig(2)+bsdim(2)*zv;
                wpcov = mean(rho)*exp(er1)*exp(er2);
                g = psyrat_diffrel('bp',idvc,'bt',trvc,'er_var',[er1 er2],...
                    'obs',[ntr ntr],'est','gen','wp_cov',wpcov,'CI',.95);
                testCase.verifyEqual(st.G.pt(idx), g.pt, 'AbsTol', 1e-9);
            end
        end

        function testVarcompRows(testCase)
            % Both variants report (a) the population-average residual correlation
            % (rho_e,12 = mean(rho)), (b) the NEW between-person SD of the residual
            % correlation (SD(rho_e,s) = std(rho)), and (c) the person scale-RE rows.
            rho = localRhoVec();
            for REL = {localFakeOneFacet(rho), localFakeTwoFacet(rho)}
                summ = psyrat_dynrel_summary(REL{1},'CI',.95,'ngrid',8);
                vc = summ.strata(1).varcomp;
                testCase.verifyEqual(localVal(vc,'rho_e,12'), mean(rho), 'AbsTol', 1e-10);
                testCase.verifyEqual(localVal(vc,'SD(rho_e,s)'), std(rho), 'AbsTol', 1e-10);
                sd = localSdId();
                testCase.verifyEqual(localVal(vc,'sigma_delta_p,1'), sd(3), 'AbsTol', 1e-10);
                testCase.verifyEqual(localVal(vc,'sigma_delta_p,2'), sd(4), 'AbsTol', 1e-10);
            end
        end

        function testFigureBuildsWithOverlay(testCase)
            REL = localFakeTwoFacet(localRhoVec());
            summ = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',8,'reltype',3);
            tab = summ.strata(1).ssrel_table;
            testCase.verifyClass(tab,'table');
            testCase.verifyEqual(height(tab),3);
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
blocks = {'bp',localBlk(9,7,2),'bpi',localBlk(1.2,0.9,0.2),...
    'bpo',localBlk(0.7,0.5,0.15),'bt',localBlk(1.5,1.0,0.3),...
    'bo',localBlk(0.8,0.6,0.1),'boi',localBlk(0.4,0.3,0.05)};
end

function z = localZ();      z = [-1; 0; 1];                 end
function e = localErss1();  e = [log(3); log(3.2); log(2.8)]; end
function e = localErss2();  e = [log(2.5); log(2.6); log(2.4)]; end
function t = localTrls1();  t = [20; 18; 22];               end
function t = localTrls2();  t = [20; 16; 24];               end
function s = localSdId();   s = [2.5 2.0 0.4 0.3];          end
function r = localRhoVec(); r = [0.3 0.5 0.6];              end % per-subject rho

function errvc = localErrvc(rhobar,nd)
% Typical-person residual covariance at z = 0: the population-AVERAGE correlation
% with the z = 0 per-event SDs (mirrors the generated-quantities convention).
s1 = 3; s2 = 2.5;
M = zeros(1,2,2);
M(1,1,1) = s1^2; M(1,2,2) = s2^2;
M(1,1,2) = rhobar*s1*s2; M(1,2,1) = rhobar*s1*s2;
errvc = repmat(M,nd,1);
end

function o = localCommonOut(rho,nd)
% shared REL.out fields (constant draws): slopes, person block, per-subject
% residuals, the PER-SUBJECT correlation rho_ss + its population summaries
% (rescor = mean(rho), err_varcov at mean(rho)), dimz, ntrials. (The fabricated
% fixtures use a fixed 3 subjects; the per-subject arrays match.)
nsub = numel(rho);
o = struct();
o.b_sigma     = {num2cell(repmat([log(3) log(2.5)],nd,1))};
o.b_sigma_dim = {num2cell(repmat([0.2 -0.1],nd,1))};
o.b_dim       = {num2cell(repmat([1 0.5],nd,1))};
o.er_var_ss1  = {num2cell(repmat(localErss1()',nd,1))};
o.er_var_ss2  = {num2cell(repmat(localErss2()',nd,1))};
sd = localSdId();
o.sd_id = {num2cell(repmat(sd,nd,1))};
Lid = zeros(nd,4,4);
for k = 1:nd; Lid(k,:,:) = eye(4); end
o.L_id = {num2cell(Lid)};
o.rho_ss     = {num2cell(repmat(rho(:)',nd,1))};      % draws x NSUB
o.rescor     = {num2cell(repmat(mean(rho),nd,1))};    % draws x 1 (mean rho)
o.rescor_mu  = {num2cell(repmat(atanh(mean(rho)),nd,1))};
o.sd_rescor  = {num2cell(repmat(0.5,nd,1))};
o.err_varcov = {num2cell(localErrvc(mean(rho),nd))};
o.dimz = {localZ()};
o.ntrials = 20;
o.labels = {'none'};
assert(nsub == 3, 'fixture assumes 3 subjects');
end

function REL = localFakeOneFacet(rho)
nd = 30; nsub = 3;
REL.analysis = 'ic_diff_dynrel_sserrvar_rho';
REL.ndim = 1; REL.dim_names = {'dim1'};
o = localCommonOut(rho,nd);
o.id_varcov  = {num2cell(repmat(localIdvc(),nd,1))};
o.trl_varcov = {num2cell(repmat(localTrvc(),nd,1))};
id = (101:100+nsub)';
ssinfo = table(id,(1:nsub)',localZ(),localTrls1(),localTrls2(),...
    round((localTrls1()+localTrls2())/2),...
    'VariableNames',{'id','id2','z1','trls1','trls2','trls'});
o.ssinfo = {ssinfo};
REL.out = o;
end

function REL = localFakeTwoFacet(rho)
nd = 30; nsub = 3;
REL.analysis = 'ic_diff_dynrel_sserrvar_trt_rho';
REL.ndim = 1; REL.dim_names = {'dim1'};
o = localCommonOut(rho,nd);
bl = localTrtBlocks();
o.id_varcov  = {num2cell(repmat(bl{2},  nd, 1))};
o.tid_varcov = {num2cell(repmat(bl{4},  nd, 1))};
o.oid_varcov = {num2cell(repmat(bl{6},  nd, 1))};
o.trl_varcov = {num2cell(repmat(bl{8},  nd, 1))};
o.occ_varcov = {num2cell(repmat(bl{10}, nd, 1))};
o.to_varcov  = {num2cell(repmat(bl{12}, nd, 1))};
o.nocc = 2;
id = (101:100+nsub)';
occs = [2; 2; 2];
ssinfo = table(id,(1:nsub)',localZ(),localTrls1(),localTrls2(),...
    round((localTrls1()+localTrls2())/2),occs,occs,...
    'VariableNames',{'id','id2','z1','trls1','trls2','trls','occs1','occs2'});
o.ssinfo = {ssinfo};
REL.out = o;
end

function v = localVal(vc, sym)
v = vc.estimate(strcmp(vc.symbol, sym));
end
