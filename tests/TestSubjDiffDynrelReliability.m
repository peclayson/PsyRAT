classdef TestSubjDiffDynrelReliability < PsyRATTestBase
    % Deterministic unit tests for the SUBJECT-LEVEL dynamic difference-score
    % reliability (Rast & Clayson, analysis 13). No CmdStan: psyrat_dynrel_summary
    % is driven with a fabricated ic_diff_dynrel_sserrvar result whose posterior
    % draws are constant (so credible intervals collapse to the point estimate),
    % and the per-participant table it produces is checked against the closed-form
    % gain-minus-loss coefficients evaluated at EACH SUBJECT's own standardized
    % dimension value z_s, using that subject's own per-event log-residual:
    %   sigma^2_p_delta = bp1 + bp2 - 2*bp_cov
    %   wp_e,s          = exp(er_var_ss_e[s] + b_sigma_dim[e].z_s)^2
    %   phi_rho,s   = uni / (uni + wp1/n1 + wp2/n2)                 (generalizability)
    %   phi_delta,s = uni / (uni + wp1/n1 + wp2/n2
    %                          + bt1/n1 + bt2/n2 - 2*bt_cov/harm(n)) (dependability)
    % The per-person scale random effect is what makes this a subject-level
    % quantity; psyrat_ssrel_diff is reused unchanged.

    methods (Test)

        function testPerSubjectMatchesClosedForm(testCase)
            REL = localFakeSubjREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 12);
            tab = summ.strata(1).ssrel_table;

            [idvc, trvc, uni] = localTruth();
            bsdim = [0.2 -0.1];                 % event-specific scale slopes
            erss1 = localErss1(); erss2 = localErss2();
            zsub = localZ(); n1 = localTrls1(); n2 = localTrls2();
            bt1 = trvc(1,1,1); bt2 = trvc(1,2,2); btc = trvc(1,2,1);

            for s = 1:numel(zsub)
                wp1 = exp(erss1(s) + bsdim(1)*zsub(s))^2;
                wp2 = exp(erss2(s) + bsdim(2)*zsub(s))^2;
                relerr = wp1/n1(s) + wp2/n2(s);
                harm = 2/(1/n1(s) + 1/n2(s));
                abserr = relerr + bt1/n1(s) + bt2/n2(s) - 2*btc/harm;
                genCf = uni / (uni + relerr);
                depCf = uni / (uni + abserr);
                iccgCf = uni / (uni + wp1 + wp2);
                iccdCf = uni / (uni + wp1 + wp2 + bt1 + bt2 - 2*btc);

                row = find(tab.id2 == s, 1);
                testCase.verifyEqual(tab.gen_pt(row), genCf, 'AbsTol', 1e-10);
                testCase.verifyEqual(tab.dep_pt(row), depCf, 'AbsTol', 1e-10);
                testCase.verifyEqual(tab.iccg_pt(row), iccgCf, 'AbsTol', 1e-10);
                testCase.verifyEqual(tab.iccd_pt(row), iccdCf, 'AbsTol', 1e-10);
                % constant draws: CI collapses to the point estimate
                testCase.verifyEqual(tab.gen_ll(row), tab.gen_pt(row), 'AbsTol', 1e-9);
            end
        end

        function testDependabilityNotAboveGeneralizability(testCase)
            REL = localFakeSubjREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 12);
            tab = summ.strata(1).ssrel_table;
            testCase.verifyLessThanOrEqual(tab.dep_pt, tab.gen_pt + 1e-12);
        end

        function testScaleVarcompRows(testCase)
            % The joint person block reports the per-event scale SDs (sd_id[3:4])
            % and, with an identity person correlation, zero location-scale and
            % scale-scale correlations.
            REL = localFakeSubjREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 12);
            vc = summ.strata(1).varcomp;

            sd = localSdId();   % [loc1 loc2 scale1 scale2]
            testCase.verifyEqual(localVal(vc, 'sigma_delta_p,1'), sd(3), 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc, 'sigma_delta_p,2'), sd(4), 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc, 'rho_p1,sigma1'), 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc, 'rho_p2,sigma2'), 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(localVal(vc, 'rho_sigma1,sigma2'), 0, 'AbsTol', 1e-10);
            % the difference universe-score variance row is still present
            testCase.verifyTrue(any(strcmp(vc.symbol, 'sigma^2_p_delta')));
        end

        function testFigureBuildsWithOverlay(testCase)
            % The figure builds headlessly and overlays the per-participant
            % points on the typical-person reference curve.
            REL = localFakeSubjREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 10);
            fh = psyrat_dynrelplot(summ, 'gcoeff', 2);
            closer = onCleanup(@() close(fh));
            testCase.verifyClass(fh, 'matlab.ui.Figure');
        end

        function testSummaryProducesSubjectTable(testCase)
            % The subject-level variant produces a per-participant table AND keeps
            % the typical-person reference surface (same uniform struct).
            REL = localFakeSubjREL();
            summ = psyrat_dynrel_summary(REL, 'CI', .95, 'ngrid', 10);
            testCase.verifyEqual(summ.strata(1).label, 'none');
            testCase.verifySize(summ.strata(1).G.pt, [1 10]);
            tab = summ.strata(1).ssrel_table;
            testCase.verifyClass(tab, 'table');
            testCase.verifyEqual(height(tab), 3);
            for col = {'gen_pt','dep_pt','iccg_pt','iccd_pt','sem_gen_pt','sem_dep_pt'}
                testCase.verifyTrue(any(strcmp(tab.Properties.VariableNames, col{1})));
            end
        end

    end
end

% --- fixtures -----------------------------------------------------------------

function [idvc, trvc, uni] = localTruth()
% Single-draw person/trial event variance-covariances (location).
idvc = zeros(1,2,2); idvc(1,1,1) = 9; idvc(1,2,2) = 7; idvc(1,1,2) = 2; idvc(1,2,1) = 2;
trvc = zeros(1,2,2); trvc(1,1,1) = 1.5; trvc(1,2,2) = 1.0; trvc(1,1,2) = 0.3; trvc(1,2,1) = 0.3;
uni = idvc(1,1,1) + idvc(1,2,2) - 2*idvc(1,2,1); % = 12
end

function z = localZ();      z = [-1; 0; 1];                 end
function e = localErss1();  e = [log(3); log(3.2); log(2.8)]; end
function e = localErss2();  e = [log(2.5); log(2.6); log(2.4)]; end
function t = localTrls1();  t = [20; 18; 22];               end
function t = localTrls2();  t = [20; 16; 24];               end
function s = localSdId();   s = [2.5 2.0 0.4 0.3];          end

function REL = localFakeSubjREL()
% Minimal ic_diff_dynrel_sserrvar result with constant posterior draws so the
% point estimates equal the closed form. 3 subjects, 1 dimension, identity person
% correlation (zero location-scale couplings).
nd = 30; nsub = 3;
REL.analysis = 'ic_diff_dynrel_sserrvar';
REL.ndim = 1; REL.dim_names = {'dim1'};

[idvc1, trvc1] = localTruth();
idvc = repmat(idvc1, nd, 1); trvc = repmat(trvc1, nd, 1);

o = struct();
o.id_varcov  = {num2cell(idvc)};
o.trl_varcov = {num2cell(trvc)};
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

% per-subject lookup aligned to id2
id = (101:100+nsub)';
ssinfo = table(id, (1:nsub)', localZ(), localTrls1(), localTrls2(), ...
    round((localTrls1()+localTrls2())/2), ...
    'VariableNames', {'id','id2','z1','trls1','trls2','trls'});
o.ssinfo = {ssinfo};

o.dimz = {localZ()};
o.ntrials = 20;
o.labels = {'none'};
REL.out = o;
end

function v = localVal(vc, sym)
v = vc.estimate(strcmp(vc.symbol, sym));
end
