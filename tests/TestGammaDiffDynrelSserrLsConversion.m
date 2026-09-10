classdef TestGammaDiffDynrelSserrLsConversion < matlab.unittest.TestCase
    %Offline checks for the LOCATION-SCALE gamma SUBJECT-LEVEL DYNAMIC DIFFERENCE
    %converter (analysis 13, gammascale = 2): psyrat_ssrel_diffdynrel_gamma_ls.
    %
    %Increment 8c. No CmdStan required - every test drives the converter directly
    %from synthetic posterior draws, so what is pinned is the conversion
    %arithmetic and the routing, NOT that the estimand is correct. Live recovery
    %(TestGammaDiffDynrelSserrLsRecovery) is what checks correctness.
    %
    %Organized around the claims the derivation in
    %documentation/gamma_scale_submodel_decision.md section J actually makes,
    %as amended by S17-FIX (the pooled per-participant residual):
    %  J.2  the CONDITIONAL part of the per-participant residual is a pure
    %       pass-through of the fitted per-person log residual SD, conditioned
    %       at that participant's own z
    %  J.3  the SIGNAL is a population quantity at the participant's z, NOT at
    %       that participant's own amplitude
    %  J.4  psyrat_gamma_crosscov is reused unchanged, and legitimately
    %  J.5  the conditional term carries no trial-averaging factor
    %  J.6  RESOLVED by S17-FIX: the per-participant error POOLS the
    %       mean-surface person x trial term exp(2*alpha_k(z))*K_k (per event,
    %       population-built) with the participant's own conditional variance,
    %       matching the group surface's sigma_pi_e2 construction and the
    %       owner's reference implementation. The gap against the plotted group
    %       curve remains two-signed (s_sd marginalization + own trial counts),
    %       and both directions stay pinned below.

    methods (Test)

        function testPerParticipantTableShapeOneDimension(testCase)
            %Contract check: one row per participant, in ssinfo row order, with
            %the (ll, pt, ul) triple on both coefficients.
            d = localDraws();
            ss = localSsinfo();
            [gen,dep] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',localErSs(1), 'er_ss2',localErSs(2), ...
                'idtable',ss, 'ndim',1, 'CI',0.95);

            testCase.verifyEqual(height(gen), height(ss));
            testCase.verifyEqual(height(dep), height(ss));
            testCase.verifyEqual(gen.id2(:)', ss.id2(:)', ...
                'row order / id2 must be preserved, not the remapped 1 used per call');
            for f = {'rel_pt','rel_ll','rel_ul'}
                testCase.verifyTrue(ismember(f{1}, gen.Properties.VariableNames));
            end
            testCase.verifyTrue(all(isfinite(gen.rel_pt)));
            testCase.verifyTrue(all(gen.rel_pt > 0 & gen.rel_pt < 1));
            %dependability uses absolute error, so it can never exceed
            %generalizability for the same participant
            testCase.verifyTrue(all(dep.rel_pt <= gen.rel_pt + 1e-12));
        end

        function testTwoDimensionsUsesTheInteractionColumn(testCase)
            %ndim = 2 takes KDIM = 6 and the design row [z1; z2; z1*z2]. A
            %converter that dropped the interaction term would still run, so the
            %check is that changing ONLY the interaction slope moves the answer.
            d = localDraws(6);
            ss = localSsinfo(true);
            [g0,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',localErSs(1), 'er_ss2',localErSs(2), ...
                'idtable',ss, 'ndim',2, 'CI',0.95);

            %columns 5/6 are the z1*z2 slopes for events 1/2
            bd = localGet(d,'b_dim'); bd(:,5) = bd(:,5) + 0.5;
            d2 = localOverride(d,'b_dim',bd);
            [g1,~] = psyrat_ssrel_diffdynrel_gamma_ls(d2{:}, ...
                'er_ss1',localErSs(1), 'er_ss2',localErSs(2), ...
                'idtable',ss, 'ndim',2, 'CI',0.95);

            %participants with z1*z2 == 0 are unmoved; at least one other must move
            interacts = abs(ss.z1 .* ss.z2) > 1e-12;
            testCase.verifyTrue(any(interacts), 'fixture must exercise z1*z2 ~= 0');
            testCase.verifyTrue(any(abs(g1.rel_pt(interacts) - ...
                g0.rel_pt(interacts)) > 1e-8), ...
                'the z1*z2 interaction slope must reach the coefficient');
        end

        function testResidualIsThePooledTermAtTheParticipantsOwnZ(testCase)
            %J.2/J.5/J.6 (S17-FIX). The per-participant residual variance per
            %event must be exactly
            %  exp(2*alpha_k(z_s))*K_k + exp(2*(er_ss_k + b_sigma_dim_k . X(z_s)))
            %- the population mean-surface person x trial term POOLED with the
            %participant's own conditional variance. No scaled-chi-square factor
            %2, and NO exp(2*sd_trl^2) trial-averaging factor on the conditional
            %term. The closed form writes K out from the section-J identity, so
            %this is an independent oracle, not a restatement of the converter.
            %
            %Checked through the coefficient rather than by reaching inside: solve
            %the reported generalizability for the relative error it implies, and
            %compare against the closed form built from the signal components.
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            [gen,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);

            [sigU, relerr] = localClosedForm(d, ss, e1, e2, 1);
            expected = sigU ./ (sigU + relerr);
            testCase.verifyEqual(gen.rel_pt, expected, 'AbsTol', 1e-9, ...
                ['per-participant G must equal sigma2_U / (sigma2_U + ' ...
                'sigma2_A,s/n1 + sigma2_B,s/n2) with each event''s residual ' ...
                'pooled: exp(2*alpha(z))*K + the participant''s own ' ...
                'conditional variance (S17-FIX)']);
        end

        function testResidualCarriesNoTrialAveragingFactor(testCase)
            %J.5, stated as a MUTATION check. Restoring exp(2*sd_trl^2) - the
            %log-nu path's trial-averaging factor - must change the answer, so a
            %future edit that "restores it for symmetry" cannot pass silently.
            %
            %Also pins the direction: inflating the residual lowers the
            %coefficient. Note it does NOT lower it by any constant factor, which
            %is why the shipped comment's old "by a constant exp(2*s_i^2)" wording
            %was wrong for the coefficient (it is right only for the residual).
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            [gen,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);

            %the mutant: add the trial-averaging factor to the stored log-SD,
            %i.e. multiply the residual VARIANCE by exp(2*sd_trl^2)
            sd_trl = localGet(d,'sd_trl');
            m1 = e1 + sd_trl(:,1).^2;
            m2 = e2 + sd_trl(:,2).^2;
            [mut,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',m1, 'er_ss2',m2, 'idtable',ss, 'ndim',1, 'CI',0.95);

            testCase.verifyTrue(all(abs(mut.rel_pt - gen.rel_pt) > 1e-6), ...
                'the trial-averaging factor must be observable in the output');
            testCase.verifyTrue(all(mut.rel_pt < gen.rel_pt), ...
                'inflating the residual must LOWER the coefficient');

            %and the mutant's effect is NOT a constant factor on the coefficient
            ratio = mut.rel_pt ./ gen.rel_pt;
            testCase.verifyTrue(max(ratio) - min(ratio) > 1e-6, ...
                ['the coefficient is a ratio, so it cannot move by a constant ' ...
                'factor - the corrected S17 wording depends on this']);
        end

        function testSignalIsNotConditionedOnTheParticipantsAmplitude(testCase)
            %J.3. sigma_p2 and sigma_i2 are BETWEEN-person variances evaluated at
            %the population alpha_k(z_s). The converter never receives a
            %per-participant log-MEAN effect at all, which is the structural
            %guarantee; this test pins the observable consequence, that two
            %participants sharing a z share a numerator and differ only through
            %their own residuals.
            d = localDraws();
            ss = localSsinfo();
            %participants 2 and 5 are constructed to share z1 (see localSsinfo)
            testCase.assertEqual(ss.z1(2), ss.z1(5), 'AbsTol', 0, ...
                'fixture must contain two participants at the same z');

            e1 = localErSs(1); e2 = localErSs(2);
            [sigU, ~] = localClosedForm(d, ss, e1, e2, 1);
            testCase.verifyEqual(sigU(2), sigU(5), 'AbsTol', 1e-12, ...
                'the signal must depend only on z, not on the participant');

            %give participant 5 a much larger residual: only their coefficient moves
            e1b = e1; e1b(:,5) = e1b(:,5) + 0.75;
            [g0,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);
            [g1,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1b, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);
            moved = abs(g1.rel_pt - g0.rel_pt) > 1e-9;
            testCase.verifyEqual(find(moved)', 5, ...
                'only the perturbed participant''s coefficient may move');
        end

        function testCrosscovArgumentsAreMeanSubmodelOnly(testCase)
            %J.4. Every cor_id entry touching the person SCALE slots (3-4) is
            %unread by this converter: the residual never touches mu, so those
            %correlations cannot enter. Only cor_id(1,2) - the cross-event MEAN
            %correlation - may move the answer.
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            [g0,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);

            cor = localGet(d,'cor_id');
            for pair = {[1 3],[2 4],[3 4],[1 4],[2 3]}
                p = pair{1};
                c2 = cor;
                c2(:,p(1),p(2)) = -0.9; c2(:,p(2),p(1)) = -0.9;
                d2 = localOverride(d,'cor_id',c2);
                [g1,~] = psyrat_ssrel_diffdynrel_gamma_ls(d2{:}, ...
                    'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);
                testCase.verifyEqual(g1.rel_pt, g0.rel_pt, 'AbsTol', 1e-12, ...
                    sprintf(['cor_id(%d,%d) involves a scale slot and must NOT ' ...
                    'enter the subject-level coefficient'], p(1), p(2)));
            end

            %positive control: cor_id(1,2) MUST move it, so the loop above cannot
            %pass by the converter ignoring cor_id wholesale
            c3 = cor; c3(:,1,2) = -0.4; c3(:,2,1) = -0.4;
            d3 = localOverride(d,'cor_id',c3);
            [g2,~] = psyrat_ssrel_diffdynrel_gamma_ls(d3{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);
            testCase.verifyTrue(any(abs(g2.rel_pt - g0.rel_pt) > 1e-8), ...
                'cor_id(1,2) is the difference numerator and must move the answer');
        end

        function testGapAgainstTheGroupCurveIsNotOneSigned(testCase)
            %Pinned in BOTH directions, which is the whole point - and STILL
            %two-signed after S17-FIX, because the pooling did not touch either
            %mechanism behind the flip.
            %
            %Both sides now carry the mean-surface person x trial term
            %exp(2*alpha(z))*K, and at the MATCHED trial counts this fixture
            %plants (localFixedZSsinfo mirrors the curve's obs) it drops out of
            %the error-term comparison exactly. What remains: the GROUP surface
            %marginalizes the conditional term over the person residual-SD
            %population (exp(2*log_sigma + 2*s_sd^2)) where a participant uses
            %their OWN realized w_s, and - in production plots, though
            %neutralized here - the curve divides by a reference n' where a
            %participant divides by their own counts (which also divides the
            %pooled term differently on the two sides; that remainder is
            %sign-aligned with the counts mechanism).
            %
            %An earlier version of this test asserted that per-participant values
            %sit ABOVE the curve, full stop. That is FALSE, and it passed only
            %because it planted w = 0 for everyone. The comparison reduces to the
            %error term, and a participant whose own scale effect exceeds s_sd^2
            %lands BELOW. Both halves are asserted here so the conditional nature
            %cannot quietly revert to an unconditional claim.
            d = localDraws();
            zfix = 0.4;
            b_sigma = localGet(d,'b_sigma');
            sd_id = localGet(d,'sd_id');
            s_sd = [sd_id(1,3) sd_id(1,4)];

            grp = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',zfix, ...
                'obs',[20 30], 'CI',0.95);

            %(1) a participant AT the population residual (w = 0) sits ABOVE
            ss = localFixedZSsinfo(zfix, 1);
            [above,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',b_sigma(:,1), 'er_ss2',b_sigma(:,2), ...
                'idtable',ss, 'ndim',1, 'CI',0.95);
            testCase.verifyGreaterThan(above.rel_pt, grp.G.pt, ...
                ['at w = 0 the per-participant coefficient must exceed the ' ...
                'group surface: both carry the pooled mean-surface term, but ' ...
                'the group residual additionally marginalizes the conditional ' ...
                'term over s_sd']);

            %(2) a participant with a LARGE own scale effect sits BELOW. w is
            %chosen well past the flip threshold rather than marginally past it,
            %so the assertion does not depend on the fixture's exact numbers.
            wbig = 3 * max(s_sd);
            [below,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',b_sigma(:,1) + wbig, 'er_ss2',b_sigma(:,2) + wbig, ...
                'idtable',ss, 'ndim',1, 'CI',0.95);
            testCase.verifyLessThan(below.rel_pt, grp.G.pt, ...
                ['a participant whose own residual-SD effect exceeds the ' ...
                'population value must fall BELOW the group curve - the gap ' ...
                'is not one-signed, and S17 must not be restated as if it were']);
        end

        function testFewerOwnTrialsAlsoFlipsTheGapAgainstTheCurve(testCase)
            %The SECOND mechanism behind the same finding, and the one most likely
            %to be forgotten: the group surface is evaluated at a reference n'
            %while each participant divides by their OWN per-event counts. Even at
            %an exactly average residual, a participant with materially fewer
            %trials than the reference falls below the curve.
            %
            %Pinned separately from the residual mechanism so that fixing one
            %cannot mask the other.
            d = localDraws();
            zfix = 0.4;
            b_sigma = localGet(d,'b_sigma');

            grp = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',zfix, ...
                'obs',[20 30], 'CI',0.95);

            ssFew = localFixedZSsinfo(zfix, 1);
            ssFew.trls1(1) = 4; ssFew.trls2(1) = 6;   % far below the reference n'
            [few,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',b_sigma(:,1), 'er_ss2',b_sigma(:,2), ...
                'idtable',ssFew, 'ndim',1, 'CI',0.95);

            testCase.verifyLessThan(few.rel_pt, grp.G.pt, ...
                ['with w = 0 but far fewer trials than the reference n'', the ' ...
                'participant must fall BELOW the curve - trial counts flip the ' ...
                'comparison independently of the residual (S17)']);
        end

        function testUnequalPerEventTrialCountsAreHonoured(testCase)
            %The difference design is routinely unbalanced and each participant
            %brings their own two counts. Swapping one participant's counts must
            %move only that participant.
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            [g0,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss, 'ndim',1, 'CI',0.95);

            ss2 = ss; ss2.trls1(3) = ss.trls1(3) * 3;
            [g1,~] = psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1, 'er_ss2',e2, 'idtable',ss2, 'ndim',1, 'CI',0.95);

            moved = abs(g1.rel_pt - g0.rel_pt) > 1e-9;
            testCase.verifyEqual(find(moved)', 3);
            testCase.verifyTrue(g1.rel_pt(3) > g0.rel_pt(3), ...
                'more trials on one event must raise that participant''s G');
        end

        function testShortErSsMatrixFailsLoudly(testCase)
            %psyrat_ssrel_diff's own range check falls back to the POPULATION
            %residual for an unindexable id2, silently. The converter must catch
            %that first, or a truncated draw matrix would turn a subject-level
            %table into a population one with no warning.
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            testCase.verifyError(@() psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1(:,1:2), 'er_ss2',e2, 'idtable',ss, 'ndim',1, ...
                'CI',0.95), 'varargin:er_ss');
        end

        function testInputValidation(testCase)
            d = localDraws();
            ss = localSsinfo();
            e1 = localErSs(1); e2 = localErSs(2);
            common = {'er_ss1',e1,'er_ss2',e2,'idtable',ss,'CI',0.95};

            testCase.verifyError(@() psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                common{:}, 'ndim',3), 'varargin:ndim');

            %sd_id must be the 4-column person block
            sdid = localGet(d,'sd_id');
            dbad = localOverride(d,'sd_id', sdid(:,1:3));
            testCase.verifyError(@() psyrat_ssrel_diffdynrel_gamma_ls(dbad{:}, ...
                common{:}, 'ndim',1), 'varargin:sd_id');

            %KDIM must match ndim
            testCase.verifyError(@() psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                common{:}, 'ndim',2), 'varargin:b_dim');

            %a missing required argument is loud
            testCase.verifyError(@() psyrat_ssrel_diffdynrel_gamma_ls(d{:}, ...
                'er_ss1',e1,'er_ss2',e2,'ndim',1,'CI',0.95), 'varargin:idtable');
        end

    end
end

% -------------------------------------------------------------------------

function args = localDraws(kdim)
%Deterministic synthetic posterior draws in the case-12/13 location-scale layout.
%Fixed values, not random: these tests pin arithmetic, so a seed would add
%nothing and would make failures harder to read. Kept deliberately identical to
%TestGammaDiffDynrelLsConversion's fixture so the two increments' numbers are
%directly comparable.
if nargin < 1, kdim = 2; end
nd = 200;
o = ones(nd,1);

b       = [1.75*o, 1.60*o];          % per-event log-mean intercepts
b_sigma = [-0.25*o, -0.10*o];        % per-event log residual-SD intercepts
if kdim == 2
    b_dim       = [0.20*o, 0.05*o];
    b_sigma_dim = [0.15*o, -0.05*o];
else
    b_dim       = [0.20*o, 0.05*o, 0.10*o, 0.02*o, 0.04*o, 0.01*o];
    b_sigma_dim = [0.15*o, -0.05*o, 0.06*o, 0.03*o, 0.02*o, 0.01*o];
end
sd_id  = [0.32*o, 0.29*o, 0.24*o, 0.21*o];  % [m1 m2 logsig1 logsig2]
sd_trl = [0.11*o, 0.13*o];

R = eye(4);
R(1,2) = 0.45; R(2,1) = 0.45;
R(1,3) = 0.30; R(3,1) = 0.30;
R(2,4) = 0.25; R(4,2) = 0.25;
R(3,4) = 0.35; R(4,3) = 0.35;
cor_id = repmat(reshape(R, [1 4 4]), [nd 1 1]);

args = {'b',b, 'b_sigma',b_sigma, 'b_dim',b_dim, 'b_sigma_dim',b_sigma_dim, ...
    'sd_id',sd_id, 'sd_trl',sd_trl, 'cor_id',cor_id, 'cor_i',0.20*o};
end

function ss = localSsinfo(twodim)
%Five participants with unequal per-event trial counts. Participants 2 and 5
%share z1 deliberately, so a test can check that the signal depends on z alone.
if nargin < 1, twodim = false; end
id2 = (1:5)';
z1  = [-1.2; 0.4; 0.0; 1.1; 0.4];
trls1 = [18; 22; 15; 30; 25];
trls2 = [27; 31; 20; 41; 33];
trls  = round((trls1 + trls2)/2);
if twodim
    z2 = [0.7; -0.5; 0.9; 0.3; -0.5];
    ss = table(id2, id2, z1, z2, trls1, trls2, trls, 'VariableNames', ...
        {'id','id2','z1','z2','trls1','trls2','trls'});
else
    ss = table(id2, id2, z1, trls1, trls2, trls, 'VariableNames', ...
        {'id','id2','z1','trls1','trls2','trls'});
end
end

function ss = localFixedZSsinfo(zfix, n)
%n participants all sitting at the same z, with the trial counts the group
%surface is evaluated at (obs = [20 30]) so the curve comparison isolates the
%residual rather than confounding it with replication.
id2 = (1:n)';
ss = table(id2, id2, repmat(zfix,n,1), repmat(20,n,1), repmat(30,n,1), ...
    repmat(25,n,1), 'VariableNames', ...
    {'id','id2','z1','trls1','trls2','trls'});
end

function e = localErSs(evt)
%Per-subject log residual SD at z = 0 for one event: the population intercept
%plus a fixed, distinct person offset per participant.
nd = 200;
w = [-0.30, 0.10, 0.00, 0.25, -0.15];   % person scale effects
base = [-0.25, -0.10];                   % b_sigma per event
e = repmat(base(evt) + w, nd, 1);
end

function [sigU, relerr] = localClosedForm(args, ss, e1, e2, ~)
%Independent closed form for the one-dimension case, written from the section-J
%formulas rather than from the converter, so it is a cross-check and not a copy
%of the implementation's control flow.
b = localGet(args,'b'); b_sigma = localGet(args,'b_sigma');
b_dim = localGet(args,'b_dim'); b_sigma_dim = localGet(args,'b_sigma_dim');
sd_id = localGet(args,'sd_id'); sd_trl = localGet(args,'sd_trl');
cor_id = localGet(args,'cor_id'); cor_i = localGet(args,'cor_i');

%the fixture is constant across draws, so take draw 1 and work in scalars
sp = [sd_id(1,1) sd_id(1,2)]; si = [sd_trl(1,1) sd_trl(1,2)];
cp = cor_id(1,1,2); ci = cor_i(1);

n = height(ss);
sigU = zeros(n,1); relerr = zeros(n,1);
for k = 1:n
    z = ss.z1(k);
    aA = b(1,1) + b_dim(1,1)*z;   aB = b(1,2) + b_dim(1,2)*z;
    %population person components (NOT at the participant's own amplitude)
    sp2A = exp(2*aA + si(1)^2)*exp(sp(1)^2)*(exp(sp(1)^2)-1);
    sp2B = exp(2*aB + si(2)^2)*exp(sp(2)^2)*(exp(sp(2)^2)-1);
    EYA = exp(aA + 0.5*(sp(1)^2 + si(1)^2));
    EYB = exp(aB + 0.5*(sp(2)^2 + si(2)^2));
    covP = EYA*EYB*(exp(cp*sp(1)*sp(2)) - 1);
    sigU(k) = sp2A + sp2B - 2*covP;
    %the participant's own conditional residual at their own z
    sA = exp(2*(e1(1,k) + b_sigma_dim(1,1)*z));
    sB = exp(2*(e2(1,k) + b_sigma_dim(1,2)*z));
    %the pooled mean-surface person x trial term per event (S17-FIX), written
    %from the section-J identity K = exp(sp^2+si^2)(exp(sp^2)-1)(exp(si^2)-1)
    KA = exp(sp(1)^2 + si(1)^2)*(exp(sp(1)^2)-1)*(exp(si(1)^2)-1);
    KB = exp(sp(2)^2 + si(2)^2)*(exp(sp(2)^2)-1)*(exp(si(2)^2)-1);
    vA = exp(2*aA)*KA;
    vB = exp(2*aB)*KB;
    relerr(k) = (vA + sA)/ss.trls1(k) + (vB + sB)/ss.trls2(k);
end
%b_sigma and ci are unused on this path by construction (the relative error has
%no trial-covariance term); reference them so a reader sees that is deliberate.
assert(~isempty(b_sigma) && ~isempty(ci));
end

function v = localGet(args, name)
idx = find(strcmp(args, name), 1);
v = args{idx + 1};
end

function args = localOverride(args, name, value)
idx = find(strcmp(args, name), 1);
args{idx + 1} = value;
end
