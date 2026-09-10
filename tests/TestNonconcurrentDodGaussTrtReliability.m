classdef TestNonconcurrentDodGaussTrtReliability < matlab.unittest.TestCase
    %TESTNONCONCURRENTDODGAUSSTRTRELIABILITY Estimand tests for the GAUSSIAN
    % arm of the TWO-FACET person-specific dynamic NONCONCURRENT
    % difference-of-differences (analysis 29, ic_dodiff_dynrel_sserrvar_trt,
    % family = 'gaussian'). The one-facet twin (analysis 28), the shared
    % disclosure surfaces, and the vetted-kernel tie live in
    % TestNonconcurrentDodGaussReliability; this class pins what the second
    % facet adds: the seven-channel closed forms, the CES/CE/CS decomposition
    % mapping onto reltype 3/1/2, the occasion-n' selector semantics, and the
    % two-facet floor disclosure.
    %
    % Same evidence posture as the one-facet class: no frozen fixture lane
    % exists for this arm (the reference bundle is gamma end to end); live
    % CmdStan planted-truth recovery is the load-bearing end-to-end gate
    % (SCIENTIFIC_FORMULA_AUDIT.md section 26, owner ruling 3 - extended to
    % the two-facet analysis exactly as the gamma-29 ruling was).
    %
    % Tolerances: AbsTol 1e-10 for longhand closed-form transcriptions;
    % AbsTol 1e-12 for same-arithmetic routes.

    methods (Test)

        function testPerParticipantCesMatchesClosedForm(testCase)
            %Constant draws, reltype 3 (CES): every channel transcribed
            %longhand. Under the identity link pio_e is diagonal, so its
            %harmonic composite is exactly sum_q c_q^2*sigma_q^2/(ni_q*no_q).
            fx = localGaussTrtConstFixture();
            ssrel = localRunGaussTrtSsrel(fx, 3);
            w = fx.cvec(:);
            ch = localChannels(fx);

            for c = 1:fx.P
                lh = localLonghand(fx, ch, w, c);
                testCase.verifyEqual(ssrel.uni_pt(c), lh.u, 'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.relerr_pt(c), lh.rel_ces, ...
                    'AbsTol', 1e-10, sprintf('CES relerr row %d', c));
                testCase.verifyEqual(ssrel.abserr_pt(c), lh.abs_ces, ...
                    'AbsTol', 1e-10, sprintf('CES abserr row %d', c));
                testCase.verifyEqual(ssrel.gen_pt(c), ...
                    lh.u / (lh.u + lh.rel_ces), 'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.dep_pt(c), ...
                    lh.u / (lh.u + lh.abs_ces), 'AbsTol', 1e-10);
                testCase.verifyEqual(ssrel.residual_var_dod_person_pt(c), ...
                    lh.resvar, 'AbsTol', 1e-10, ...
                    'the four-term residual identity must hold');
                %constant draws collapse the intervals
                testCase.verifyEqual(ssrel.gen_ll(c), ssrel.gen_pt(c), ...
                    'AbsTol', 1e-12);
                testCase.verifyEqual(ssrel.invalid_frac(c), 0);
            end
            %the identity-link expected score carries no facet-variance
            %correction: expected_dod_person is the additive predictor's
            %contrast, with no +0.5*v_other lognormal term
            c = 1;
            mu_own = zeros(1,4);
            for q = 1:4
                mu_own(q) = fx.mu_ss.(sprintf('c%d',q))(1,c) + ...
                    fx.b_dim(1,q,1) * fx.idtable.z1(c);
            end
            testCase.verifyEqual(ssrel.expected_dod_person_pt(1), ...
                mu_own * w, 'AbsTol', 1e-10);
        end

        function testReltypeSelectsTheMatchingDecomposition(testCase)
            %CE moves the scaled person x occasion term into the universe;
            %CS moves the scaled person x trial term; CES is pure person
            %variance. Pinned longhand per reltype, and the three tables must
            %genuinely differ.
            fx = localGaussTrtConstFixture();
            w = fx.cvec(:);
            ch = localChannels(fx);
            tCE = localRunGaussTrtSsrel(fx, 1);
            tCS = localRunGaussTrtSsrel(fx, 2);
            tCES = localRunGaussTrtSsrel(fx, 3);

            for c = 1:fx.P
                lh = localLonghand(fx, ch, w, c);
                testCase.verifyEqual(tCE.uni_pt(c), lh.u + lh.po_h, ...
                    'AbsTol', 1e-10, 'CE universe must add the scaled po term');
                testCase.verifyEqual(tCE.relerr_pt(c), lh.pi_h + lh.pio_h, ...
                    'AbsTol', 1e-10);
                testCase.verifyEqual(tCE.abserr_pt(c), ...
                    lh.pi_h + lh.pio_h + lh.i_h + lh.io_h, 'AbsTol', 1e-10);
                testCase.verifyEqual(tCS.uni_pt(c), lh.u + lh.pi_h, ...
                    'AbsTol', 1e-10, 'CS universe must add the scaled pi term');
                testCase.verifyEqual(tCS.relerr_pt(c), lh.po_h + lh.pio_h, ...
                    'AbsTol', 1e-10);
                testCase.verifyEqual(tCS.abserr_pt(c), ...
                    lh.po_h + lh.pio_h + lh.o_h + lh.io_h, 'AbsTol', 1e-10);
                testCase.verifyEqual(tCES.uni_pt(c), lh.u, 'AbsTol', 1e-10, ...
                    'CES universe must be pure person variance');
            end
            testCase.verifyNotEqual(tCE.gen_pt(1), tCES.gen_pt(1));
            testCase.verifyNotEqual(tCS.gen_pt(1), tCES.gen_pt(1));
            testCase.verifyEqual(tCE.reltype(1), 1);
            testCase.verifyEqual(tCS.reltype(1), 2);
            testCase.verifyEqual(tCES.reltype(1), 3);
        end

        function testTwoFacetFloorsAreSeparate(testCase)
            %B16, both facets: a zero trial cell and a zero occasion cell are
            %floored to 1 for the composite, disclosed in SEPARATE columns,
            %and the true zeros stay in the table.
            fx = localGaussTrtConstFixture();
            fx.idtable.trls2(2) = 0;
            fx.idtable.occs3(2) = 0;
            ssrel = localRunGaussTrtSsrel(fx, 3);
            testCase.verifyEqual(ssrel.trls2(2), 0);
            testCase.verifyEqual(ssrel.occs3(2), 0);
            testCase.verifyEqual(ssrel.n_floored(2), 1);
            testCase.verifyEqual(ssrel.n_floored_occ(2), 1);
            testCase.verifyEqual(ssrel.n_floored(1), 0);
            testCase.verifyEqual(ssrel.n_floored_occ(1), 0);
            testCase.verifyTrue(isfinite(ssrel.gen_pt(2)));

            %corrupt counts still error loudly on either facet
            fxbad = localGaussTrtConstFixture();
            fxbad.idtable.occs1(1) = NaN;
            testCase.verifyError(@() localRunGaussTrtSsrel(fxbad, 3), ...
                'varargin:n');
        end

        function testSummaryDispatchGaussianTrt(testCase)
            %The summary routes a Gaussian analysis-29 REL through the
            %Gaussian trt calculators; the occasion-n' selector keeps its
            %semantics (default 1; 'observed' = the stored per-cell medians).
            REL = localSyntheticGaussTrtRel(localGaussTrtConstFixture());
            w = warning('off', 'dynrel:dodinvalidity');
            cleaner = onCleanup(@() warning(w)); %#ok<NASGU>

            %default: occasion n' = 1 per cell
            summ1 = psyrat_dynrel_summary(REL);
            st = summ1.strata(1);
            testCase.assertTrue(istable(st.ssrel_table));
            testCase.verifyEqual(st.nocc, 1);
            testCase.verifyEqual(numel(st.obs_ev), 4);
            testCase.verifyEqual(summ1.dod_labels.inference_framework, ...
                'four_gaussian_margins_all_residual_covariances_fixed_zero_v1');

            %'observed' opts into the stored per-cell medians
            summ2 = psyrat_dynrel_summary(REL, 'nocc', 'observed');
            testCase.verifyEqual(summ2.strata(1).nocc, ...
                round(mean(max(1, round(REL.out.nocc(:,1)')))));
            %the two surfaces must differ: occasion terms scale with n'_o
            testCase.verifyNotEqual(summ1.strata(1).G.pt(1), ...
                summ2.strata(1).G.pt(1));

            %reltype reaches the table
            summ3 = psyrat_dynrel_summary(REL, 'reltype', 1);
            testCase.verifyEqual(summ3.strata(1).ssrel_table.reltype(1), 1);

            %fail-loud arms hold for the trt string too
            RELlognu = REL;
            RELlognu.family = 'gamma';
            RELlognu.gammascale = 1;
            testCase.verifyError(@() psyrat_dynrel_summary(RELlognu), ...
                'dynrel:dodgammals');
        end

    end
end

% ---------------------------------------------------------------------------
% Fixtures: constant draws on Gaussian (identity-link) scales, all-distinct
% per-cell and per-facet values so channel swaps fail rather than passing by
% symmetry. Seeds/values distinct from the one-facet class.
% ---------------------------------------------------------------------------

function fx = localGaussTrtConstFixture()
nd = 4;
fx.nd = nd;
fx.P = 3;
fx.cvec = [1 -1 -1 1];
fx.ndim = 1;
fx.b = repmat([-4.8 -1.9 -4.1 -3.4], nd, 1);
fx.b_sigma = repmat([0.61 1.01 0.64 0.66], nd, 1);
fx.b_dim = repmat(reshape([0.6 -0.25 0.4 -0.1], [1 4 1]), nd, 1);
fx.b_sigma_dim = repmat(reshape([-0.08 0.12 -0.05 0.07], [1 4 1]), nd, 1);
fx.sd_id = repmat([2.0 1.7 2.2 1.9 0.30 0.34 0.28 0.32], nd, 1);
fx.cor_id = repmat(reshape(localCorr8(), [1 8 8]), nd, 1);
fx.sd_trl = repmat([0.80 0.70 0.90 0.75], nd, 1);
fx.cor_trl = repmat(reshape([1 .25 .17 -.14; .25 1 -.19 .21; ...
    .17 -.19 1 .16; -.14 .21 .16 1], [1 4 4]), nd, 1);
fx.sd_occ = repmat([0.45 0.38 0.52 0.41], nd, 1);
fx.cor_occ = repmat(reshape([1 -.22 .19 .13; -.22 1 .15 .20; ...
    .19 .15 1 -.17; .13 .20 -.17 1], [1 4 4]), nd, 1);
fx.sd_tid = repmat([0.40 0.47 0.35 0.44], nd, 1);
fx.cor_tid = repmat(reshape([1 .18 -.16 .23; .18 1 .21 -.11; ...
    -.16 .21 1 .14; .23 -.11 .14 1], [1 4 4]), nd, 1);
fx.sd_oid = repmat([0.34 0.28 0.39 0.31], nd, 1);
fx.cor_oid = repmat(reshape([1 .24 .12 -.18; .24 1 -.13 .16; ...
    .12 -.13 1 .25; -.18 .16 .25 1], [1 4 4]), nd, 1);
fx.sd_to = repmat([0.30 0.38 0.26 0.34], nd, 1);
fx.cor_to = repmat(reshape([1 -.15 .22 .11; -.15 1 .17 .19; ...
    .22 .17 1 -.21; .11 .19 -.21 1], [1 4 4]), nd, 1);

uloc = [0.9 -0.5 0.3 -0.7; -1.1 0.6 -0.2 0.8; 0.2 0.1 -0.4 0.5];
ulsig = [0.15 -0.10 0.05 -0.20; -0.12 0.18 -0.06 0.09; 0.05 -0.03 0.11 -0.14];
for q = 1:4
    fx.mu_ss.(sprintf('c%d',q)) = repmat(fx.b(1,q) + uloc(:,q)', nd, 1);
    fx.er_ss.(sprintf('c%d',q)) = repmat(fx.b_sigma(1,q) + ulsig(:,q)', nd, 1);
end

fx.idtable = table({'s1';'s2';'s3'}, (1:3)', [-0.8; 0.1; 1.2], ...
    [12; 8; 20], [9; 15; 6], [14; 10; 18], [11; 13; 7], ...
    [2; 3; 2], [3; 2; 2], [2; 2; 3], [3; 3; 2], ...
    'VariableNames', {'id','id2','z1','trls1','trls2','trls3','trls4', ...
    'occs1','occs2','occs3','occs4'});
fx.nprime = [30 30 30 30];
fx.nprime_o = [2 2 2 2];
end

function R = localCorr8()
L = [0.55 0.30; 0.48 -0.35; -0.42 0.28; 0.36 0.22; ...
    0.25 0.31; -0.28 0.24; 0.22 -0.26; 0.30 0.20];
S = L * L' + diag([0.60 0.68 0.75 0.83 0.71 0.77 0.82 0.74]);
R = diag(1 ./ sqrt(diag(S))) * S * diag(1 ./ sqrt(diag(S)));
end

function ch = localChannels(fx)
%The six signal channel covariance matrices from the first (constant) draw.
cov4 = @(sd, cor) diag(sd(1,:)) * squeeze(cor(1,:,:)) * diag(sd(1,:));
ch.p = cov4(fx.sd_id(:,1:4), fx.cor_id(:,1:4,1:4));
ch.i = cov4(fx.sd_trl, fx.cor_trl);
ch.o = cov4(fx.sd_occ, fx.cor_occ);
ch.pi = cov4(fx.sd_tid, fx.cor_tid);
ch.po = cov4(fx.sd_oid, fx.cor_oid);
ch.io = cov4(fx.sd_to, fx.cor_to);
end

function lh = localLonghand(fx, ch, w, c)
%Longhand per-participant channels at the participant's own z and ACTUAL
%counts: trial terms / n_i, occasion terms / n_o, product terms / n_i*n_o;
%the diagonal residual reduces to the plain weighted sum.
z = fx.idtable.z1(c);
ni = [fx.idtable.trls1(c) fx.idtable.trls2(c) ...
    fx.idtable.trls3(c) fx.idtable.trls4(c)];
no = [fx.idtable.occs1(c) fx.idtable.occs2(c) ...
    fx.idtable.occs3(c) fx.idtable.occs4(c)];
sig = zeros(1,4);
for q = 1:4
    sig(q) = exp(fx.er_ss.(sprintf('c%d',q))(1,c) + ...
        fx.b_sigma_dim(1,q,1) * z);
end
lh.u = w' * ch.p * w;
lh.i_h = localHarmonicLonghand(ch.i, w, ni);
lh.o_h = localHarmonicLonghand(ch.o, w, no);
lh.pi_h = localHarmonicLonghand(ch.pi, w, ni);
lh.po_h = localHarmonicLonghand(ch.po, w, no);
lh.io_h = localHarmonicLonghand(ch.io, w, ni .* no);
lh.pio_h = sum((w'.^2) .* sig.^2 ./ (ni .* no));
lh.rel_ces = lh.pi_h + lh.po_h + lh.pio_h;
lh.abs_ces = lh.rel_ces + lh.i_h + lh.o_h + lh.io_h;
lh.resvar = sum(sig.^2);
end

function v = localHarmonicLonghand(S, w, n)
v = 0;
for q = 1:4
    v = v + w(q)^2 * S(q,q) / n(q);
end
for a = 1:3
    for b = (a+1):4
        h = 2 / (1/n(a) + 1/n(b));
        v = v + 2 * w(a) * w(b) * S(a,b) / h;
    end
end
end

function ssrel = localRunGaussTrtSsrel(fx, reltype)
ssrel = psyrat_ssrel_dodiffdynrel_trt('b', fx.b, 'b_sigma', fx.b_sigma, ...
    'b_dim', fx.b_dim, 'b_sigma_dim', fx.b_sigma_dim, ...
    'sd_id', fx.sd_id, 'cor_id', fx.cor_id, ...
    'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
    'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
    'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
    'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
    'sd_to', fx.sd_to, 'cor_to', fx.cor_to, ...
    'mu_ss1', fx.mu_ss.c1, 'mu_ss2', fx.mu_ss.c2, ...
    'mu_ss3', fx.mu_ss.c3, 'mu_ss4', fx.mu_ss.c4, ...
    'er_ss1', fx.er_ss.c1, 'er_ss2', fx.er_ss.c2, ...
    'er_ss3', fx.er_ss.c3, 'er_ss4', fx.er_ss.c4, ...
    'idtable', fx.idtable, 'nprime', fx.nprime, ...
    'nprime_o', fx.nprime_o, 'reltype', reltype, ...
    'ndim', fx.ndim, 'CI', 0.95, 'cvec', fx.cvec);
end

function REL = localSyntheticGaussTrtRel(fx)
%A synthetic Gaussian analysis-29 REL in the engine's storage idiom. The
%gammascale stamp is 1 deliberately (the engine stamps it on every run).
REL = struct();
REL.analysis = 'ic_dodiff_dynrel_sserrvar_trt';
REL.family = 'gaussian';
REL.gammascale = 1;
REL.ndim = 1;
REL.dim_names = {'dim1'};
out = struct();
out.labels = {'none'};
out.b = {num2cell(fx.b)};
out.b_sigma = {num2cell(fx.b_sigma)};
out.b_dim = {num2cell(fx.b_dim)};
out.b_sigma_dim = {num2cell(fx.b_sigma_dim)};
out.sd_id = {num2cell(fx.sd_id)};
out.sd_trl = {num2cell(fx.sd_trl)};
out.cor_id = {num2cell(fx.cor_id)};
out.cor_trl = {num2cell(fx.cor_trl)};
out.sd_occ = {num2cell(fx.sd_occ)};
out.cor_occ = {num2cell(fx.cor_occ)};
out.sd_tid = {num2cell(fx.sd_tid)};
out.cor_tid = {num2cell(fx.cor_tid)};
out.sd_oid = {num2cell(fx.sd_oid)};
out.cor_oid = {num2cell(fx.cor_oid)};
out.sd_to = {num2cell(fx.sd_to)};
out.cor_to = {num2cell(fx.cor_to)};
for q = 1:4
    out.(sprintf('mu_ss%d',q)) = {num2cell(fx.mu_ss.(sprintf('c%d',q)))};
    out.(sprintf('er_var_ss%d',q)) = {num2cell(fx.er_ss.(sprintf('c%d',q)))};
end
out.ssinfo = {fx.idtable};
out.dimz = {fx.idtable.z1};
out.ntrials = [13; 10; 14; 11];
out.nocc = [2; 2; 2; 3];
out.elabels = {'E1';'E2';'E3';'E4'};
out.glabels = {'none'};
out.dod_labels = struct(...
    'inference_framework',...
    'four_gaussian_margins_all_residual_covariances_fixed_zero_v1',...
    'cut_posterior_status',...
    'no_dependence_module_exists_because_all_constituents_are_assumed_nonconcurrent',...
    'residual_dependence_structure',...
    'all six pairwise conditional residual covariances fixed exactly to zero',...
    'scale_parameterization','person_specific_log_residual_sd',...
    'estimand_version','signed_four_component_nonconcurrent_dod_v1');
REL.out = out;
end
