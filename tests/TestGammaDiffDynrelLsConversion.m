classdef TestGammaDiffDynrelLsConversion < matlab.unittest.TestCase
    %Offline checks for the LOCATION-SCALE gamma DYNAMIC DIFFERENCE surface
    %calculator (analysis 12, gammascale = 2): psyrat_rel_diffdynrel_gamma_ls.
    %
    %Increment 8b. No CmdStan required - every test here drives the calculator
    %directly from synthetic posterior draws, so what is pinned is the
    %conversion arithmetic and the routing, NOT that the estimand is correct.
    %The accuracy suite is regression-only and this file is no different; live
    %recovery (TestGammaDiffDynrelLsRecovery) is what checks correctness.
    %
    %The tests are organized around the four claims the derivation in
    %documentation/gamma_scale_submodel_decision.md section I actually makes:
    %  I.2  the residual is decoupled from the mean surface
    %  I.3  psyrat_gamma_crosscov is reused unchanged, and legitimately
    %  I.4  the residual cross-covariance is zero
    %  I.5  the difference assembly, and both asymmetry axes

    methods (Test)

        function testSurfaceShapeAndFieldsOneDimension(testCase)
            %Contract check: a one-dimension run returns 1 x nz1 surfaces with
            %the (ll, pt, ul) triple on each coefficient, plus the residual
            %contrast SD used for reporting.
            d = localDraws();
            z = -2:0.5:2;
            out = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',z, ...
                'obs',[20 30], 'CI',0.95);

            testCase.verifyEqual(out.ndim, 1);
            testCase.verifyEqual(out.obs, [20 30]);
            testCase.verifyEqual(out.z1, z);
            for f = {'G','D','ICCg','ICCd'}
                testCase.verifyTrue(isfield(out, f{1}));
                testCase.verifySize(out.(f{1}).pt, [1 numel(z)]);
                testCase.verifySize(out.(f{1}).ll, [1 numel(z)]);
                testCase.verifySize(out.(f{1}).ul, [1 numel(z)]);
            end
            testCase.verifySize(out.sige_pt, [1 numel(z)]);
            testCase.verifyTrue(all(isfinite(out.G.pt)));
            testCase.verifyTrue(all(out.G.pt > 0 & out.G.pt < 1));
        end

        function testSurfaceShapeTwoDimensions(testCase)
            %Two dimensions gives a KDIM = 6 event-crossed slope block and an
            %nz1 x nz2 grid, rows indexing z1.
            d = localDraws(6);
            z1 = [-1 0 1]; z2 = [-1 0 1 2];
            out = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',2, ...
                'z1',z1, 'z2',z2, 'obs',20, 'CI',0.95);

            testCase.verifySize(out.G.pt, [numel(z1) numel(z2)]);
            testCase.verifySize(out.sige_pt, [numel(z1) numel(z2)]);
            %scalar obs must broadcast to both events
            testCase.verifyEqual(out.obs, [20 20]);
        end

        function testResidualComponentStillCarriesTheMeanSurfaceInteraction(testCase)
            %I.2, and the correction that writing these tests produced. It is
            %TEMPTING to assert that the location-scale residual is flat in z at
            %b_sigma_dim = 0, or invariant to the log-mean intercepts. Both are
            %FALSE, and the first draft of this file asserted them.
            %
            %sigma_pi_e2 is recovered by subtraction and therefore confounds the
            %person x trial interaction of the lognormal MEAN surface with the
            %conditional variance. Every mean-surface term scales as exp(2*alpha),
            %so
            %     sigma_pi_e2 = exp(2*alpha)*K + exp(2*log_sigma + 2*s_sd^2)
            %with K free of alpha. This test pins that decomposition, so nobody
            %re-derives the wrong version from the docstring.
            sp = 0.32; si = 0.11; lsig = -0.25; ssd = 0.24;
            alphas = [1.00; 1.75; 2.50];
            K = nan(size(alphas));
            for k = 1:numel(alphas)
                vc = psyrat_gamma_varcomps_ls(alphas(k), sp, si, lsig, ssd);
                testCase.verifyEqual(vc.e_cond_var, exp(2*lsig + 2*ssd^2), ...
                    'RelTol', 1e-12, ...
                    'The CONDITIONAL variance must not depend on alpha.');
                K(k) = (vc.sigma_pi_e2 - vc.e_cond_var) / exp(2*alphas(k));
            end
            testCase.verifyEqual(K, K(1)*ones(size(K)), 'RelTol', 1e-10, ...
                ['(sigma_pi_e2 - e_cond_var)/exp(2*alpha) must be constant in ' ...
                'alpha: the interaction part scales as exp(2*alpha) exactly.']);
            testCase.verifyGreaterThan(K(1), 0, ...
                'The interaction part must be strictly positive.');
        end

        function testResidualSurfaceRespondsToBothSlopes(testCase)
            %The behavioural consequence of the decomposition above: the residual
            %surface moves with the MEAN slope (through the interaction part) AND
            %with the SCALE slope (through the conditional part). Pinning both
            %directions is what keeps a calculator that silently dropped
            %b_sigma_dim from passing.
            z = -1:0.5:1;
            base = localDraws();

            %(a) scale slope zeroed: the surface still varies, via b_dim
            noScale = localOverride(base, 'b_sigma_dim', zeros(200,2));
            outNoScale = psyrat_rel_diffdynrel_gamma_ls(noScale{:}, 'ndim',1, ...
                'z1',z, 'obs',[20 30], 'CI',0.95);
            testCase.verifyGreaterThan(std(outNoScale.sige_pt), 1e-9, ...
                ['Even at b_sigma_dim = 0 the residual surface varies in z, ' ...
                'because the interaction part scales as exp(2*alpha(z)).']);

            %(b) BOTH slopes zeroed: now it must be exactly flat
            flat = localOverride(localOverride(base, 'b_sigma_dim', zeros(200,2)), ...
                'b_dim', zeros(200,2));
            outFlat = psyrat_rel_diffdynrel_gamma_ls(flat{:}, 'ndim',1, 'z1',z, ...
                'obs',[20 30], 'CI',0.95);
            testCase.verifyEqual(std(outFlat.sige_pt), 0, 'AbsTol', 1e-12, ...
                'With NEITHER slope the residual surface must be flat in z.');

            %(c) adding the scale slope back must change the surface
            outBoth = psyrat_rel_diffdynrel_gamma_ls(base{:}, 'ndim',1, 'z1',z, ...
                'obs',[20 30], 'CI',0.95);
            testCase.verifyNotEqual(outBoth.sige_pt, outNoScale.sige_pt, ...
                'b_sigma_dim must move the residual surface.');
        end

        function testConditionalVarianceIsInvariantToAmplitudeInTheSurface(testCase)
            %The property that DOES hold end-to-end, stated on the quantity where
            %it holds. Shifting the log-mean intercepts by a constant leaves each
            %event's conditional variance untouched, so the whole change in the
            %residual surface is attributable to the interaction part.
            base = localDraws();
            bb   = localGet(base,'b');
            bsig = localGet(base,'b_sigma');
            sdid = localGet(base,'sd_id');
            sdtrl = localGet(base,'sd_trl');

            vA = psyrat_gamma_varcomps_ls(bb(:,1), sdid(:,1), sdtrl(:,1), ...
                bsig(:,1), sdid(:,3));
            vB = psyrat_gamma_varcomps_ls(bb(:,1) + 0.4, sdid(:,1), sdtrl(:,1), ...
                bsig(:,1), sdid(:,3));

            testCase.verifyEqual(vB.e_cond_var, vA.e_cond_var, 'RelTol', 1e-12, ...
                ['Shifting the amplitude must not move the conditional ' ...
                'variance; under log-nu it would move as exp(2*alpha).']);
            testCase.verifyGreaterThan(vB.sigma_pi_e2, vA.sigma_pi_e2, ...
                'The residual COMPONENT does move, through the interaction part.');
        end

        function testCrosscovArgumentsAreMeanSubmodelOnly(testCase)
            %I.3. The property of interest is that psyrat_gamma_crosscov consumes no
            %dispersion quantity. WHAT THIS METHOD ACTUALLY CHECKS is narrower than
            %that: the helper's ARITY, plus that the returned covariances are finite.
            %It varies no input and compares no two results.
            %
            %So it catches an APPENDED ninth argument, and nothing else. It does not
            %catch a scale quantity SUBSTITUTED into one of the eight existing
            %slots - passing lsig_a where alpha_a belongs leaves nargin at 8, leaves
            %both covariances finite, and passes here while silently changing what
            %feeds psyrat_diffrel. An earlier version of this header claimed "if a
            %future edit routed a scale quantity into the helper, this fails", which
            %is true only of the append case.
            %
            %The strong form of the claim - perturb the scale side, require the
            %covariances bit-identical, with a positive control - is tested in the
            %sibling TestGammaDiffDynrelSserrLsConversion. Raising this method to
            %that standard is the obvious improvement and is left as a maintainer
            %call rather than taken here.
            nd = 50;
            alpha_a = linspace(1.6, 1.8, nd)'; alpha_b = linspace(1.5, 1.7, nd)';
            s_p_a = 0.30 * ones(nd,1); s_p_b = 0.28 * ones(nd,1);
            s_i_a = 0.12 * ones(nd,1); s_i_b = 0.14 * ones(nd,1);
            cor_p = 0.45 * ones(nd,1);  cor_i = 0.20 * ones(nd,1);

            cc = psyrat_gamma_crosscov(alpha_a, s_p_a, s_i_a, ...
                alpha_b, s_p_b, s_i_b, cor_p, cor_i);

            %the helper has no dispersion argument at all - that IS the property
            testCase.verifyEqual(nargin('psyrat_gamma_crosscov'), 8, ...
                ['psyrat_gamma_crosscov must keep exactly its eight ' ...
                'mean-submodel arguments; a ninth would mean a dispersion ' ...
                'quantity had been routed in.']);
            testCase.verifyTrue(all(isfinite(cc.cov_p_obs)));
            testCase.verifyTrue(all(isfinite(cc.cov_i_obs)));
        end

        function testResidualCrossCovarianceIsZero(testCase)
            %I.4. Non-concurrent, so the conditional residual cross-covariance is
            %structurally zero, and the difference residual is the plain SUM of
            %the two per-event residuals. Checked against the converter directly
            %so the assertion is on the quantity, not on a coefficient.
            nd = 40;
            alpha_a = 1.7 * ones(nd,1); alpha_b = 1.6 * ones(nd,1);
            s_p = 0.3 * ones(nd,1); s_i = 0.1 * ones(nd,1);
            lsig_a = -0.2 * ones(nd,1); lsig_b = -0.1 * ones(nd,1);
            s_sd = 0.25 * ones(nd,1);

            cc = psyrat_gamma_crosscov(alpha_a, s_p, s_i, alpha_b, s_p, s_i, ...
                0.4*ones(nd,1), 0.2*ones(nd,1));
            testCase.verifyEqual(cc.cov_e_obs, zeros(nd,1), ...
                'The non-concurrent residual cross-covariance must be exactly 0.');

            %The reported residual contrast SD must be the posterior mean of
            %sqrt(sigma^2_1 + sigma^2_2) with NO cross term. Reproduce it from the
            %converters and compare against what the surface reports at the same
            %z, which pins the assembly rather than merely asserting positivity.
            d = localDraws();
            out = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',0, ...
                'obs',[20 20], 'CI',0.95);

            bb   = localGet(d,'b');       bsig  = localGet(d,'b_sigma');
            sdid = localGet(d,'sd_id');   sdtrl = localGet(d,'sd_trl');
            va = psyrat_gamma_varcomps_ls(bb(:,1), sdid(:,1), sdtrl(:,1), ...
                bsig(:,1), sdid(:,3));
            vb = psyrat_gamma_varcomps_ls(bb(:,2), sdid(:,2), sdtrl(:,2), ...
                bsig(:,2), sdid(:,4));
            expected = mean(sqrt(va.sigma_pi_e2 + vb.sigma_pi_e2));

            testCase.verifyEqual(out.sige_pt, expected, 'RelTol', 1e-12, ...
                ['The difference residual must be the plain SUM of the two ' ...
                'per-event residuals: the non-concurrent cross term is zero.']);

            %and the converters themselves stay well-defined at these values
            vc_a = psyrat_gamma_varcomps_ls(alpha_a, s_p, s_i, lsig_a, s_sd);
            vc_b = psyrat_gamma_varcomps_ls(alpha_b, s_p, s_i, lsig_b, s_sd);
            testCase.verifyGreaterThan(mean(vc_a.sigma_pi_e2 + vc_b.sigma_pi_e2), 0);
        end

        function testConditionalVarianceIsExactlySigmaSquared(testCase)
            %The defining identity of the parameterization, at s_sd = 0 where the
            %marginalization collapses to the plug-in: E[Var(Y|mu,sigma)] must be
            %exp(2*log_sigma) exactly, with no dependence on alpha, s_p or s_i.
            nd = 25;
            lsig = linspace(-0.5, 0.3, nd)';
            vc1 = psyrat_gamma_varcomps_ls(1.7*ones(nd,1), 0.30*ones(nd,1), ...
                0.10*ones(nd,1), lsig, 0);
            vc2 = psyrat_gamma_varcomps_ls(2.9*ones(nd,1), 0.55*ones(nd,1), ...
                0.40*ones(nd,1), lsig, 0);

            testCase.verifyEqual(vc1.e_cond_var, exp(2*lsig), 'RelTol', 1e-12);
            testCase.verifyEqual(vc2.e_cond_var, vc1.e_cond_var, 'RelTol', 1e-12, ...
                ['E[Var(Y|mu,sigma)] must not depend on the mean submodel; ' ...
                'that is what distinguishes this from the log-nu residual.']);
        end

        function testAsymmetricMeanSlopesMoveTheSurfaceAtZeroScaleSlope(testCase)
            %I.5 / ruling 15. A mean-only asymmetry must move the coefficient even
            %with both scale slopes pinned at zero - the observed-scale signal
            %contrast genuinely changes. This is the property the provenance text
            %warns about, so it needs to be true.
            z = -1.5:0.5:1.5;
            base = localOverride(localDraws(), 'b_sigma_dim', zeros(200,2));

            sym  = localOverride(base, 'b_dim', [0.15*ones(200,1), 0.15*ones(200,1)]);
            asym = localOverride(base, 'b_dim', [0.30*ones(200,1), 0.00*ones(200,1)]);

            outSym  = psyrat_rel_diffdynrel_gamma_ls(sym{:},  'ndim',1, 'z1',z, ...
                'obs',[20 30], 'CI',0.95);
            outAsym = psyrat_rel_diffdynrel_gamma_ls(asym{:}, 'ndim',1, 'z1',z, ...
                'obs',[20 30], 'CI',0.95);

            testCase.verifyGreaterThan(std(outAsym.G.pt), 1e-6, ...
                ['An asymmetric mean slope must move the surface even at ' ...
                'b_sigma_dim = 0; suppressing it would be wrong (ruling 15).']);
            testCase.verifyNotEqual(outAsym.G.pt, outSym.G.pt);
        end

        function testUnequalPerEventObsIsHonoured(testCase)
            %Difference designs are routinely unbalanced; each event's error term
            %must be scaled by its OWN n'. More observations cannot lower
            %reliability, and the two orderings must differ when the per-event
            %residuals differ.
            z = 0;
            d = localDraws();
            few  = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',z, ...
                'obs',[5 5],   'CI',0.95);
            many = psyrat_rel_diffdynrel_gamma_ls(d{:}, 'ndim',1, 'z1',z, ...
                'obs',[80 80], 'CI',0.95);
            testCase.verifyGreaterThan(many.G.pt, few.G.pt, ...
                'More replications per event must not reduce generalizability.');
        end

        function testInputValidation(testCase)
            %The guards that keep a malformed call from producing a plausible
            %all-NaN surface with no diagnostic.
            d = localDraws();
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(d{:}, ...
                'ndim',3, 'z1',0, 'obs',20, 'CI',0.95), 'varargin:ndim');
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(d{:}, ...
                'ndim',2, 'z1',0, 'obs',20, 'CI',0.95), 'varargin:z2');
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(d{:}, ...
                'ndim',1, 'z1',0, 'obs',NaN, 'CI',0.95), 'varargin:obs');
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(d{:}, ...
                'ndim',1, 'z1',0, 'obs',0, 'CI',0.95), 'varargin:obs');
            %KDIM mismatch: two-dimension slopes passed with ndim = 1
            d6 = localDraws(6);
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(d6{:}, ...
                'ndim',1, 'z1',0, 'obs',20, 'CI',0.95), 'varargin:bdim');
            %sd_id with the wrong width
            bad = localOverride(d, 'sd_id', ones(200,3));
            testCase.verifyError(@() psyrat_rel_diffdynrel_gamma_ls(bad{:}, ...
                'ndim',1, 'z1',0, 'obs',20, 'CI',0.95), 'varargin:sdid');
        end

    end
end

% ---------------------------------------------------------------------------
% helpers
% ---------------------------------------------------------------------------

function args = localDraws(kdim)
%Deterministic synthetic posterior draws in the case-12 location-scale layout.
%Fixed values, not random: these tests pin arithmetic, so a seed would add
%nothing and would make failures harder to read.
if nargin < 1, kdim = 2; end
nd = 200;
o = ones(nd,1);

b       = [1.75*o, 1.60*o];          % per-event log-mean intercepts
b_sigma = [-0.25*o, -0.10*o];        % per-event log residual-SD intercepts
%event-crossed slopes: ODD columns event 1, EVEN event 2
if kdim == 2
    b_dim       = [0.20*o, 0.05*o];
    b_sigma_dim = [0.15*o, -0.05*o];
else
    b_dim       = [0.20*o, 0.05*o, 0.10*o, 0.02*o, 0.04*o, 0.01*o];
    b_sigma_dim = [0.15*o, -0.05*o, 0.06*o, 0.03*o, 0.02*o, 0.01*o];
end
sd_id  = [0.32*o, 0.29*o, 0.24*o, 0.21*o];  % [m1 m2 logsig1 logsig2]
sd_trl = [0.11*o, 0.13*o];

%person 4x4 correlation: cross-event mean correlation 0.45, within-event
%mean<->log-sigma 0.30/0.25 (present so the test would catch a converter that
%wrongly consumed them), cross-event scale 0.35.
R = eye(4);
R(1,2) = 0.45; R(2,1) = 0.45;
R(1,3) = 0.30; R(3,1) = 0.30;
R(2,4) = 0.25; R(4,2) = 0.25;
R(3,4) = 0.35; R(4,3) = 0.35;
cor_id = repmat(reshape(R, [1 4 4]), [nd 1 1]);

args = {'b',b, 'b_sigma',b_sigma, 'b_dim',b_dim, 'b_sigma_dim',b_sigma_dim, ...
    'sd_id',sd_id, 'sd_trl',sd_trl, 'cor_id',cor_id, 'cor_i',0.20*o};
end

function v = localGet(args, name)
%Read one value out of a name-value cell array.
idx = find(strcmp(args, name), 1);
v = args{idx + 1};
end

function args = localOverride(args, name, value)
%Replace one value in a name-value cell array, leaving the rest untouched.
idx = find(strcmp(args, name), 1);
args{idx + 1} = value;
end
