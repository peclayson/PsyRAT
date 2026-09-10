classdef TestGammaDodLsConversion < matlab.unittest.TestCase
    %TESTGAMMADODLSCONVERSION Offline checks for the LOCATION-SCALE gamma
    % four-event difference-of-differences path (analysis 9, gammascale = 2).
    %
    % Sibling of TestGammaDodConversion, which covers the log-nu twin. These
    % tests are INDEPENDENT of the production extractor
    % (psyrat_gamma_extract_dodiff_ls), which is a local function of
    % psyrat_computevarcomp.m and cannot be called directly. They exercise the
    % public helpers it composes -- psyrat_gamma_varcomps_ls per cell,
    % psyrat_gamma_crosscov per cell pair, psyrat_dodiffrel for the
    % c = [1 -1 -1 1] contrast -- against an independent transcription, and pin
    % the extractor's own contribution (the gamma_disp stamp that routes the
    % provenance writer) through the dispatch tests at the end.
    %
    % WHAT THE INCREMENT TURNS ON, and what these tests are therefore built to
    % discriminate: analysis 9 is a ONE-FACET design replicated over four cells,
    % so it needs NO new converter -- the existing one-facet
    % psyrat_gamma_varcomps_ls serves every cell diagonal. The load-bearing
    % claims are (a) the mean-surface components are bit-identical across the two
    % parameterizations, which is what licenses reusing psyrat_gamma_crosscov
    % unchanged across all six pairs, and (b) the residual carries no amplitude
    % term, which is the whole point of the parameterization. Both are asserted
    % below in forms that FAIL if the wrong converter is wired in.

    methods (Test)

        % ---------------------------------------------------------------
        % Assembly parity
        % ---------------------------------------------------------------

        function testFourCellAssemblyMatchesIndependentDodMath(testCase)
            % Random pseudo-posterior draws of the per-cell log-scale parameters
            % under the LOCATION-SCALE parameterization; assemble via the public
            % helpers -> psyrat_dodiffrel and compare G, D, both ICCs and the
            % universe-score variance against an independent transcription of the
            % DoD math that computes the observed-scale components from scratch.
            rng(24601, 'twister');
            nd = 3000; C = 4;
            pairs = nchoosek(1:C, 2);
            obs   = [18 22 15 20];
            cvec  = [1 -1 -1 1];

            alpha = 2.5  + 0.30 .* randn(nd, C);
            sd_p  = abs(0.35 + 0.10 .* randn(nd, C));
            sd_i  = abs(0.18 + 0.06 .* randn(nd, C));
            lsig  = 0.60 + 0.25 .* randn(nd, C);   % per-cell log residual-SD
            s_sd  = abs(0.45 + 0.10 .* randn(nd, C));
            cor_p = 0.30 + 0.50 .* rand(nd, size(pairs,1));
            cor_i = 0.00 + 0.40 .* rand(nd, size(pairs,1));

            [idv, tlv, bsig] = localAssembleLs(alpha, sd_p, sd_i, lsig, s_sd, ...
                cor_p, cor_i, pairs);
            gen = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', obs, 'est', 'gen', 'ci', .95, 'cvec', cvec);
            dep = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', obs, 'est', 'dep', 'ci', .95, 'cvec', cvec);

            ref = localRefDodLs(alpha, sd_p, sd_i, lsig, s_sd, cor_p, cor_i, ...
                pairs, obs, cvec);

            testCase.verifyEqual(gen.pt,         ref.G_pt,      'AbsTol', 1e-9);
            testCase.verifyEqual(dep.pt,         ref.D_pt,      'AbsTol', 1e-9);
            testCase.verifyEqual(gen.ll,         ref.G_ll,      'AbsTol', 1e-9);
            testCase.verifyEqual(gen.ul,         ref.G_ul,      'AbsTol', 1e-9);
            testCase.verifyEqual(dep.ll,         ref.D_ll,      'AbsTol', 1e-9);
            testCase.verifyEqual(dep.ul,         ref.D_ul,      'AbsTol', 1e-9);
            testCase.verifyEqual(gen.icc_pt,     ref.ICCrel_pt, 'AbsTol', 1e-9);
            testCase.verifyEqual(dep.icc_pt,     ref.ICCabs_pt, 'AbsTol', 1e-9);
            testCase.verifyEqual(gen.uni_var_pt, ref.uni_pt,    'AbsTol', 1e-9);
        end

        % ---------------------------------------------------------------
        % The claims that license reusing the existing helpers
        % ---------------------------------------------------------------

        function testMeanSurfaceMatricesAreParameterizationInvariant(testCase)
            % THE structural claim of this increment. Given the SAME log-mean
            % submodel draws, the assembled 4x4 person and trial (co)variance
            % matrices must be bit-identical under the two parameterizations,
            % because both are built from the log-mean submodel alone. This is
            % why psyrat_gamma_crosscov is reused unchanged across all six pairs
            % and why no psyrat_gamma_varcomps_dod_ls was written.
            %
            % Bit-identity (0 tolerance), not a tolerance: these are the same
            % arithmetic on the same inputs, not two routes that happen to agree.
            rng(4242, 'twister');
            nd = 500; C = 4; pairs = nchoosek(1:C,2);

            alpha = 1.6  + 0.20 .* randn(nd, C);
            sd_p  = abs(0.40 + 0.08 .* randn(nd, C));
            sd_i  = abs(0.20 + 0.05 .* randn(nd, C));
            cor_p = 0.20 + 0.40 .* rand(nd, size(pairs,1));
            cor_i = 0.10 + 0.30 .* rand(nd, size(pairs,1));

            % Scale submodels chosen to be genuinely different, so an accidental
            % agreement cannot come from the residuals coinciding.
            lsig  = 0.55 + 0.20 .* randn(nd, C);
            s_sd  = abs(0.40 + 0.10 .* randn(nd, C));
            nu    = exp(log(9) + 0.30 .* randn(nd, C));
            s_v   = abs(1.00 + 0.20 .* randn(nd, C));
            covpv = 0.05 .* randn(nd, C);

            [idvLs, tlvLs, bsigLs] = localAssembleLs(alpha, sd_p, sd_i, lsig, ...
                s_sd, cor_p, cor_i, pairs);
            [idvNu, tlvNu, bsigNu] = localAssembleNu(alpha, sd_p, sd_i, nu, ...
                s_v, covpv, cor_p, cor_i, pairs);

            testCase.verifyEqual(idvLs, idvNu, ...
                'Person (co)variance matrix must not depend on the scale parameterization.');
            testCase.verifyEqual(tlvLs, tlvNu, ...
                'Trial (co)variance matrix must not depend on the scale parameterization.');

            % Non-vacuous: the RESIDUAL must differ, otherwise the test above
            % would pass on two identical assemblies and prove nothing.
            testCase.verifyGreaterThan(max(abs(bsigLs(:) - bsigNu(:))), 1e-3, ...
                ['The two parameterizations were configured to give the SAME ' ...
                'residual, so the invariance check above is vacuous.']);
        end

        function testAllSixCrossCellCovariancesAreScaleIndependent(testCase)
            % psyrat_gamma_crosscov is reused UNCHANGED for every one of the six
            % cell pairs. Its whole argument list is mean-submodel quantities, so
            % it has no dispersion argument to get wrong -- checked pair by pair
            % rather than inferred from the assembled matrix, so a bug that
            % cancelled in the assembly would still be caught.
            rng(777, 'twister');
            nd = 200; C = 4; pairs = nchoosek(1:C,2);
            testCase.verifyEqual(size(pairs,1), 6, ...
                'A four-cell design has exactly six cross-cell pairs.');

            alpha = 1.5 + 0.15 .* randn(nd, C);
            sd_p  = abs(0.38 + 0.06 .* randn(nd, C));
            sd_i  = abs(0.19 + 0.04 .* randn(nd, C));
            cor_p = 0.25 + 0.40 .* rand(nd, 6);
            cor_i = 0.05 + 0.30 .* rand(nd, 6);

            for pp = 1:6
                j = pairs(pp,1); k = pairs(pp,2);
                cc = psyrat_gamma_crosscov(alpha(:,j), sd_p(:,j), sd_i(:,j), ...
                                           alpha(:,k), sd_p(:,k), sd_i(:,k), ...
                                           cor_p(:,pp), cor_i(:,pp));
                % Closed form for jointly lognormal cell means; the residual
                % cross-covariance is structurally 0 for a non-concurrent design.
                ey = exp(alpha(:,j) + 0.5.*(sd_p(:,j).^2 + sd_i(:,j).^2)) ...
                   .* exp(alpha(:,k) + 0.5.*(sd_p(:,k).^2 + sd_i(:,k).^2));
                expP = ey .* (exp(cor_p(:,pp).*sd_p(:,j).*sd_p(:,k)) - 1);
                expI = ey .* (exp(cor_i(:,pp).*sd_i(:,j).*sd_i(:,k)) - 1);
                testCase.verifyEqual(cc.cov_p_obs, expP, 'RelTol', 1e-12, ...
                    sprintf('pair %d-%d person covariance', j, k));
                testCase.verifyEqual(cc.cov_i_obs, expI, 'RelTol', 1e-12, ...
                    sprintf('pair %d-%d trial covariance', j, k));
                testCase.verifyEqual(cc.cov_e_obs, zeros(nd,1), ...
                    sprintf('pair %d-%d residual covariance must be 0', j, k));
            end
        end

        function testResidualCarriesNoAmplitudeUnderLocationScale(testCase)
            % The decisive discriminator, and the reason this parameterization
            % exists. Shifting every cell's log-mean intercept by a constant must
            % leave the location-scale residual EXACTLY unchanged, because
            % e_cond_var = exp(2*log_sigma + 2*s_sd^2) contains no alpha. Under
            % log-nu the same shift scales the residual by exp(2*shift), since
            % Var(Y|mu,nu) = 2*mu^2/nu is mean-coupled.
            %
            % A test that only checked the LS side would pass if someone wired
            % the log-nu converter in and separately zeroed alpha, so the log-nu
            % side is asserted to MOVE.
            alpha = log([5.0 4.5 5.2 4.8]);
            sd_p  = [0.40 0.38 0.42 0.39];
            sd_i  = [0.20 0.18 0.22 0.19];
            lsig  = [0.60 0.95 0.55 0.90];
            s_sd  = [0.45 0.50 0.42 0.48];
            nu    = [7 9 8 10];
            shift = 0.35;

            for c = 1:4
                a = psyrat_gamma_varcomps_ls(alpha(c),       sd_p(c), sd_i(c), lsig(c), s_sd(c));
                b = psyrat_gamma_varcomps_ls(alpha(c)+shift, sd_p(c), sd_i(c), lsig(c), s_sd(c));
                testCase.verifyEqual(b.e_cond_var, a.e_cond_var, ...
                    sprintf(['cell %d: the location-scale conditional variance ' ...
                    'must not depend on the log-mean intercept.'], c));

                % The log-nu twin must move, by exactly exp(2*shift).
                an = psyrat_gamma_varcomps(alpha(c),       sd_p(c), sd_i(c), nu(c));
                bn = psyrat_gamma_varcomps(alpha(c)+shift, sd_p(c), sd_i(c), nu(c));
                testCase.verifyEqual(bn.e_cond_var, an.e_cond_var*exp(2*shift), ...
                    'RelTol', 1e-12, sprintf(['cell %d: the log-nu conditional ' ...
                    'variance is mean-coupled and must scale by exp(2*shift); ' ...
                    'if it does not, the check above is vacuous.'], c));
            end
        end

        % ---------------------------------------------------------------
        % Residual encoding and the non-concurrent contract
        % ---------------------------------------------------------------

        function testPerCellResidualEncodingRoundTripsExactly(testCase)
            % psyrat_dodiffrel reads er_var as a LOG SD and squares exp(.), so the
            % extractor stores b_sigma = 0.5*log(sigma_pi_e2). Pin the round trip:
            % exp(0.5*log(v))^2 == v to machine precision, per cell.
            alpha = log([5.0 4.5 5.2 4.8]);
            sd_p  = [0.40 0.38 0.42 0.39];
            sd_i  = [0.20 0.18 0.22 0.19];
            lsig  = [0.60 0.95 0.55 0.90];
            s_sd  = [0.45 0.50 0.42 0.48];

            for c = 1:4
                vc = psyrat_gamma_varcomps_ls(alpha(c), sd_p(c), sd_i(c), lsig(c), s_sd(c));
                encoded = 0.5*log(vc.sigma_pi_e2);
                testCase.verifyEqual(exp(encoded).^2, vc.sigma_pi_e2, ...
                    'RelTol', 1e-12, sprintf('cell %d residual encoding', c));
            end
        end

        function testNonconcurrentDodResidualIsTheWeightedCellSum(testCase)
            % The residual enters psyrat_dodiffrel as a per-cell VECTOR, so the
            % unscaled error contrast variance is sum_c c_c^2 * sigma_pi_e2_c with
            % no residual cross-covariance. That zero is structural for this
            % design, not an assumption: correlating the person SCALE effects
            % across cells correlates the residual VARIANCES, not the residuals.
            C = 4; pairs = nchoosek(1:C,2);
            alpha = log([5.0 4.5 5.2 4.8]);
            sd_p  = [0.40 0.38 0.42 0.39];
            sd_i  = [0.20 0.18 0.22 0.19];
            lsig  = [0.60 0.95 0.55 0.90];
            s_sd  = [0.45 0.50 0.42 0.48];
            cvec  = [1 -1 -1 1];

            [idv, tlv, bsig] = localAssembleLs(alpha, sd_p, sd_i, lsig, s_sd, ...
                0.5*ones(1,6), 0.2*ones(1,6), pairs);
            out = psyrat_dodiffrel('bp', idv, 'bt', tlv, 'er_var', bsig, ...
                'obs', [20 20 20 20], 'est', 'gen', 'ci', .95, 'cvec', cvec);

            sig_pie2 = zeros(1,C);
            for c = 1:C
                vc = psyrat_gamma_varcomps_ls(alpha(c), sd_p(c), sd_i(c), lsig(c), s_sd(c));
                sig_pie2(c) = vc.sigma_pi_e2;
            end
            testCase.verifyEqual(out.er_contrast_var_pt, sum((cvec.^2).*sig_pie2), ...
                'RelTol', 1e-12);
        end

        function testMarginalCorrectionFactorIsExpTwoSsdSquared(testCase)
            % The extractor stores disp_factor = exp(2*s_sd^2), the ratio of the
            % MARGINAL conditional variance to the homogeneous-residual plug-in.
            % Pin both halves of that claim: the ratio itself, and that s_sd = 0
            % recovers the plug-in exactly (factor 1). This is the location-scale
            % correction; it is NOT the log-nu twin's exp(0.5*s_v^2 - 2*cov_pv),
            % and the two must never be interchanged.
            alpha = 1.7; sd_p = 0.40; sd_i = 0.20; lsig = 0.65;
            for s_sd = [0 0.25 0.45 0.70]
                vc   = psyrat_gamma_varcomps_ls(alpha, sd_p, sd_i, lsig, s_sd);
                plug = psyrat_gamma_varcomps_ls(alpha, sd_p, sd_i, lsig, 0);
                testCase.verifyEqual(vc.e_cond_var / plug.e_cond_var, ...
                    exp(2*s_sd^2), 'RelTol', 1e-12, ...
                    sprintf('marginal correction factor at s_sd = %g', s_sd));
            end
            homog = psyrat_gamma_varcomps_ls(alpha, sd_p, sd_i, lsig, 0);
            testCase.verifyEqual(homog.e_cond_var, exp(2*lsig), 'RelTol', 1e-12, ...
                'At s_sd = 0 the conditional variance is the plug-in exp(2*log_sigma).');
        end

        % ---------------------------------------------------------------
        % Provenance dispatch -- how the local-function extractor is pinned
        % ---------------------------------------------------------------

        function testFourCellStampReachesTheDodBranch(testCase)
            % The extractor cannot be called directly, so its contribution is
            % pinned through the stamp it writes: a four-column gamma_disp with
            % is_ls + is_dod must reach the four-cell provenance branch and report
            % all four cells, not a pooled median.
            rel   = localDodLsRel();
            lines = psyrat_provenance_lines(struct('rel', rel));
            txt   = strjoin(lines, ' | ');

            testCase.verifySubstring(txt, 'LOCATION-SCALE');
            testCase.verifySubstring(txt, 's_sd) by cell');
            testCase.verifySubstring(txt, 'difference-of-differences (location-scale)');
            % All four medians present, in order, with the planted values.
            testCase.verifySubstring(txt, '0.300, 0.400, 0.500, 0.600');
            % Never the log-nu machinery: this model has no nu at all.
            testCase.verifyEmpty(strfind(txt, 's_v'), ...
                'A location-scale run must never report a person log-nu SD.');
            testCase.verifyEmpty(strfind(txt, 'exp(0.5*s_v^2'), ...
                'A location-scale run must never report the log-nu F correction.');
        end

        function testFourCellStampDoesNotReachTheTwoEventBranch(testCase)
            % The two-event branch requires exactly two columns and its text says
            % "by event". A four-cell layout must not reach it -- pooling four
            % cells into a two-event report would misdescribe the design.
            rel   = localDodLsRel();
            txt   = strjoin(psyrat_provenance_lines(struct('rel', rel)), ' | ');
            testCase.verifyEmpty(strfind(txt, 'by event'), ...
                'The four-cell layout must not be reported through the two-event branch.');
        end

        function testTwoEventStampStillReachesItsOwnBranch(testCase)
            % Positive counterpart, so the pair above cannot pass vacuously by
            % rejecting every location-scale difference run. A genuine two-event
            % is_diff stamp must still reach the two-event branch untouched.
            rel = localDodLsRel();
            rel.out.gamma_disp = struct( ...
                's_sd',        [0.30 0.40] .* ones(40,2), ...
                'rho_ps',      [-0.20 -0.10] .* ones(40,2), ...
                'disp_factor', exp(2 .* ([0.30 0.40].^2)) .* ones(40,2), ...
                'is_ls',       true, ...
                'is_diff',     true);
            txt = strjoin(psyrat_provenance_lines(struct('rel', rel)), ' | ');
            testCase.verifySubstring(txt, 'by event');
            testCase.verifyEmpty(strfind(txt, 'by cell'), ...
                'A two-event run must not be reported through the four-cell branch.');
        end

        function testUnstampedFourColumnLayoutFallsThrough(testCase)
            % Fail closed on the stamp, not on the column count. A four-column
            % LS layout WITHOUT the is_dod stamp must emit no location-scale
            % dispersion lines rather than being absorbed by any branch.
            rel = localDodLsRel();
            rel.out.gamma_disp = rmfield(rel.out.gamma_disp, 'is_dod');
            txt = strjoin(psyrat_provenance_lines(struct('rel', rel)), ' | ');
            testCase.verifyEmpty(strfind(txt, 'by cell'));
            testCase.verifyEmpty(strfind(txt, 'by event'));
            testCase.verifyEmpty(strfind(txt, 'Gamma scale parameterization: LOCATION-SCALE'));
        end

        function testStampedButMalformedLayoutFallsThrough(testCase)
            % A stamped struct whose columns do not match must fall through
            % rather than error on gd.s_sd(:,4). This is the guard that keeps a
            % future layout change from producing an index crash inside a table
            % export.
            rel = localDodLsRel();
            rel.out.gamma_disp.s_sd = rel.out.gamma_disp.s_sd(:,1:3);
            txt = strjoin(psyrat_provenance_lines(struct('rel', rel)), ' | ');
            testCase.verifyEmpty(strfind(txt, 'by cell'));
        end

        function testLogNuDodRunReachesNoLocationScaleBranch(testCase)
            % A gammascale = 1 DoD run stamps s_v / cov_pv and must reach the
            % log-nu chain, never any location-scale branch. Without this the
            % suite could pass while every gamma run printed location-scale text.
            rel = localDodLsRel();
            rel.out.gamma_disp = struct( ...
                's_v',         1.10 .* ones(40,4), ...
                'cov_pv',      0.05 .* ones(40,4), ...
                'disp_factor', exp(0.5*1.10^2 - 2*0.05) .* ones(40,4));
            txt = strjoin(psyrat_provenance_lines(struct('rel', rel)), ' | ');
            testCase.verifyEmpty(strfind(txt, 'LOCATION-SCALE'), ...
                'A log-nu run must never report the location-scale estimand.');
            testCase.verifyEmpty(strfind(txt, 'by cell'));
            % Positive counterpart, so this cannot pass by emitting nothing at
            % all: the log-nu chain must still report its own quantities.
            testCase.verifySubstring(txt, 'person log-nu SD, s_v');
            testCase.verifySubstring(txt, 'F = exp(0.5*s_v^2 - 2*cov)');
        end

    end
end

% =========================================================================
function rel = localDodLsRel()
% Minimal REL carrying what psyrat_provenance_lines reads for a location-scale
% gamma difference-of-differences run. s_sd values are distinct per cell so a
% pooled median could not accidentally reproduce them.
s_sd = [0.30 0.40 0.50 0.60];
rel = struct();
rel.family   = 'gamma';
rel.analysis = 'ic_dodiff';
rel.seed     = 12345;
rel.out.gamma_disp = struct( ...
    's_sd',        s_sd  .* ones(40,4), ...
    'rho_ps',      [-0.20 -0.10 -0.30 -0.15] .* ones(40,4), ...
    'disp_factor', exp(2 .* (s_sd.^2)) .* ones(40,4), ...
    'is_ls',       true, ...
    'is_dod',      true);
end

% =========================================================================
function [id_varcov, trl_varcov, b_sigma] = localAssembleLs(alpha, sd_p, sd_i, ...
    log_sigma, s_sd, cor_p, cor_i, pairs)
% Replicate the production psyrat_gamma_extract_dodiff_ls assembly using the
% public helpers: per-cell diagonals from psyrat_gamma_varcomps_ls, cross-cell
% off-diagonals from psyrat_gamma_crosscov, per-cell residual encoded as
% b_sigma = 0.5*log(sigma_pi_e2).
nd = size(alpha,1); C = size(alpha,2);
id_varcov  = zeros(nd,C,C);
trl_varcov = zeros(nd,C,C);
b_sigma    = zeros(nd,C);
for c = 1:C
    vc = psyrat_gamma_varcomps_ls(alpha(:,c), sd_p(:,c), sd_i(:,c), ...
        log_sigma(:,c), s_sd(:,c));
    id_varcov(:,c,c)  = vc.sigma_p2;
    trl_varcov(:,c,c) = vc.sigma_i2;
    b_sigma(:,c)      = 0.5 .* log(vc.sigma_pi_e2);
end
[id_varcov, trl_varcov] = localFillPairs(id_varcov, trl_varcov, alpha, sd_p, ...
    sd_i, cor_p, cor_i, pairs);
end

% =========================================================================
function [id_varcov, trl_varcov, b_sigma] = localAssembleNu(alpha, sd_p, sd_i, ...
    nu, s_v, cov_pv, cor_p, cor_i, pairs)
% The log-nu twin of localAssembleLs, used only by the invariance test.
nd = size(alpha,1); C = size(alpha,2);
id_varcov  = zeros(nd,C,C);
trl_varcov = zeros(nd,C,C);
b_sigma    = zeros(nd,C);
for c = 1:C
    vc = psyrat_gamma_varcomps(alpha(:,c), sd_p(:,c), sd_i(:,c), nu(:,c), ...
        s_v(:,c), cov_pv(:,c));
    id_varcov(:,c,c)  = vc.sigma_p2;
    trl_varcov(:,c,c) = vc.sigma_i2;
    b_sigma(:,c)      = 0.5 .* log(vc.sigma_pi_e2);
end
[id_varcov, trl_varcov] = localFillPairs(id_varcov, trl_varcov, alpha, sd_p, ...
    sd_i, cor_p, cor_i, pairs);
end

% =========================================================================
function [id_varcov, trl_varcov] = localFillPairs(id_varcov, trl_varcov, ...
    alpha, sd_p, sd_i, cor_p, cor_i, pairs)
% Shared six-pair cross-cell fill. Deliberately common to both assemblies: the
% invariance test's claim is precisely that this step takes no scale argument.
for pp = 1:size(pairs,1)
    j = pairs(pp,1); k = pairs(pp,2);
    cc = psyrat_gamma_crosscov(alpha(:,j), sd_p(:,j), sd_i(:,j), ...
                               alpha(:,k), sd_p(:,k), sd_i(:,k), ...
                               cor_p(:,pp), cor_i(:,pp));
    id_varcov(:,j,k)  = cc.cov_p_obs; id_varcov(:,k,j)  = cc.cov_p_obs;
    trl_varcov(:,j,k) = cc.cov_i_obs; trl_varcov(:,k,j) = cc.cov_i_obs;
end
end

% =========================================================================
function ref = localRefDodLs(alpha, sd_p, sd_i, log_sigma, s_sd, cor_p, cor_i, ...
    pairs, obs, cvec)
% Independent transcription of the four-cell DoD math under the LOCATION-SCALE
% parameterization. The mean-surface components are the shared lognormal
% expressions; only the conditional variance differs from the log-nu reference in
% TestGammaDodConversion, where it is (2/nu)*exp(2*alpha + 2*ssum). Here it is
% exp(2*log_sigma + 2*s_sd^2) -- no alpha, no ssum.
nd = size(alpha,1); C = size(alpha,2); w = cvec(:); n = obs(:);
np = size(pairs,1);

sig_p2 = zeros(nd,C); sig_i2 = zeros(nd,C); sig_pie2 = zeros(nd,C);
for c = 1:C
    vp = sd_p(:,c).^2; vi = sd_i(:,c).^2; vt = vp + vi;
    base = exp(2.*alpha(:,c) + vt);
    sig_p2(:,c)   = base .* (exp(vp) - 1);
    sig_i2(:,c)   = base .* (exp(vi) - 1);
    var_mu        = base .* (exp(vt) - 1);
    e_cond        = exp(2.*log_sigma(:,c) + 2.*s_sd(:,c).^2);
    sig_pie2(:,c) = max(var_mu + e_cond - sig_p2(:,c) - sig_i2(:,c), 0);
end

cov_p = zeros(nd,np); cov_i = zeros(nd,np);
for pp = 1:np
    j = pairs(pp,1); k = pairs(pp,2);
    ey = exp(alpha(:,j) + 0.5.*(sd_p(:,j).^2 + sd_i(:,j).^2)) ...
       .* exp(alpha(:,k) + 0.5.*(sd_p(:,k).^2 + sd_i(:,k).^2));
    cov_p(:,pp) = ey .* (exp(cor_p(:,pp).*sd_p(:,j).*sd_p(:,k)) - 1);
    cov_i(:,pp) = ey .* (exp(cor_i(:,pp).*sd_i(:,j).*sd_i(:,k)) - 1);
end

G = zeros(nd,1); D = zeros(nd,1); ICCr = zeros(nd,1); ICCa = zeros(nd,1); U = zeros(nd,1);
for d = 1:nd
    Vp = diag(sig_p2(d,:)); Vi = diag(sig_i2(d,:)); Vpie = diag(sig_pie2(d,:));
    for pp = 1:np
        j = pairs(pp,1); k = pairs(pp,2);
        Vp(j,k)=cov_p(d,pp); Vp(k,j)=cov_p(d,pp);
        Vi(j,k)=cov_i(d,pp); Vi(k,j)=cov_i(d,pp);
    end
    u        = w'*Vp*w;
    trial_ab = localQformDivN(Vi, w, n);
    rel_err  = localQformDivN(Vpie, w, n);
    U(d)   = u;
    G(d)   = u / (u + rel_err);
    D(d)   = u / (u + rel_err + trial_ab);
    ICCr(d)= u / (u + w'*Vpie*w);
    ICCa(d)= u / (u + w'*Vpie*w + w'*Vi*w);
end

qg = quantile(G,[.025 .975]); qd = quantile(D,[.025 .975]);
ref = struct('G_pt',mean(G),'G_ll',qg(1),'G_ul',qg(2), ...
             'D_pt',mean(D),'D_ll',qd(1),'D_ul',qd(2), ...
             'ICCrel_pt',mean(ICCr),'ICCabs_pt',mean(ICCa),'uni_pt',mean(U));
end

% =========================================================================
function v = localQformDivN(V, w, n)
% Reference qform_div_n: diagonal / n_j, off-diagonal / harmonic-mean(n_j,n_k).
C = numel(w); v = 0;
for j = 1:C, v = v + (w(j)^2) * V(j,j) / n(j); end
for j = 1:C-1
    for k = j+1:C
        nharm = 2 / ((1/n(j)) + (1/n(k)));
        v = v + 2*w(j)*w(k)*V(j,k)/nharm;
    end
end
end
