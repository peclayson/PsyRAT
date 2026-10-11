classdef TestGammaDynrelSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL dynamic/conditional reliability designs: analysis 26
    % (ic_dynrel_sserrvar, one facet) and analysis 27 (ic_dynrel_sserrvar_trt,
    % trial + occasion), both with family='gamma'.
    %
    % WHY THIS EXISTS. The offline tests prove the calculators
    % (TestGammaDynrelSserrReliability) and the summary wiring
    % (TestDynrelSubjNondiffReliability), but neither exercises the ESTIMATOR:
    % that psyrat_computevarcomp initializes the gamma REL.out slots for
    % ind_bs/ind_nu/ssinfo, extracts those per-person effects from the fit, and
    % appends ssinfo per stratum. Before this work a gamma run of analysis 26/27
    % ERRORED on the first stratum ("Reference to non-existent field 'ssinfo'"),
    % because the {end+1} append was guarded on the analysis but not on the
    % family. This test is what pins that plumbing.
    %
    % These are also the FIRST live recovery tests for analyses 26/27 in either
    % family: the development tracker records that their Gaussian recovery passed
    % during development but was deliberately not re-run after the salvage.
    %
    % WHAT IS ASSERTED. The planted per-person dispersion ordering must come back:
    % the recovered per-participant conditional coefficient must correlate with
    % the truth computed from the planted (u_p, v_p) at each participant's own z.
    % That is the assertion that fails if ind_bs/ind_nu are misaligned against
    % id2, or if the conditional conversion is wrong. The population surface is
    % not asserted here - it is the marginal estimand and is covered by
    % TestGammaDynrelRecovery.
    %
    % PLANTED TRUTH (2026-09-13 correction, bounded 3A / G89 / Q12 prerequisite;
    % see the frozen 2C S20 recovery reference (local development evidence)
    % and its scientific contract). The per-participant residual is the accepted
    % population-pooled target, not the observation variance alone:
    %   analysis 26 (one facet):  E_p = PI(z) + C_p,
    %       C_p = 2 exp(2(alpha_z + u_p) + 2 s_I^2 - (lognu_z + v_p)),
    %       G_p = P / (P + E_p / n_i);
    %   analysis 27 (two facets): E_p = R(z) + C_p,
    %       C_p = 2 exp(2(alpha_z + u_p) + 2 W - (lognu_z + v_p)),
    %       W = s_O^2 + s_T^2 + s_PO^2 + s_PT^2 + s_OT^2,
    %       G_p = P / (P + PI/n_i + PO/n_o + E_p/(n_i n_o)),  n_o = 1 (CES).
    % P, PI, PO and R are the population mean components at the participant's
    % own z, formed from explicit support-union products of q_A = expm1(v_A)
    % (no production converter or reliability kernel is called to build the
    % expected values). The earlier truth omitted PI / R and is superseded.
    %
    % GENERATOR. Observations are exact continuous gamma draws,
    % gamrnd(nu_p/2, 2 mu/nu_p) with nu_p = exp(lognu_p) unrounded, so the
    % fitted continuous log-nu likelihood is the generator's own model and the
    % planted (lognu0, b_nu, s_V) are exact parameter truths. The earlier
    % integer-df chi-square generator (max(2, round(nu))) is superseded; its
    % logs remain historical evidence for that older DGP. Tolerances, design
    % sizes, seeds and sampler settings are unchanged.
    %
    % This correction is a bounded prerequisite for 3A / G89 / Q12. It does not
    % by itself establish repeated-dataset uncertainty coverage.
    %
    % Dependency setup and the R-hat gate both use assumptions, so a missing
    % CmdStan, a missing gamrnd, or a non-converged fit reports Incomplete,
    % never a false PASS. When PSYRAT_2C_RECOVERY_CAPTURE_DIR is set, each fit
    % (REL, truth, data, raw CmdStan output) is captured through
    % PsyRAT2CRecoveryCapture BEFORE the convergence assumption, so a skipped
    % run still leaves its evidence, and the class teardown skips its scratch
    % cleanup so the raw output also survives a failed capture or test.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            % When 2C capture is requested, leave the raw scratch untouched so it
            % survives even a failed capture or a failed test; root's external
            % artifact guard archives/restores the owned scratch afterwards.
            if ~isempty(PsyRAT2CRecoveryCapture.captureDir())
                fprintf('[gamma dynrel sserr recovery] %s set; scratch cleanup skipped.\n', ...
                    PsyRAT2CRecoveryCapture.EnvVar);
                return;
            end
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)

        function testOneFacetSubjectLevelGammaDynrelRecovers(testCase)
            % Continuous-df gamma draws need gamrnd (Statistics Toolbox); the
            % static sibling TestGammaSserrRecovery guards the same way.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            alpha0 = log(5.0);  b    = 0.25;    % log-mean intercept + dimension slope
            lognu0 = log(20);   b_nu = -0.35;   % log-nu intercept + dimension slope
            sP = 0.40;   % person log-mean SD
            sI = 0.25;   % trial main-effect log-mean SD
            sV = 0.30;   % person log-nu SD
            rho = 0.35;  % corr(person log-mean, person log-nu)

            nsub = 45; ntrl = 30;

            [dataTbl, truth] = localSimulateGammaDynrelSserr(nsub, ntrl, ...
                alpha0, b, lognu0, b_nu, sP, sI, sV, rho);

            rel = localRunDynrelSserrGamma(testCase, dataTbl, 'recov_gamma_dynrelss');

            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Opt-in evidence capture (no-op unless PSYRAT_2C_RECOVERY_CAPTURE_DIR
            % is set), taken BEFORE the convergence assumption below.
            PsyRAT2CRecoveryCapture.capture('dynrel26_lognu_sserr', rel, truth, dataTbl, ...
                fullfile(localCmdStanOutputDir(testCase.projectRoot()), 'recov_gamma_dynrelss'), ...
                'meta', struct('test', 'TestGammaDynrelSserrRecovery/testOneFacetSubjectLevelGammaDynrelRecovers', ...
                'gammascale', 1, 'chains', 2, 'warmup', 750, 'sampling', 750, 'seed', 12345, ...
                'nsub', nsub, 'ntrl', ntrl, 'design', 'one-facet persons x trials'));

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma dynrel sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'subject-level dynrel fit non-converged (max R-hat=%.4f >= 1.1 / ', ...
                'n_eff gate); recovery not validated in this environment.'], maxRhat));

            % ---- The estimator stored what the subject-level path needs ----
            % (this is precisely what was missing before: a gamma analysis-26 run
            % had no ind_bs/ind_nu/ssinfo at all)
            for f = {'ind_bs','ind_nu','ssinfo'}
                testCase.verifyTrue(isfield(rel.out, f{1}), sprintf( ...
                    ['REL.out.%s missing on a gamma analysis-26 run - the ', ...
                    'subject-level slots are not being initialized/extracted.'], f{1}));
            end
            indbs = cell2mat(rel.out.ind_bs{1});   % draws x NSUB
            indnu = cell2mat(rel.out.ind_nu{1});
            testCase.verifySize(indbs, [size(rel.out.mu,1) nsub]);
            testCase.verifySize(indnu, [size(rel.out.mu,1) nsub]);
            testCase.verifyEqual(height(rel.out.ssinfo{1}), nsub);

            % ---- The per-participant table is produced and is per-participant ----
            summ = psyrat_dynrel_summary(rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;
            testCase.verifyNotEmpty(T, ...
                'No per-participant table: analysis 26 degraded to analysis 11.');
            testCase.verifyEqual(height(T), nsub);
            testCase.verifyTrue(all(isfinite(T.gen_pt)));

            % ---- Recovery: the per-participant ordering matches the truth ----
            r = corr_(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel sserr recovery] per-participant G corr = %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.70, sprintf( ...
                ['Recovered per-participant coefficients correlate only %.3f with ', ...
                'the planted truth. Either the conditional conversion is wrong or ', ...
                'ind_bs/ind_nu are misaligned against id2.'], r));

            % SCALE. Asserting on the coefficient itself is nearly useless here:
            % G saturates near 1 at these trial counts, so an upper bound can never
            % fire and a lower bound tolerates a large residual error. Assert on the
            % implied RESIDUAL instead, which is unsaturated and is what a wrong
            % conditional conversion actually gets wrong.
            %
            % Invert the coefficient for each participant: from
            % G = bp / (bp + err/n) we get err/n = bp*(1-G)/G, so the ratio of
            % implied to planted error variance is comparable without needing bp.
            errRatio = localImpliedErrRatio(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel sserr recovery] implied error-variance ratio = %.3f\n', errRatio);
            testCase.verifyGreaterThan(errRatio, 0.60, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- the conditional residual is too small.'], errRatio));
            testCase.verifyLessThan(errRatio, 1.70, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- the conditional residual is too large (a person-facet leak into ', ...
                'V_w would show up here).'], errRatio));

            % genuine between-person spread (else it has collapsed to the surface)
            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3, ...
                ['No between-person spread in the per-participant coefficients; ', ...
                'the conditional conversion has degenerated to the marginal one.']);
        end

        function testTwoFacetSubjectLevelGammaDynrelRecovers(testCase)
            % The analysis-27 twin. Smaller design (a crossed two-facet gamma fit
            % is expensive) - the assertion is the same: the plumbing exists and
            % the per-participant ordering is recovered.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            alpha0 = log(5.0);  b    = 0.20;
            lognu0 = log(20);   b_nu = -0.30;
            sP = 0.35; sO = 0.20; sT = 0.15; sPO = 0.12; sPT = 0.10; sOT = 0.08;
            sV = 0.30; rho = 0.35;

            nsub = 30; nocc = 3; ntrl = 10;

            [dataTbl, truth] = localSimulateGammaDynrelSserrTrt(nsub, nocc, ntrl, ...
                alpha0, b, lognu0, b_nu, sP, sO, sT, sPO, sPT, sOT, sV, rho);

            rel = localRunDynrelSserrTrtGamma(testCase, dataTbl, 'recov_gamma_dynrelsstrt');

            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar_trt');
            testCase.verifyEqual(rel.family, 'gamma');

            % Opt-in evidence capture (no-op unless PSYRAT_2C_RECOVERY_CAPTURE_DIR
            % is set), taken BEFORE the convergence assumption below.
            PsyRAT2CRecoveryCapture.capture('dynrel27_lognu_sserr', rel, truth, dataTbl, ...
                fullfile(localCmdStanOutputDir(testCase.projectRoot()), 'recov_gamma_dynrelsstrt'), ...
                'meta', struct('test', 'TestGammaDynrelSserrRecovery/testTwoFacetSubjectLevelGammaDynrelRecovers', ...
                'gammascale', 1, 'chains', 2, 'warmup', 750, 'sampling', 750, 'seed', 12345, ...
                'nsub', nsub, 'nocc', nocc, 'ntrl', ntrl, 'design', 'two-facet persons x occasions x trials'));

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma dynrel trt sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'two-facet subject-level dynrel fit non-converged (max R-hat=%.4f); ', ...
                'recovery not validated in this environment.'], maxRhat));

            for f = {'ind_bs','ind_nu','ssinfo'}
                testCase.verifyTrue(isfield(rel.out, f{1}), sprintf( ...
                    'REL.out.%s missing on a gamma analysis-27 run.', f{1}));
            end

            summ = psyrat_dynrel_summary(rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;
            testCase.verifyNotEmpty(T, ...
                'No per-participant table: analysis 27 degraded to analysis 14.');
            testCase.verifyEqual(height(T), nsub);
            testCase.verifyTrue(all(isfinite(T.gen_pt)));

            r = corr_(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel trt sserr recovery] per-participant G corr = %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.60, sprintf( ...
                ['Recovered per-participant coefficients correlate only %.3f with ', ...
                'the planted truth.'], r));

            % SCALE, on the unsaturated implied residual rather than the saturated
            % coefficient - the correlation above is rank-invariant and cannot see a
            % constant multiplicative error in V_w.
            errRatio = localImpliedErrRatio(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel trt sserr recovery] implied error-variance ratio = %.3f\n', errRatio);
            testCase.verifyGreaterThan(errRatio, 0.55, sprintf( ...
                'Implied per-participant error variance is %.2fx the planted value.', errRatio));
            testCase.verifyLessThan(errRatio, 1.80, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- a person-facet leak into V_w would show up here.'], errRatio));

            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3);
        end

    end
end

function [tbl, truth] = localSimulateGammaDynrelSserr(nsub, ntrl, ...
    alpha0, b, lognu0, b_nu, sP, sI, sV, rho)
% One-facet dynamic gamma table. The dimension value z is a BETWEEN-person
% covariate (one per participant, the Rast & Clayson design), carried in a column
% the toolbox standardizes. The person block (u_p, v_p) is correlated; nu_p is
% subject-constant AND z-conditioned via b_nu. Observations are exact continuous
% gamma draws, Y = gamrnd(nu_p/2, 2*mu/nu_p) with E[Y] = mu, Var[Y] = 2*mu^2/nu_p,
% so the continuous log-nu likelihood is the generator's own model.
rng(12345);
z1 = randn(nsub,1);
z2 = randn(nsub,1);
u_p = sP * z1;
v_p = sV * (rho * z1 + sqrt(1 - rho^2) * z2);

% raw (unstandardized) dimension values; psyrat_zscore_subjectlevel standardizes
zraw = randn(nsub,1);
zstd = (zraw - mean(zraw)) ./ std(zraw);

u_t = sI * randn(ntrl,1);

nrow = nsub * ntrl;
ids  = cell(nrow,1);
meas = zeros(nrow,1);
dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    lognu_p = lognu0 + b_nu * zstd(p) + v_p(p);
    nuDraw = exp(lognu_p);                      % continuous df, no rounding/clipping
    for t = 1:ntrl
        r = r + 1;
        ids{r}  = sprintf('S%02d', p);
        dim1(r) = zraw(p);
        logmu = alpha0 + b * zstd(p) + u_p(p) + u_t(t);
        mu = exp(logmu);
        meas(r) = gamrnd(nuDraw / 2, 2 * mu / nuDraw);   % (mu/nu) * chi2(nu), noninteger nu
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id','meas','dim1'});

% ---- closed-form per-participant truth, at each participant's own z ----
% Independent S20 one-facet target (S20_RECOVERY_REFERENCE.md, "Independent
% one-facet target"): with mu_pi = exp(alpha_z + U_p + U_i) and independent
% centered normal log effects of variances vP = sP^2, vI = sI^2,
%   B  = exp(2 alpha_z + vP + vI),
%   P  = B expm1(vP),            PI = B expm1(vP) expm1(vI),
%   C_p = 2 exp(2(alpha_z + u_p) + 2 vI - (lognu_z + v_p)),
%   E_p = PI + C_p,              G_p = P / (P + E_p / n_i).
% PI is a population quantity at the person's z: it is NOT recomputed from
% alpha_z + u_p, and the exp(2 vI) inside C_p is a second-moment factor that
% does not already contain it. Explicit product terms only; no production
% converter is called.
G  = zeros(nsub,1);
Pt = zeros(nsub,1); PIt = zeros(nsub,1); Ct = zeros(nsub,1); Et = zeros(nsub,1);
for p = 1:nsub
    alpha_z = alpha0 + b * zstd(p);
    lognu_z = lognu0 + b_nu * zstd(p);
    B_true      = exp(2*alpha_z + sP^2 + sI^2);
    P_true      = B_true * expm1(sP^2);
    piMean_true = B_true * expm1(sP^2) * expm1(sI^2);
    C_p   = 2 * exp(2*(alpha_z + u_p(p)) + 2*sI^2 - (lognu_z + v_p(p)));
    sigE2 = piMean_true + C_p;
    G(p)  = P_true / (P_true + sigE2 / ntrl);
    Pt(p) = P_true; PIt(p) = piMean_true; Ct(p) = C_p; Et(p) = sigE2;
end
truth.G_true = G;
truth.P_true = Pt; truth.PI_true = PIt; truth.C_p_true = Ct; truth.E_true = Et;
truth.ntrl = ntrl;
truth.u_p = u_p; truth.v_p = v_p; truth.zstd = zstd;
truth.residual_definition = 'E_p = PI(z) + C_p (population mean interaction plus own observation variance); G_p = P/(P + E_p/ntrl)';
end

function [tbl, truth] = localSimulateGammaDynrelSserrTrt(nsub, nocc, ntrl, ...
    alpha0, b, lognu0, b_nu, sP, sO, sT, sPO, sPT, sOT, sV, rho)
% Two-facet (trial + occasion) dynamic gamma table, same construction as the
% one-facet simulator plus the five crossed log-mean facets. Observations are
% exact continuous gamma draws with unrounded nu_p (see the one-facet simulator).
rng(12345);
a1 = randn(nsub,1); a2 = randn(nsub,1);
u_p = sP * a1;
v_p = sV * (rho * a1 + sqrt(1 - rho^2) * a2);

zraw = randn(nsub,1);
zstd = (zraw - mean(zraw)) ./ std(zraw);

u_o  = sO  * randn(nocc,1);
u_t  = sT  * randn(ntrl,1);
u_po = sPO * randn(nsub,nocc);
u_pt = sPT * randn(nsub,ntrl);
u_ot = sOT * randn(ntrl,nocc);

nrow = nsub * nocc * ntrl;
ids  = cell(nrow,1);
time = zeros(nrow,1);
meas = zeros(nrow,1);
dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    lognu_p = lognu0 + b_nu * zstd(p) + v_p(p);
    nuDraw = exp(lognu_p);                      % continuous df, no rounding/clipping
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            dim1(r) = zraw(p);
            logmu = alpha0 + b * zstd(p) + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p,o) + u_pt(p,t) + u_ot(t,o);
            mu = exp(logmu);
            meas(r) = gamrnd(nuDraw / 2, 2 * mu / nuDraw);   % (mu/nu) * chi2(nu), noninteger nu
        end
    end
end
tbl = table(ids, meas, time, dim1, 'VariableNames', {'id','meas','time','dim1'});

% ---- closed-form per-participant truth (n_i = ntrl, n_o = 1: the summary's
% documented D-study default; NOT the observed number of occasions) ----
% Q8 two-facet target (SCIENTIFIC_CONTRACT.md; S20_RECOVERY_REFERENCE.md,
% "Separate two-facet Q8 recovery reference"): P, PI, PO and R = K_POI are the
% population mean components at the participant's own z from the support-union
% products (localSupportUnionComponents), and
%   W   = sO^2 + sT^2 + sPO^2 + sPT^2 + sOT^2,
%   C_p = 2 exp(2(alpha_z + u_p) + 2 W - (lognu_z + v_p)),
%   E_p = R + C_p,
%   G_p = P / (P + PI/n_i + PO/n_o + E_p/(n_i n_o))      (CES, relative error).
% PI keeps its own 1/n_i slot; R is the highest-order remainder and is added
% ONCE to the residual. No production converter or reliability kernel is used.
W  = sO^2 + sT^2 + sPO^2 + sPT^2 + sOT^2;
n_i = ntrl; n_o = 1;
G  = zeros(nsub,1);
Pt = zeros(nsub,1); PIt = zeros(nsub,1); POt = zeros(nsub,1);
Rt = zeros(nsub,1); Ct = zeros(nsub,1); Et = zeros(nsub,1);
for p = 1:nsub
    alpha_z = alpha0 + b * zstd(p);
    lognu_z = lognu0 + b_nu * zstd(p);
    K = localSupportUnionComponents(alpha_z, sP^2, sO^2, sT^2, sPO^2, sPT^2, sOT^2);
    C_p = 2 * exp(2*(alpha_z + u_p(p)) + 2*W - (lognu_z + v_p(p)));
    E_p = K.R + C_p;
    G(p) = K.P / (K.P + K.PI/n_i + K.PO/n_o + E_p/(n_i*n_o));
    Pt(p) = K.P; PIt(p) = K.PI; POt(p) = K.PO; Rt(p) = K.R; Ct(p) = C_p; Et(p) = E_p;
end
truth.G_true = G;
truth.P_true = Pt; truth.PI_true = PIt; truth.PO_true = POt; truth.R_true = Rt;
truth.C_p_true = Ct; truth.E_true = Et;
truth.ntrl = ntrl; truth.nocc_dstudy = n_o;
truth.u_p = u_p; truth.v_p = v_p; truth.zstd = zstd;
truth.residual_definition = 'E_p = R(z) + C_p (population highest-order remainder plus own observation variance); G_p = P/(P + PI/n_i + PO/n_o + E_p/(n_i n_o)), n_i = ntrl, n_o = 1';
end

function K = localSupportUnionComponents(alpha_z, vP, vO, vI, vPO, vPI, vOI)
% Independent support-union products (SCIENTIFIC_CONTRACT.md, "Equations fixed
% before numerical edits"). With log mean alpha(z) + U_P + U_O + U_I + U_PO +
% U_PI + U_OI, independent centered normal factor families of variances v_A,
% V = sum(v_A), B = exp(2 alpha + V) and q_A = expm1(v_A), each observed-scale
% mean component is
%   K_S = B * sum_{T : union of supports of A in T equals S} prod_{A in T} q_A
% over all nonempty subsets T of the six factors, and R = K_{POI}. The 63
% subsets are enumerated explicitly rather than simplified algebraically, so a
% simplification slip cannot hide; no production routine is called.
v = [vP vO vI vPO vPI vOI];
q = expm1(v);
B = exp(2*alpha_z + sum(v));
% support of each factor over the three facets (P, O, I)
supp = logical([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]);
acc = zeros(1,8);   % indexed by 1 + P + 2*O + 4*I of the union mask
for m = 1:63
    inT = logical(bitget(m, 1:6));
    u = any(supp(inT,:), 1);
    idx = 1 + u(1) + 2*u(2) + 4*u(3);
    acc(idx) = acc(idx) + prod(q(inT));
end
% every nonempty subset has exactly one union, so the components sum to
% prod(1+q) - 1; a mismatch means the enumeration itself is wrong
assert(abs(sum(acc) - (prod(1+q) - 1)) <= 1e-12 * max(1, prod(1+q)), ...
    'localSupportUnionComponents: subset enumeration does not sum to prod(1+q)-1.');
K.B  = B;
K.P  = B * acc(2);   % union {P}
K.O  = B * acc(3);   % union {O}
K.PO = B * acc(4);   % union {P,O}
K.I  = B * acc(5);   % union {I}
K.PI = B * acc(6);   % union {P,I}
K.OI = B * acc(7);   % union {O,I}
K.R  = B * acc(8);   % union {P,O,I} = highest-order remainder
end

function ratio = localImpliedErrRatio(Ghat, Gtrue)
% Ratio of implied to planted per-participant error variance, recovered by
% inverting the coefficient. From G = bp/(bp + err/n), err/n = bp*(1-G)/G, so
% (1-G)/G is proportional to the error variance at fixed bp and n. Taking the mean
% ratio across participants gives a scale check that does NOT saturate as G -> 1,
% unlike a bound on G itself.
oddsHat  = (1 - Ghat(:))  ./ Ghat(:);
oddsTrue = (1 - Gtrue(:)) ./ Gtrue(:);
ratio = mean(oddsHat) / mean(oddsTrue);
end

function r = corr_(a, b)
% Pearson correlation without the Statistics Toolbox.
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function rel = localRunDynrelSserrGamma(testCase, dataTbl, tag)
% dynrel=2 (a Dimension 1 column) + sserrvar=2 routes to analysis 26.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    ... %pin log-nu: this design is movable, and since rulings 33-36 the
    ... %engine default is location-scale. Without this the test would
    ... %silently change estimand and check log-nu targets against an
    ... %LS fit. The location-scale twin has its own *LsRecovery test.
    'gammascale', 1, ...
    'dynrel', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function rel = localRunDynrelSserrTrtGamma(testCase, dataTbl, tag)
% dynrel=2 + a time facet + sserrvar=2 routes to analysis 27.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    ... %pin log-nu: this design is movable, and since rulings 33-36 the
    ... %engine default is location-scale. Without this the test would
    ... %silently change estimand and check log-nu targets against an
    ... %LS fit. The location-scale twin has its own *LsRecovery test.
    'gammascale', 1, ...
    'dynrel', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
% Remove the per-run output directories and any stray CmdStan build/output files.
for tag = {'recov_gamma_dynrelss','recov_gamma_dynrelsstrt'}
    p = fullfile(baseDir, tag{1});
    if isfolder(p)
        try
            rmdir(p, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                ['Could not remove recovery artifact directory ''%s''. ', ...
                'How to fix: close any process using it and delete it manually.'], p);
        end
    end
end

pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
for i = 1:numel(pats)
    files = dir(fullfile(baseDir, pats{i}));
    for j = 1:numel(files)
        if ~files(j).isdir
            delete(fullfile(baseDir, files(j).name));
        end
    end
end

dirs = dir(fullfile(baseDir, 'cmdstan_*'));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
for i = 1:numel(dirs)
    try
        rmdir(fullfile(baseDir, dirs(i).name), 's');
    catch
        % best-effort cleanup
    end
end
end
