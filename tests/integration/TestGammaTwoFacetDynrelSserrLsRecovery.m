classdef TestGammaTwoFacetDynrelSserrLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE gamma SUBJECT-LEVEL
    % two-facet dynamic/conditional reliability design: analysis 27
    % (ic_dynrel_sserrvar_trt) with family='gamma', gammascale=2, dynrel=2,
    % sserrvar=2. This is the "Dynamic 27 / LS" cell of the bounded 2C route
    % matrix (the 2C validation-plan review, a local development record),
    % which had no non-difference recovery test before this class.
    % TestGammaDynrelTrtLsRecovery fits analysis 14 only, and
    % TestGammaDynrelLsRecovery supplies sserrvar=2 only for one facet (26).
    %
    % GENERATIVE MODEL (independently specified; no production converter or
    % reliability kernel is used to build expected values). For person p,
    % occasion o, trial i, with one standardized between-person covariate z_p:
    %   log mu_poi   = alpha0 + b z_p + U_p + U_o + U_i + U_po + U_pi + U_oi
    %   log sigma_p  = ell0 + b_sigma z_p + w_p
    %   y_poi ~ Gamma(shape = (mu/sigma_p)^2, scale = sigma_p^2/mu)
    % so, conditional on ALL planted effects (the six log-mean effects and w_p),
    % E[y] = mu and Var[y | all effects] = sigma_p^2 (an exact continuous gamma
    % draw; (mu/nu) chi2(nu) with nu = 2 (mu/sigma)^2, noninteger). Var[y | p]
    % alone is larger, since it also carries the facet effects. The six log-mean
    % factor families are independent centered normals with their own grouping
    % (U_o shared across persons and trials, U_pi by person x trial, ...). The
    % person block (U_p, w_p) is correlated at rhoPW.
    %
    % PLANTED TRUTH (SCIENTIFIC_CONTRACT.md "Equations fixed before numerical
    % edits" and the Q8 target). At the participant's own z, with v_A the six
    % log-mean variances, V = sum(v_A), B = exp(2 alpha(z) + V), q_A = expm1(v_A):
    %   K_S = B sum_{T: union(supp(A), A in T) = S} prod_{A in T} q_A
    % gives the population components P, O, I, PO, PI, OI and the remainder
    % R = K_POI (explicit enumeration of all 63 subsets). Location-scale own
    % observation variance C_p = exp(2 (ell(z) + w_p)), E_p = R + C_p. Under the
    % CES allocation at the subject's trial count n_i and the default n_o = 1:
    %   G_p = P / (P + PI/n_i + PO/n_o + E_p/(n_i n_o))
    %   D_p = P / (P + PI/n_i + PO/n_o + E_p/(n_i n_o) + I/n_i + O/n_o + OI/(n_i n_o))
    %
    % PREDECLARED ASSERTIONS (fixed before any root execution; not to be tuned
    % to the result):
    %   - route: rel.analysis 'ic_dynrel_sserrvar_trt', family gamma,
    %     gammascale 2, engine cmdstan; summary routed as location-scale gamma;
    %     REL.out.ind_sdlog / ssinfo present with draws x NSUB / NSUB rows;
    %   - every per-participant coefficient (gen and dep) finite; NSUB rows;
    %   - identity/ordering: table ids match the simulated labels, id2 is a
    %     permutation of 1:NSUB, and each row's z1 equals the planted z of the
    %     participant it names (AbsTol 1e-6);
    %   - Pearson corr(G_hat, G_true) > .60;
    %   - median over participants of |E_hat - E_true| / E_true < .80, where
    %     E_hat is a DIAGNOSTIC transform: the residual implied by inverting
    %     the table's point coefficient gen_pt with the planted P, PI, PO (all
    %     estimation error attributed to the residual). It is not the posterior
    %     mean of the residual and is not read from any posterior residual draw;
    %   - R_true > 0 for every participant (the fixture discriminates the
    %     pooled target from the observation-only one).
    % Coverage of G_true / D_true by the per-participant 95% intervals and
    % coverage of the planted ell0, b_sigma, b, s_P, s_W by their posterior
    % intervals are RECORDED (printed and returned in the capture), not
    % asserted: root adjudicates inference coverage rather than this test
    % inventing a pass rule. Observation-only vs pooled residual discrimination
    % is likewise recorded diagnostically.
    %
    % GATES. Dependency setup, gamrnd availability, a nonempty convergence
    % table and psyrat_checkconv's converged flag are assumptions, so a missing
    % CmdStan / Statistics Toolbox or a non-converged fit reports Incomplete,
    % never a false PASS. The available diagnostics are printed verbatim; a
    % pass of this early gate is not final-inference clearance. When
    % PSYRAT_2C_RECOVERY_CAPTURE_DIR is set, REL/truth/data and the raw CmdStan
    % output are captured through PsyRAT2CRecoveryCapture BEFORE the
    % convergence assumption, and the class teardown skips its scratch cleanup
    % so the raw output also survives a failed capture or test.
    %
    % Sampler settings mirror the sibling recovery classes exactly: 2 chains,
    % 750 warmup, 750 sampling, seed 12345. No new sampler settings.

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
                fprintf('[gamma LS dynrel trt sserr recovery] %s set; scratch cleanup skipped.\n', ...
                    PsyRAT2CRecoveryCapture.EnvVar);
                return;
            end
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)

        function testTwoFacetSubjectLevelLsDynrelRecovers(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous gamma observations.');

            % ---- Planted truth (log expected-score / log residual-SD scales) ----
            c = localTruthConstants();

            [dataTbl, truth] = localSimulateGammaDynrelSserrTrtLs(c);
            nsub = c.nsub;

            % the fixture must discriminate R + C_p from C_p alone
            testCase.verifyTrue(all(truth.R_true > 0), ...
                'Planted highest-order remainder R is not positive; the fixture cannot discriminate the pooled target.');

            rel = localRunDynrelSserrTrtGammaLs(testCase, dataTbl, 'recov_gamma_dynrelsstrt_ls');

            % ---- Route ----
            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar_trt');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Opt-in evidence capture, BEFORE any convergence assumption ----
            PsyRAT2CRecoveryCapture.capture('dynrel27_ls_sserr', rel, truth, dataTbl, ...
                fullfile(localCmdStanOutputDir(testCase.projectRoot()), 'recov_gamma_dynrelsstrt_ls'), ...
                'meta', struct('test', 'TestGammaTwoFacetDynrelSserrLsRecovery/testTwoFacetSubjectLevelLsDynrelRecovers', ...
                'gammascale', 2, 'chains', 2, 'warmup', 750, 'sampling', 750, 'seed', 12345, ...
                'rng_seed', c.rngSeed, 'nsub', c.nsub, 'nocc', c.nocc, 'ntrl', c.ntrl, ...
                'design', 'two-facet persons x occasions x trials, location-scale'));

            % ---- Convergence gate (assumption) with the raw diagnostics printed ----
            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma LS dynrel trt sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            fprintf('[gamma LS dynrel trt sserr recovery] REL.out.conv fields available: %s\n', ...
                strjoin(fieldnames(rel.out.conv)', ', '));
            localPrintConvTable(rel.out.conv.data);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma LS ', ...
                'two-facet subject-level dynrel fit non-converged (max R-hat=%.4f >= 1.1 / ', ...
                'n_eff gate); recovery not validated in this environment.'], maxRhat));

            % ---- Subject-level slots exist and are sized for this design ----
            for f = {'ind_sdlog', 'ssinfo'}
                testCase.verifyTrue(isfield(rel.out, f{1}) && ~isempty(rel.out.(f{1})), sprintf( ...
                    'REL.out.%s missing/empty on a gamma LS analysis-27 run.', f{1}));
            end
            indsd = cell2mat(rel.out.ind_sdlog{1});   % draws x NSUB (person log-sigma w_p)
            ndraws = size(indsd, 1);
            testCase.verifySize(indsd, [ndraws nsub]);
            testCase.verifyEqual(height(rel.out.ssinfo{1}), nsub);

            % ---- Summary routes as location-scale gamma and is per-participant ----
            % reltype 3 (CES) and nocc 1 are the documented summary defaults;
            % passed explicitly so the reference count is pinned, not inherited.
            summ = psyrat_dynrel_summary(rel, 'CI', 0.95, 'reltype', 3, 'nocc', 1);
            testCase.verifyTrue(isfield(summ, 'isgamma_ls') && logical(summ.isgamma_ls), ...
                'The summary did not route this run as location-scale gamma.');
            T = summ.strata(1).ssrel_table;
            testCase.assertNotEmpty(T, ...
                'No per-participant table: analysis 27 degraded to analysis 14.');
            testCase.verifyEqual(height(T), nsub);
            testCase.assertTrue(all(ismember({'id', 'id2', 'z1', 'gen_pt', 'gen_ll', 'gen_ul', ...
                'dep_pt', 'dep_ll', 'dep_ul'}, T.Properties.VariableNames)), ...
                'Per-participant table lacks an expected id/id2/z1/gen/dep column.');
            testCase.verifyTrue(all(isfinite(T.gen_pt)) && all(isfinite(T.gen_ll)) && all(isfinite(T.gen_ul)), ...
                'Non-finite per-participant generalizability coefficient(s).');
            testCase.verifyTrue(all(isfinite(T.dep_pt)) && all(isfinite(T.dep_ll)) && all(isfinite(T.dep_ul)), ...
                'Non-finite per-participant dependability coefficient(s).');
            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3, ...
                'Per-participant coefficients are constant; the read-out is not conditional.');

            % ---- Identity / ordering: map each table row to its planted participant ----
            testCase.verifyTrue(isequal(sort(T.id2(:)'), 1:nsub), ...
                'ssinfo.id2 is not a permutation of 1:NSUB.');
            rowLabels = string(T.id(:));
            [tf, pIdx] = ismember(rowLabels, string(truth.ids(:)));
            testCase.assertTrue(all(tf), ...
                'Per-participant table ids do not all match the simulated participant labels.');
            testCase.verifyEqual(sort(pIdx(:)'), 1:nsub, ...
                'Per-participant table does not name each simulated participant exactly once.');
            testCase.verifyEqual(T.z1(:), truth.zstd(pIdx), 'AbsTol', 1e-6, ...
                ['Per-participant z1 differs from the planted standardized z of the named ', ...
                'participant; the toolbox standardization or the id alignment is off.']);
            if ismember('trls', T.Properties.VariableNames)
                fprintf('[gamma LS dynrel trt sserr recovery] subject trial counts (min/median/max) = %g / %g / %g (planted %d)\n', ...
                    min(T.trls), median(T.trls), max(T.trls), c.ntrl);
            end

            Gtrue = truth.G_true(pIdx);
            Dtrue = truth.D_true(pIdx);
            Etrue = truth.E_true(pIdx);
            Ctrue = truth.C_p_true(pIdx);
            Ptrue = truth.P_true(pIdx);
            PItrue = truth.PI_true(pIdx);
            POtrue = truth.PO_true(pIdx);

            % ---- Recovery: ordering ----
            rG = corr_(T.gen_pt, Gtrue);
            rD = corr_(T.dep_pt, Dtrue);
            fprintf('[gamma LS dynrel trt sserr recovery] per-participant corr: G %.3f, D %.3f\n', rG, rD);
            testCase.verifyGreaterThan(rG, 0.60, sprintf( ...
                ['Recovered per-participant generalizability coefficients correlate only %.3f ', ...
                'with the planted truth. Either the conditional conversion is wrong or ', ...
                'ind_sdlog is misaligned against id2.'], rG));

            % ---- Recovery: pooled residual scale (predeclared) ----
            % Invert G = P/(P + PI/n_i + PO/n_o + E/(n_i n_o)) at the planted P, PI,
            % PO and the reference counts to get the residual each POINT
            % coefficient (gen_pt) implies; compare with the planted pooled
            % E = R + C_p. This is a diagnostic transform of the point
            % coefficient, not the posterior mean of the residual.
            n_i = c.ntrl; n_o = 1;
            Ehat = (n_i * n_o) .* (Ptrue .* (1 - T.gen_pt) ./ T.gen_pt - PItrue ./ n_i - POtrue ./ n_o);
            relDevPooled = abs(Ehat - Etrue) ./ Etrue;
            medRelDev = median(relDevPooled);
            fprintf('[gamma LS dynrel trt sserr recovery] implied pooled residual: median rel. abs. deviation = %.3f (min %.3f, max %.3f)\n', ...
                medRelDev, min(relDevPooled), max(relDevPooled));
            testCase.verifyLessThan(medRelDev, 0.80, sprintf( ...
                ['Median relative absolute deviation of the implied pooled residual from ', ...
                'the planted R + C_p is %.3f (predeclared bound .80).'], medRelDev));
            % RECORD (2026-09-21, owner rulings 2-4; the findings ledger's
            % LS27 entry; kept BELOW the assertion so its line number stays :215).
            % The assertion above FAILED on its only recorded run: 1.798368833778583
            % against .80 on a converged fit (2C acceptance, 2026-09-13). One
            % prospective dataset in the same regime replicated it (A: converged,
            % C0 1.6003); a second in a different regime (R/E .86) showed the same
            % component inflation, would have passed C0 (.3578) and did not
            % converge (R-hat 1.1351, ESS 12); and a no-fit re-evaluation of every
            % 2C gamma two-facet capture showed the same mechanism: with three or four
            % occasion levels the raw-scale components are exponentials of squared
            % log-scale SDs whose posteriors are right-skewed, so their posterior
            % MEANS (what production reports) are tail-dominated, and the point
            % coefficient carries a small downward bias that the inversion above
            % amplifies. Not a formula error (wherever production was compared with
            % an independent enumeration it agreed to 2e-12 or better). Ruling: the
            % assertion, its bound and this class's
            % assumptions stay exactly as predeclared and the FAIL stays on record
            % as the adjudicated expected outcome of the design
            % (RELEASE_VERIFICATION_CHECKLIST.md section 6); a convergence-filtered
            % Incomplete on this class is recorded the same way, never counted as a
            % pass. The gamma two-facet designs (14, 25, 27; location-scale and
            % log-nu) are RESTRICTED and disclosed at selection, on view and in
            % every export header (psyrat_gamma_twofacet_restriction). No threshold
            % change, no relabel, no skip. The one-facet gamma designs (11, 26) do
            % not show the pattern and are outside the label.

            % ---- Diagnostic (recorded, not asserted): pooled vs observation-only ----
            relDevObsOnly = abs(Ehat - Ctrue) ./ Ctrue;
            fprintf(['[gamma LS dynrel trt sserr recovery] discrimination: R_true min/median/max = %.4g / %.4g / %.4g; ', ...
                'R/(R+C_p) min %.3f median %.3f max %.3f; median rel. dev. of implied residual ', ...
                'vs pooled %.3f vs observation-only %.3f\n'], ...
                min(truth.R_true), median(truth.R_true), max(truth.R_true), ...
                min(truth.R_true(pIdx) ./ Etrue), median(truth.R_true(pIdx) ./ Etrue), ...
                max(truth.R_true(pIdx) ./ Etrue), medRelDev, median(relDevObsOnly));

            % ---- Coverage (recorded, not asserted) ----
            covG = mean(Gtrue >= T.gen_ll & Gtrue <= T.gen_ul);
            covD = mean(Dtrue >= T.dep_ll & Dtrue <= T.dep_ul);
            fprintf(['[gamma LS dynrel trt sserr recovery] per-participant 95%% interval coverage of ', ...
                'the planted coefficient: G %.3f (%d/%d), D %.3f (%d/%d); median interval width G %.4f, D %.4f\n'], ...
                covG, round(covG * nsub), nsub, covD, round(covD * nsub), nsub, ...
                median(T.gen_ul - T.gen_ll), median(T.dep_ul - T.dep_ll));

            % ---- Parameter recovery (recorded, not asserted) ----
            % Analysis 27 LS keeps the raw log-sigma intercept/effects; read the raw
            % fields the LS14 sibling reads and report interval coverage. Absent
            % fields are reported as unavailable rather than assumed.
            localReportParam('ell0 (pop_sdlog)', localField(rel.out, 'pop_sdlog', 1), c.ell0);
            localReportParam('b_sigma',          localField(rel.out, 'b_sigma', 1), c.b_sigma);
            localReportParam('b (mean slope)',   localField(rel.out, 'b', 1), c.b);
            gro = localField(rel.out, 'gro_sds', 1);
            if ~isempty(gro) && size(gro, 2) >= 2
                localReportParam('s_P (gro_sds 1)', gro(:,1), c.sP);
                localReportParam('s_W (gro_sds 2)', gro(:,2), c.sW);
            else
                fprintf('[gamma LS dynrel trt sserr recovery] gro_sds unavailable or not draws x 2; s_P / s_W coverage not recorded\n');
            end
            % planted w_p versus posterior-median per-person log-sigma offset
            wHat = median(indsd, 1)';                 % NSUB x 1, column c = id2
            wHatByRow = wHat(T.id2(:));
            fprintf('[gamma LS dynrel trt sserr recovery] corr(planted w_p, posterior-median ind_sdlog) = %.3f\n', ...
                corr_(truth.w_p(pIdx), wHatByRow));
        end

    end
end

% =========================================================================
function c = localTruthConstants()
% Planted constants (independently specified for this cell; brief of 2026-09-13).
c = struct();
c.alpha0  = log(5.5);   c.b       =  0.20;   % log-mean intercept + dimension slope
c.ell0    = log(2);     c.b_sigma = -0.30;   % log residual-SD intercept + dimension slope
c.sP  = 0.35;   % person
c.sO  = 0.22;   % occasion
c.sT  = 0.15;   % trial
c.sPO = 0.12;   % person x occasion
c.sPT = 0.10;   % person x trial
c.sOT = 0.08;   % occasion x trial
c.sW  = 0.30;   % person log residual-SD SD
c.rhoPW = 0.35; % corr(person log-mean, person log residual-SD)
c.nsub = 40; c.nocc = 4; c.ntrl = 12;   % 1920 rows
c.rngSeed = 71027;
end

% =========================================================================
function [tbl, truth] = localSimulateGammaDynrelSserrTrtLs(c)
% Balanced crossed Persons x Occasions x Trials long table from the generative
% model in the class header. Returns the per-participant planted truth at each
% participant's own standardized z.
rng(c.rngSeed);
nsub = c.nsub; nocc = c.nocc; ntrl = c.ntrl;

% person block: (u_p, w_p) correlated at rhoPW through a shared standard normal
a1 = randn(nsub, 1); a2 = randn(nsub, 1);
u_p = c.sP * a1;
w_p = c.sW * (c.rhoPW * a1 + sqrt(1 - c.rhoPW^2) * a2);

% between-person covariate: raw values go in the table; the toolbox standardizes
% once across participants (center, sample SD), which is what zstd reproduces
zraw = randn(nsub, 1);
zstd = (zraw - mean(zraw)) ./ std(zraw);

% independent log-mean factor families at their own grouping
u_o  = c.sO  * randn(nocc, 1);          % occasion: shared across persons and trials
u_t  = c.sT  * randn(ntrl, 1);          % trial:    shared across persons and occasions
u_po = c.sPO * randn(nsub, nocc);       % person x occasion
u_pt = c.sPT * randn(nsub, ntrl);       % person x trial
u_ot = c.sOT * randn(nocc, ntrl);       % occasion x trial

nrow = nsub * nocc * ntrl;
ids  = cell(nrow, 1);
time = zeros(nrow, 1);
meas = zeros(nrow, 1);
dim1 = zeros(nrow, 1);
labels = cell(nsub, 1);
r = 0;
for p = 1:nsub
    labels{p} = sprintf('S%02d', p);
    sigma_p = exp(c.ell0 + c.b_sigma * zstd(p) + w_p(p));   % own residual SD at own z
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = labels{p};
            time(r) = o;
            dim1(r) = zraw(p);
            mu = exp(c.alpha0 + c.b * zstd(p) + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p, o) + u_pt(p, t) + u_ot(o, t));
            % Gamma(shape, scale): mean = shape*scale = mu, var = shape*scale^2 = sigma_p^2
            shape = (mu / sigma_p)^2;
            scale = sigma_p^2 / mu;
            meas(r) = gamrnd(shape, scale);
        end
    end
end
tbl = table(ids, meas, time, dim1, 'VariableNames', {'id', 'meas', 'time', 'dim1'});

% ---- per-participant planted truth at own z; n_i = ntrl, n_o = 1 (CES) ----
n_i = ntrl; n_o = 1;
G = zeros(nsub, 1); D = zeros(nsub, 1);
Pt = zeros(nsub, 1); Ot = zeros(nsub, 1); It = zeros(nsub, 1);
POt = zeros(nsub, 1); PIt = zeros(nsub, 1); OIt = zeros(nsub, 1);
Rt = zeros(nsub, 1); Ct = zeros(nsub, 1); Et = zeros(nsub, 1);
for p = 1:nsub
    alpha_z = c.alpha0 + c.b * zstd(p);
    ell_z   = c.ell0 + c.b_sigma * zstd(p);
    K = localSupportUnionComponents(alpha_z, c.sP^2, c.sO^2, c.sT^2, c.sPO^2, c.sPT^2, c.sOT^2);
    C_p = exp(2 * (ell_z + w_p(p)));          % location-scale own observation variance
    E_p = K.R + C_p;                          % population remainder + own observation variance
    relErr = K.PI / n_i + K.PO / n_o + E_p / (n_i * n_o);
    absAdd = K.I / n_i + K.O / n_o + K.OI / (n_i * n_o);
    G(p) = K.P / (K.P + relErr);
    D(p) = K.P / (K.P + relErr + absAdd);
    Pt(p) = K.P; Ot(p) = K.O; It(p) = K.I;
    POt(p) = K.PO; PIt(p) = K.PI; OIt(p) = K.OI;
    Rt(p) = K.R; Ct(p) = C_p; Et(p) = E_p;
end
truth.G_true = G; truth.D_true = D;
truth.P_true = Pt; truth.O_true = Ot; truth.I_true = It;
truth.PO_true = POt; truth.PI_true = PIt; truth.OI_true = OIt;
truth.R_true = Rt; truth.C_p_true = Ct; truth.E_true = Et;
truth.u_p = u_p; truth.w_p = w_p; truth.zstd = zstd; truth.zraw = zraw;
truth.ids = labels;
truth.ntrl = ntrl; truth.nocc_data = nocc; truth.nocc_dstudy = n_o;
truth.constants = c;
truth.residual_definition = ['E_p = R(z) + C_p, C_p = exp(2(ell0 + b_sigma z + w_p)); ', ...
    'CES: G_p = P/(P + PI/n_i + PO/n_o + E_p/(n_i n_o)), ', ...
    'D_p = P/(P + PI/n_i + PO/n_o + E_p/(n_i n_o) + I/n_i + O/n_o + OI/(n_i n_o)), n_i = ntrl, n_o = 1'];
end

% =========================================================================
function K = localSupportUnionComponents(alpha_z, vP, vO, vI, vPO, vPI, vOI)
% Independent support-union products (SCIENTIFIC_CONTRACT.md, "Equations fixed
% before numerical edits"). With log mean alpha(z) + U_P + U_O + U_I + U_PO +
% U_PI + U_OI, independent centered normal factor families of variances v_A,
% V = sum(v_A), B = exp(2 alpha + V) and q_A = expm1(v_A), each observed-scale
% mean component is
%   K_S = B * sum_{T : union of supports of A in T equals S} prod_{A in T} q_A
% over all nonempty subsets T of the six factors, and R = K_{POI}. The 63
% subsets are enumerated explicitly rather than simplified algebraically, so a
% simplification slip cannot hide; no production routine is called. (Duplicated
% in TestGammaDynrelSserrRecovery so each class stays self-contained.)
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

% =========================================================================
function r = corr_(a, b)
% Pearson correlation without the Statistics Toolbox.
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / sqrt((a' * a) * (b' * b));
end

% =========================================================================
function x = localField(out, name, stratum)
% Draws of a REL.out field for one stratum, or [] when absent. Handles the two
% storage conventions seen in the sibling LS tests: a cell per stratum holding
% draws x k (b, b_sigma, gro_sds) and a draws x nstrata numeric (pop_sdlog).
x = [];
if ~isstruct(out) || ~isfield(out, name) || isempty(out.(name))
    return;
end
v = out.(name);
if iscell(v)
    if numel(v) >= stratum && ~isempty(v{stratum})
        x = v{stratum};
        if iscell(x), x = cell2mat(x); end
    end
elseif isnumeric(v)
    if size(v, 2) >= stratum
        x = v(:, stratum);
    end
end
end

% =========================================================================
function localReportParam(label, draws, trueVal)
% Print posterior median and a sorted-draw 95% interval with the planted value
% and whether it is covered. Recorded only; no assertion. Base MATLAB only.
if isempty(draws)
    fprintf('[gamma LS dynrel trt sserr recovery] %s: field unavailable; coverage not recorded (true %.4f)\n', ...
        label, trueVal);
    return;
end
d = sort(draws(:, 1));
n = numel(d);
ll = d(max(1, floor(0.025 * n)));
ul = d(min(n, ceil(0.975 * n)));
md = median(d);
fprintf('[gamma LS dynrel trt sserr recovery] %s: median %.4f, 95%% [%.4f, %.4f], true %.4f, covered=%d\n', ...
    label, md, ll, ul, trueVal, trueVal >= ll && trueVal <= ul);
end

% =========================================================================
function localPrintConvTable(convData)
% Print the header names and per-column extrema of REL.out.conv.data, using
% the same shape convention psyrat_conv_maxrhat / psyrat_checkconv accept: a
% 1xN cell of per-fit cell tables, or a single Mx3-or-4 cell table
% {name, n_eff, r_hat, ...} whose row 1 is the header and rows 2:end are data.
% Anything else is reported by class/size only; the raw field is retained in
% the capture regardless. Diagnostic only; never throws.
pre = '[gamma LS dynrel trt sserr recovery]';
try
    if ~iscell(convData) || isempty(convData)
        fprintf('%s conv data is %s (size %s); raw field retained, no extrema printed\n', ...
            pre, class(convData), mat2str(size(convData)));
        return;
    end
    if size(convData, 1) == 1 && all(cellfun(@iscell, convData))
        tables = convData;            % 1xN container of per-fit tables
    else
        tables = {convData};          % a single table
    end
    for k = 1:numel(tables)
        tbl = tables{k};
        if ~iscell(tbl) || size(tbl, 1) < 2 || size(tbl, 2) < 3
            fprintf('%s conv table %d: empty/header-only/malformed (size %s); raw field retained\n', ...
                pre, k, mat2str(size(tbl)));
            continue;
        end
        hdr = cellfun(@(x) char(string(x)), tbl(1, :), 'UniformOutput', false);
        fprintf('%s conv table %d: %d data rows; columns: %s\n', ...
            pre, k, size(tbl, 1) - 1, strjoin(hdr, ', '));
        for j = 2:size(tbl, 2)
            col = tbl(2:end, j);
            isnum = cellfun(@(v) isnumeric(v) && isscalar(v), col);
            vals = cell2mat(col(isnum));
            if isempty(vals)
                fprintf('%s   %s: no numeric entries\n', pre, hdr{j});
            else
                fprintf('%s   %s: min %.6g, max %.6g, n finite %d, n NaN %d, n Inf %d (of %d rows)\n', ...
                    pre, hdr{j}, min(vals, [], 'omitnan'), max(vals, [], 'omitnan'), ...
                    sum(isfinite(vals)), sum(isnan(vals)), sum(isinf(vals)), numel(col));
            end
        end
    end
catch errPrint
    fprintf('%s conv table print failed: %s (raw field retained in the capture)\n', pre, errPrint.message);
end
end

% =========================================================================
function rel = localRunDynrelSserrTrtGammaLs(testCase, dataTbl, tag)
% dynrel=2 (a Dimension 1 column) + a time facet + sserrvar=2 + gammascale=2
% routes to analysis 27 under the location-scale scale submodel. Same sampler
% settings as every sibling recovery class.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
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

% =========================================================================
function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

% =========================================================================
function localCleanupArtifacts(baseDir)
% Remove ONLY this class's per-run output directory and the stray CmdStan
% build/output files the wrapper leaves at the artifact root (the same stray
% patterns the sibling classes clean). Captured evidence lives under
% PSYRAT_2C_RECOVERY_CAPTURE_DIR and is never touched here.
p = fullfile(baseDir, 'recov_gamma_dynrelsstrt_ls');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            ['Could not remove recovery artifact directory ''%s''. ', ...
            'How to fix: close any process using it and delete it manually.'], p);
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
end
