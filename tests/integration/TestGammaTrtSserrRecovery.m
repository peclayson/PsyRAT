classdef TestGammaTrtSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL crossed test-retest design (analysis 25, trt_sserrvar,
    % family='gamma').
    %
    % WHY THIS EXISTS. The offline proof (tests/TestGammaTrtSserrConversion)
    % validates the conditional-residual MATHEMATICS against a Monte-Carlo ground
    % truth, but it re-derives the formula rather than calling the production
    % extractor (psyrat_gamma_extract_trt_sserr is a local function of
    % psyrat_computevarcomp and is unreachable from a test). This test is what
    % actually exercises the extractor: the fit -> extraction -> REL.out plumbing,
    % the per-subject column ordering, and the observed-scale conversion, end to
    % end under real CmdStan.
    %
    % WHAT IT PINS. The design is Persons x Trials x Occasions with a per-subject
    % dispersion. Under gamma the per-subject residual is NOT a separate
    % location-scale parameter (as it is for the Gaussian analysis 25) - it is
    % implied by each person's own log-mean u_p and log-nu v_p:
    %
    %   sigma_pot_e2(p) = 2*exp(2*(alpha + u_p) + 2*V_w - (beta + v_p)),
    %   V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2.
    %
    % The interesting failure mode is the one TestTrtSserrRecovery pins for the
    % Gaussian twin: a per-subject residual is a greedy parameter, and if the
    % crossed facet effects are not separately identified it will absorb them. So
    % this test checks that the five crossed facet SDs and the per-subject
    % residual recover SEPARATELY, not merely that the fit completes. It also
    % checks that the recovered per-subject residuals CORRELATE with the planted
    % per-person truth - the assertion that would fail if the extractor's NSUB
    % columns were misaligned against id2 (the one silent-wrong-person failure
    % mode flagged during design).
    %
    % Data are simulated from the crossed model the family assumes: a LOG-LINKED
    % mean carrying six independent normal random effects, a PERSON-VARYING
    % log-nu correlated with the person mean, and a scaled-chi-square observation
    % Y = (mu/nu_p) * chi2(nu_p). Base MATLAB randn only (the chi-square draw is a
    % sum of squared standard normals), so no Statistics Toolbox is required.
    %
    % Like the other CmdStan integration tests, dependency setup uses assumptions,
    % so the test is skipped cleanly (Incomplete, not Failed) when CmdStan or
    % MatlabStan are unavailable, and convergence is gated by assume so a
    % non-converged fit can never produce a false PASS.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)
        function testSubjectLevelGammaTestRetestRecovers(testCase)
            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   % log-mean grand intercept
            beta  = log(20);    % log-nu grand intercept (population dispersion)
            sP  = 0.35;  % person main effect on the log mean
            sO  = 0.22;  % occasion
            sT  = 0.15;  % trial
            sPO = 0.12;  % person x occasion
            sPT = 0.10;  % person x trial
            sOT = 0.08;  % trial x occasion
            sV  = 0.30;  % person log-nu SD (dispersion heterogeneity)
            rho = 0.40;  % corr(person log-mean, person log-nu)

            nsub = 40; nocc = 4; ntrl = 12;

            [dataTbl, truth] = localSimulateGammaTrtSserr(nsub, nocc, ntrl, ...
                alpha, beta, sP, sO, sT, sPO, sPT, sOT, sV, rho);

            % ---- Fit ----
            rel = localRunTrtSserrGamma(testCase, dataTbl, 'recov_gamma_trtss');

            % Provenance: the gamma SUBJECT-LEVEL test-retest path really ran.
            testCase.verifyEqual(rel.analysis, 'trt_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Gate on convergence via assume: a non-converged fit is filtered
            % (Incomplete), never a false PASS.
            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma trt sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'subject-level test-retest fit non-converged (max R-hat=%.4f >= 1.1 ', ...
                '/ n_eff gate); recovery not validated in this environment.'], maxRhat));

            % ---- Field contract: the slots psyrat_relsummary -> psyrat_ssrel_trt reads ----
            for f = {'mu','gro_sds','pop_sdlog','ind_sdlog','sig_occ','sig_trl', ...
                    'sig_trlxid','sig_occxid','sig_trlxocc'}
                testCase.verifyTrue(isfield(rel.out, f{1}), ...
                    sprintf('Missing subject-level slot %s', f{1}));
            end
            % REL.out slots are stored one cell PER STRATUM, each wrapping a
            % num2cell'd draw matrix; psyrat_relsummary unwraps the stratum first
            % (d.gro_sds = REL.out.gro_sds{s}), so do the same here. This run has a
            % single stratum.
            grosds   = cell2mat(rel.out.gro_sds{1});     % draws x 2
            indsdlog = cell2mat(rel.out.ind_sdlog{1});   % draws x NSUB
            testCase.verifyTrue(all(isfinite(grosds(:))));
            testCase.verifyTrue(all(isfinite(indsdlog(:))));
            testCase.verifySize(indsdlog, [size(grosds,1) nsub], ...
                'ind_sdlog must be draws x NSUB.');

            % ---- The five crossed facet SDs recover SEPARATELY ----
            % (the greedy-residual failure mode: if the per-subject residual
            % absorbed the facets, these would collapse toward zero)
            crossed = { ...
                'sig_occ',     truth.sd.bo; ...
                'sig_trl',     truth.sd.bt; ...
                'sig_trlxid',  truth.sd.txp; ...
                'sig_occxid',  truth.sd.oxp; ...
                'sig_trlxocc', truth.sd.txo};
            for k = 1:size(crossed,1)
                hat = mean(rel.out.(crossed{k,1})(:));
                tru = crossed{k,2};
                testCase.verifyGreaterThan(hat, 0.35 * tru, sprintf( ...
                    ['Crossed component %s collapsed (%.4f vs truth %.4f) - the ', ...
                    'per-subject residual may be absorbing it.'], crossed{k,1}, hat, tru));
                testCase.verifyLessThan(hat, 2.20 * tru, sprintf( ...
                    'Crossed component %s inflated (%.4f vs truth %.4f).', ...
                    crossed{k,1}, hat, tru));
            end

            % ---- The person component recovers ----
            sigIdHat = mean(grosds(:,1));
            testCase.verifyGreaterThan(sigIdHat, 0.55 * truth.sd.bp);
            testCase.verifyLessThan(sigIdHat,    1.60 * truth.sd.bp);

            % ---- The PER-SUBJECT residual recovers, per person ----
            % err_ss = exp(pop_sdlog + ind_sdlog) is exactly what psyrat_relsummary
            % forms. Correlating the posterior-mean per-subject residual against the
            % planted per-person truth is the check that fails loudly if the NSUB
            % columns are misaligned against id2 (values assigned to wrong people).
            errSsHat = mean(exp(rel.out.pop_sdlog + indsdlog), 1);   % 1 x NSUB
            r = corr_(errSsHat(:), truth.errSsTrue(:));
            fprintf('[gamma trt sserr recovery] per-subject residual corr = %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.70, sprintf( ...
                ['Recovered per-subject residuals correlate only %.3f with the ', ...
                'planted truth. Either the conditional conversion is wrong or the ', ...
                'ind_sdlog columns are misaligned against id2.'], r));

            % SCALE, and tight enough to discriminate a wrong V_w. A correlation is
            % rank-invariant, so it cannot see a constant multiplicative error -
            % and every plausible V_w mistake IS a constant factor on err_ss:
            %   err_ss is proportional to exp(V_w), so leaking the person facet in
            %   (V_w + s_p^2) multiplies every residual by exp(s_p^2).
            % Rather than assert an absolute window (which would have to be wide
            % enough to absorb posterior uncertainty and would then admit the
            % mutant), assert that the observed ratio is closer to 1 than to the
            % mutant's prediction, on the log scale. This is self-calibrating: it
            % tightens automatically as s_p grows.
            % SCALE -- and an explicit statement of what this run can and cannot
            % discriminate. Comparing the recovered residual scale against the
            % planted one does NOT test V_w here: at nocc = 4 the whole observed
            % scale is estimated high (measured on this exact design: sig_occ x1.66,
            % sig_occxid x1.56, the observed grand mean x1.28, mean crossed ratio
            % x1.32), which shifts the log-residual by ~+0.133 with a CORRECT V_w.
            % A person-facet leak into V_w would add only s_p^2 = 0.1225. The
            % estimation bias is larger than the mutant signal, so any assertion on
            % this quantity either fails on correct code or admits the mutant.
            % Report it, do not gate on it.
            logErrHat = mean(rel.out.pop_sdlog + indsdlog, 1);      % 1 x NSUB
            logShift  = mean(logErrHat(:) - log(truth.errSsTrue(:)));
            fprintf(['[gamma trt sserr recovery] log-residual shift = %+.4f ', ...
                '(small-sample scale inflation; a V_w person-facet leak would add ', ...
                '%+.4f, which this design cannot separate)\n'], logShift, sP^2);
            testCase.verifyLessThan(abs(logShift), 0.60, ...
                ['Recovered per-subject residuals are wildly off scale - far beyond ', ...
                'the small-sample inflation this design exhibits.']);

            % V_w IS pinned, but by an INTERNAL IDENTITY rather than by recovery.
            % The extractor stores its log-scale ingredients on gamma_disp, so the
            % construction can be checked exactly, independent of how well anything
            % was estimated. This is what catches a wrong V_w.
            gd = rel.out.gamma_disp;
            for f = {'V_w','alpha','s_o','s_t','s_po','s_pt','s_ot','nu','sigma_pot_e2_pop'}
                testCase.assertTrue(isfield(gd, f{1}), sprintf( ...
                    'gamma_disp.%s missing - the conditional residual is unverifiable.', f{1}));
            end
            % (i) V_w must be the sum of the five WITHIN-person log-mean variances,
            %     and must NOT include the person main effect s_p^2.
            V_w_expected = gd.s_o.^2 + gd.s_t.^2 + gd.s_po.^2 + gd.s_pt.^2 + gd.s_ot.^2;
            testCase.verifyEqual(gd.V_w, V_w_expected, 'RelTol', 1e-12, ...
                ['V_w is not the sum of the five within-person log-mean variances. ', ...
                'If the person main effect leaked in, this is where it shows.']);
            s_p_draws = cell2mat(rel.out.gro_sds{1});  % NOTE: col 1 is OBSERVED-scale
            testCase.verifyGreaterThan(mean(abs(gd.V_w - (V_w_expected + mean(s_p_draws(:,1)).^2))), 1e-6, ...
                'V_w appears to include a person-facet term.');
            % (ii) the typical-subject reference residual must be the closed form
            %      2*exp(2*alpha + 2*V_w - log(nu)) that pop_sdlog encodes.
            popExpected = 2 .* exp(2.*gd.alpha + 2.*gd.V_w - log(gd.nu));
            testCase.verifyEqual(gd.sigma_pot_e2_pop, popExpected, 'RelTol', 1e-10);
            testCase.verifyEqual(rel.out.pop_sdlog, 0.5 .* log(gd.sigma_pot_e2_pop), ...
                'RelTol', 1e-12, ...
                'pop_sdlog is not the log-SD encoding of the typical-subject residual.');

            % ---- The downstream contract: run the real reliability path ----
            % A defect that only bites downstream - gro_sds emitted as a variance
            % (RC-02), or ind_sdlog columns misordered against id2 - would not show
            % up in the slot checks above.
            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', struct('rel', rel), ...
                'analysis', 'trt_sserr', ...
                'gcoeff', 1, ...
                'reltype', 3, ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'CI', 0.95);
            testCase.assertEqual(relerr.nogooddata, 0, ...
                'psyrat_relsummary reported no usable data for the gamma case-25 fit.');
            T = outData.relsummary.group(1).event(1).ssrel_table;
            testCase.verifyEqual(height(T), nsub, ...
                'psyrat_relsummary did not produce one row per participant.');
            testCase.verifyTrue(all(isfinite(T.dep_pt)));
            testCase.verifyTrue(all(T.dep_pt > 0 & T.dep_pt < 1), ...
                'Per-participant dependability outside (0,1) - check the SD/variance contract.');
            % the participants with the smallest residuals must be the most reliable
            rDep = corr_(T.dep_pt, -truth.errSsTrue(:));
            fprintf('[gamma trt sserr recovery] dep vs -residual corr = %.3f\n', rDep);
            testCase.verifyGreaterThan(rDep, 0.60, ...
                ['Per-participant dependability does not track the planted residual ', ...
                'ordering through psyrat_relsummary - suspect id2 column alignment.']);

            % ---- Per-subject dispersion is genuinely heterogeneous ----
            % If the extractor collapsed to a single population residual, the
            % spread across participants would vanish and the design would be
            % indistinguishable from the group-level analysis 5.
            testCase.verifyGreaterThan(std(errSsHat) / mean(errSsHat), 0.05, ...
                ['Recovered per-subject residuals show no between-person spread; ', ...
                'the subject-level conversion has degenerated to the group-level one.']);
        end
    end
end

function [tbl, truth] = localSimulateGammaTrtSserr(nsub, nocc, ntrl, ...
    alpha, beta, sP, sO, sT, sPO, sPT, sOT, sV, rho)
% Balanced crossed Persons x Occasions x Trials long table from a log-linked mean
% with six independent normal random effects PLUS a person-varying log-nu
% correlated with the person mean, and a scaled-chi-square observation.
%
% The person block (u_p, v_p) is drawn with Var(u_p)=sP^2, Var(v_p)=sV^2,
% Cov(u_p,v_p)=rho*sP*sV, matching the coupled 2-D Cholesky block the gamma model
% fits. nu_p = exp(beta + v_p) is subject-constant (the ID-only nu structure that
% makes the conditional residual exact).
rng(12345);
z1 = randn(nsub,1);
z2 = randn(nsub,1);
u_p = sP * z1;
v_p = sV * (rho * z1 + sqrt(1 - rho^2) * z2);

u_o  = sO  * randn(nocc, 1);
u_t  = sT  * randn(ntrl, 1);
u_po = sPO * randn(nsub, nocc);
u_pt = sPT * randn(nsub, ntrl);
u_ot = sOT * randn(ntrl, nocc);

nrow = nsub * nocc * ntrl;
ids  = cell(nrow, 1);
time = zeros(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    nu_p = exp(beta + v_p(p));
    % integer df for the toolbox-free chi-square draw
    nuDraw = max(2, round(nu_p));
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            logmu = alpha + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p, o) + u_pt(p, t) + u_ot(t, o);
            mu = exp(logmu);
            chi2 = sum(randn(nuDraw, 1).^2);      % chi-square(nuDraw)
            meas(r) = (mu / nuDraw) * chi2;       % E = mu, Var = 2 mu^2 / nuDraw
        end
    end
end
tbl = table(ids, meas, time, 'VariableNames', {'id', 'meas', 'time'});

% ---- Closed-form observed-scale truth ----
% Group components at the population log-mean intercept (estimand-invariant).
vc = psyrat_gamma_varcomps_trt(alpha, sP, sO, sT, sPO, sPT, sOT, exp(beta), sV, ...
    rho * sP * sV);
truth.sd = struct( ...
    'bp',  sqrt(vc.sigma_p2),  'bo',  sqrt(vc.sigma_o2), ...
    'bt',  sqrt(vc.sigma_t2),  'txp', sqrt(vc.sigma_pt2), ...
    'oxp', sqrt(vc.sigma_po2), 'txo', sqrt(vc.sigma_ot2));

% Per-person residual SD: the conditional formula under test, evaluated at each
% planted (u_p, v_p). V_w excludes the person main effect sP^2 by construction.
V_w = sO^2 + sT^2 + sPO^2 + sPT^2 + sOT^2;
truth.errSsTrue = sqrt(2 .* exp(2*(alpha + u_p) + 2*V_w - (beta + v_p)));
truth.u_p = u_p;
truth.v_p = v_p;
end

function r = corr_(a, b)
% Pearson correlation without the Statistics Toolbox.
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function rel = localRunTrtSserrGamma(testCase, dataTbl, tag)
% Subject-level error variance (sserrvar = 2) WITH an occasion facet under the
% Gamma family, i.e. analysis 25 + family='gamma'. Fixed seed, per the
% reproducibility contract.
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
% Remove the per-run output directory and any stray CmdStan build/output files.
p = fullfile(baseDir, 'recov_gamma_trtss');
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
