classdef TestGammaDynrelLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE Gamma /
    % scaled-chi-square DYNAMIC/conditional reliability designs (analysis 11 and
    % its subject-level sibling 26, family='gamma', gammascale=2).
    %
    % WHY THIS FILE IS THE POINT OF THE WHOLE WORKSTREAM. The offline lane is
    % regression-only: the accuracy oracle duplicates the production formulas, so
    % a green suite says "unchanged", never "this new estimand is correct". Every
    % other check on gammascale=2 - unit tests, Stan snapshot goldens, stanc
    % parses - is compatible with a systematically wrong estimand. This is the
    % only check that is not.
    %
    % Data are simulated from the LOCATION-SCALE generative model, in which the
    % residual SD is a person-and-dimension quantity and the mean never enters it:
    %
    %   log mu_pt   = alpha   + b*z_p       + u_p + u_t
    %   log sigma_p = logsig0 + b_sigma*z_p + w_p        (w_p corr rho with u_p)
    %   nu_pt       = 2*(mu_pt/sigma_p)^2
    %   Y_pt ~ Gamma(nu/2, 2*mu/nu)   =>  E[Y]=mu_pt,  Var[Y]=sigma_p^2
    %
    % Note nu varies by TRIAL here, because mu does; sigma does not. That is the
    % exact feature that distinguishes this parameterization from the log-nu one
    % (gammascale=1), whose nu is person-constant. Under log-nu the residual is
    % 2*mu^2/nu and therefore tracks the trial-level mean; here it does not.
    %
    % THREE TESTS, ANSWERING THREE DIFFERENT QUESTIONS:
    %   1. testGroupSurfaceRecovered - does the fit recover the planted
    %      parameters, and does the reported G(z)/D(z) surface cover the truth?
    %   2. testCrossParameterizationIdentityHolds - fits the SAME data BOTH ways
    %      and checks b - b_nu/2 (log-nu) against b_sigma (location-scale). This
    %      is the acceptance criterion from
    %      documentation/gamma_scale_submodel_decision.md, and it is a genuine
    %      cross-implementation check: two different Stan models, two different
    %      MATLAB extraction paths, one algebraic identity linking them.
    %   3. testSubjectLevelPerParticipantRecovered - analysis 26, the only design
    %      that exercises the ind_sd extraction and psyrat_ssrel_dynrel_gamma_ls.
    %      It is also the only test that executes the analysis-26 arm of the
    %      gammascale=2 init block.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent),
    % and convergence-gated via assumeEqual, so a non-converged fit reports
    % Incomplete rather than a false PASS.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            base = localCmdStanOutputDir(testCase.projectRoot());
            for tag = {'recov_gamma_ls', 'recov_gamma_ls_lognu', 'recov_gamma_ls_sserr'}
                localCleanupArtifacts(base, tag{1});
            end
        end
    end

    methods (Test)

        function testGroupSurfaceRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');
            c = localTruth();
            [dataTbl, ~] = localSimulate(c);

            rel = localRunLs(testCase, dataTbl, 'recov_gamma_ls');
            testCase.verifyEqual(rel.analysis, 'ic_dynrel');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2, ...
                'The run did not record gammascale = 2.');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma LS dynrel recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'location-scale dynrel fit non-converged (max R-hat=%.4f); ', ...
                'recovery not validated in this environment.'], maxRhat));

            % ---- the location-scale slots must be the ones that got filled ----
            % This is also the live proof that the init-block fix works: before
            % it, the run died here with MATLAB:nonExistentField.
            testCase.verifyTrue(isfield(rel.out,'pop_sdlog') && ~isempty(rel.out.pop_sdlog));
            testCase.verifyTrue(isfield(rel.out,'b_sigma') && ~isempty(rel.out.b_sigma));
            testCase.verifyFalse(isfield(rel.out,'pop_lognu'), ...
                ['pop_lognu must be ABSENT under gammascale=2 - downstream ', ...
                'provenance gates on its existence.']);

            [est, draws] = localExtract(rel);
            fprintf(['[gamma LS dynrel recovery] b=%.3f (true %.2f)  b_sigma=%.3f (true %.2f)  ', ...
                's_p=%.3f (true %.2f)  s_sd=%.3f (true %.2f)  s_i=%.3f (true %.2f)  rho=%.3f (true %.2f)\n'], ...
                est.b, c.B, est.b_sigma, c.BSigma, est.s_p, c.SP, est.s_sd, c.SSd, ...
                est.s_i, c.SI, est.rho, c.Rho);

            % (1) mean dimension slope - the best identified quantity here
            testCase.verifyLessThan(abs(est.b - c.B), 0.15, sprintf( ...
                'Recovered log-mean slope %.3f far from true %.2f.', est.b, c.B));
            % (2) THE PARAMETER THIS PARAMETERIZATION EXISTS TO MAKE INTERPRETABLE.
            % Sign and magnitude both matter: b_sigma is meant to be readable as a
            % pure residual-SD effect, so a recovered value of the wrong sign
            % would invalidate the whole ruling.
            testCase.verifyGreaterThan(est.b_sigma, 0.05, sprintf( ...
                'Recovered b_sigma %.3f not clearly positive (true %.2f).', est.b_sigma, c.BSigma));
            testCase.verifyLessThan(abs(est.b_sigma - c.BSigma), 0.20, sprintf( ...
                'Recovered b_sigma %.3f far from true %.2f.', est.b_sigma, c.BSigma));
            % (3) the two person SDs and the trial SD
            testCase.verifyLessThan(abs(est.s_p - c.SP), 0.12, sprintf( ...
                'Recovered s_p %.3f far from true %.2f.', est.s_p, c.SP));
            testCase.verifyLessThan(abs(est.s_sd - c.SSd), 0.20, sprintf( ...
                'Recovered s_sd %.3f far from true %.2f.', est.s_sd, c.SSd));
            testCase.verifyLessThan(abs(est.s_i - c.SI), 0.12, sprintf( ...
                'Recovered s_i %.3f far from true %.2f.', est.s_i, c.SI));

            % (4) The recovered G(z)/D(z) surface must cover the true surface at
            % each grid z. This folds b, b_sigma, s_p, s_sd and s_i together
            % through the production calculator, so it fails if any of them, or
            % the conversion itself, is wrong.
            %
            % READ THIS BEFORE TREATING IT AS THE SHARP TEST, BECAUSE IT IS NOT.
            % At n' = 40 trials this design sits near the reliability ceiling:
            % measured on the 2026-08-03 run, true G moves only 0.974 -> 0.979
            % across the whole z grid while the posterior CIs are ~0.02 wide. So
            % coverage here is a SANITY check - it would catch a conversion that
            % returned 0.5, not one that got the z-slope modestly wrong. The
            % genuinely discriminating evidence in this file is elsewhere: direct
            % recovery of b_sigma in (2) above, the cross-implementation identity
            % in testCrossParameterizationIdentityHolds, and the partial
            % correlation in testSubjectLevelPerParticipantRecovered. If this file
            % is ever tightened, lowering n' (or raising the residual) to move the
            % coefficients off the ceiling is the change that would buy the most.
            %
            % Observed 2026-08-03: all six true values fell inside their CIs, but
            % D sat near the UPPER limit at all three z - a coherent pattern, not
            % noise, consistent with s_i recovering high (0.233 vs 0.20), since D
            % carries sigma_i^2 in its error term.
            zgrid = [-1 0 1];
            ci = 0.95;
            [trueG, trueD] = localTrueSurface(c, zgrid, ci);
            surf = psyrat_rel_dynrel_gamma_ls('alpha0',draws.alpha0,'b',draws.b, ...
                'logsig0',draws.logsig0,'b_sigma',draws.b_sigma,'sig_p',draws.s_p, ...
                'sig_sd',draws.s_sd,'sig_i',draws.s_i,'ndim',1,'z1',zgrid, ...
                'obs',c.NTrl,'CI',ci);
            for k = 1:numel(zgrid)
                fprintf('   z=%+d  G true %.3f in [%.3f, %.3f]   D true %.3f in [%.3f, %.3f]\n', ...
                    zgrid(k), trueG(k), surf.G.ll(k), surf.G.ul(k), ...
                    trueD(k), surf.D.ll(k), surf.D.ul(k));
                testCase.verifyGreaterThanOrEqual(trueG(k), surf.G.ll(k), sprintf( ...
                    'True G(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', zgrid(k), trueG(k), surf.G.ll(k)));
                testCase.verifyLessThanOrEqual(trueG(k), surf.G.ul(k), sprintf( ...
                    'True G(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', zgrid(k), trueG(k), surf.G.ul(k)));
                testCase.verifyGreaterThanOrEqual(trueD(k), surf.D.ll(k), sprintf( ...
                    'True D(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', zgrid(k), trueD(k), surf.D.ll(k)));
                testCase.verifyLessThanOrEqual(trueD(k), surf.D.ul(k), sprintf( ...
                    'True D(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', zgrid(k), trueD(k), surf.D.ul(k)));
            end

            % (5) The surface must move in the direction the TRUE surface moves.
            % Both slopes are positive here and they push G in OPPOSITE
            % directions (b raises it, b_sigma lowers it), so the sign of the
            % net effect is a property of the planted truth, not a free choice -
            % which is why this asserts agreement with trueG rather than a fixed
            % direction. A hard-coded direction would be a test of the simulation
            % settings, not of the fit.
            testCase.verifyEqual(sign(surf.G.pt(end) - surf.G.pt(1)), ...
                sign(trueG(end) - trueG(1)), ...
                'Recovered G(z) trends the opposite way from the true surface.');
        end

        function testCrossParameterizationIdentityHolds(testCase)
            % THE ACCEPTANCE CRITERION from
            % documentation/gamma_scale_submodel_decision.md. Both submodels
            % describe the same likelihood through
            %   log nu = log(2) + 2*log mu - 2*log sigma,
            % so their dimension slopes are linked by  b_nu = 2*b - 2*b_sigma,
            % i.e.  b - b_nu/2 = b_sigma.
            %
            % Fitting the SAME simulated data both ways therefore gives a
            % cross-implementation check: two different Stan models, two
            % different extraction branches, two different downstream slots, and
            % one identity that must hold across them. A sign error or a factor
            % of two anywhere in the location-scale path breaks it.
            %
            % NOTE the log-nu model is MILDLY MISSPECIFIED for this data by
            % construction: it forces nu to be person-constant, while the
            % generative nu varies with the trial effect u_t. That mismatch is
            % z-independent, so it should not bias the z-SLOPE identity; the
            % tolerance below is set for the resulting extra noise, not for a
            % systematic shift. This misspecification is the whole reason the two
            % parameterizations are different estimands.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required.');
            c = localTruth();
            [dataTbl, ~] = localSimulate(c);

            relLs = localRunLs(testCase, dataTbl, 'recov_gamma_ls');
            relNu = localRunLogNu(testCase, dataTbl, 'recov_gamma_ls_lognu');

            for r = {relLs, relNu}
                rr = psyrat_checkconv(r{1});
                testCase.assumeEqual(rr.out.conv.converged, 1, ...
                    'One of the two fits did not converge; identity not testable.');
            end

            bLs      = cell2mat(relLs.out.b{1});
            bSigmaLs = cell2mat(relLs.out.b_sigma{1});
            bNu      = cell2mat(relNu.out.b_nu{1});
            bLogNu   = cell2mat(relNu.out.b{1});

            % The two models share an IDENTICAL log-mean submodel, so their
            % log-mean slopes should agree closely. Checking this first isolates
            % the identity test below: if b already disagreed across the two
            % fits, a failure downstream would not tell you which side was wrong.
            fprintf('   b (LS fit) = %.3f;  b (log-nu fit) = %.3f;  true %.2f\n', ...
                mean(bLs), mean(bLogNu), c.B);
            testCase.verifyLessThan(abs(mean(bLs) - mean(bLogNu)), 0.10, sprintf( ...
                ['The two fits disagree on the log-mean slope (%.3f vs %.3f), ', ...
                'though they share an identical mean submodel. Investigate that ', ...
                'before reading the identity check below.'], mean(bLs), mean(bLogNu)));

            % implied b_sigma from the log-nu fit, as a posterior
            impliedBSigma = bLogNu - bNu./2;

            lo = quantile(impliedBSigma, 0.025);
            hi = quantile(impliedBSigma, 0.975);
            lsMean = mean(bSigmaLs);
            fprintf(['[cross-parameterization identity] b_sigma (LS fit) = %.3f ', ...
                '[%.3f, %.3f];  b - b_nu/2 (log-nu fit) = %.3f [%.3f, %.3f];  true %.2f\n'], ...
                lsMean, quantile(bSigmaLs,0.025), quantile(bSigmaLs,0.975), ...
                mean(impliedBSigma), lo, hi, c.BSigma);

            testCase.verifyGreaterThanOrEqual(lsMean, lo, sprintf( ...
                ['b_sigma from the location-scale fit (%.3f) falls below the ', ...
                'log-nu fit''s implied 95%% interval [%.3f, %.3f]. The identity ', ...
                'b_nu = 2*b - 2*b_sigma does not hold across the two ', ...
                'implementations.'], lsMean, lo, hi));
            testCase.verifyLessThanOrEqual(lsMean, hi, sprintf( ...
                ['b_sigma from the location-scale fit (%.3f) exceeds the log-nu ', ...
                'fit''s implied 95%% interval [%.3f, %.3f].'], lsMean, lo, hi));

            % SIGN CHECK, which is what the previously-wrong criterion 1 got
            % backwards: at a positive mean slope with a smaller positive
            % b_sigma, b_nu = 2*b - 2*b_sigma must be POSITIVE and roughly
            % 2*(b - b_sigma). Recorded rather than gated tightly, because it is
            % criterion 1 restated and criterion 2 above is the sharper test.
            fprintf('   b_nu (log-nu fit) = %.3f;  2*(b - b_sigma) at truth = %.3f\n', ...
                mean(bNu), 2*(c.B - c.BSigma));
            testCase.verifyGreaterThan(mean(bNu), 0, ...
                'b_nu should be positive when the mean slope exceeds the residual-SD slope.');
        end

        function testSubjectLevelPerParticipantRecovered(testCase)
            % Analysis 26 under gammascale=2. This is the ONLY test that runs the
            % analysis-26 arm of the init block, the ind_sd extraction, and
            % psyrat_ssrel_dynrel_gamma_ls end to end against a real fit.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required.');
            c = localTruth();
            [dataTbl, truth] = localSimulate(c);

            rel = localRunLsSserr(testCase, dataTbl, 'recov_gamma_ls_sserr');
            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyTrue(isfield(rel.out,'ind_sdlog') && ~isempty(rel.out.ind_sdlog), ...
                'The subject-level location-scale run stored no per-person residual block.');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma LS subject-level recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                'Subject-level fit non-converged (max R-hat=%.4f).', maxRhat));

            summ = psyrat_dynrel_summary(rel, 'CI', 0.95);
            testCase.verifyTrue(summ.isgamma_ls, ...
                'The summary did not route this run as location-scale gamma.');
            T = summ.strata(1).ssrel_table;
            testCase.assertNotEmpty(T, 'No per-participant table was produced.');
            testCase.verifyEqual(height(T), c.NSub);

            % Each participant's planted log residual-SD offset w_p should order
            % the recovered per-participant coefficients: more residual SD, lower
            % reliability. Correlation against the PLANTED w_p, not against a
            % re-derivation, so this is a known-truth check.
            r = corr_(truth.w_p(:), T.gen_pt(:));
            fprintf('[gamma LS subject-level recovery] corr(planted w_p, recovered G_s) = %.3f\n', r);
            testCase.verifyLessThan(r, -0.70, sprintf( ...
                ['Recovered per-participant coefficients track the planted ', ...
                'person residual-SD effects only at r=%.3f (expected strongly ', ...
                'negative). The conditional conversion is suspect.'], r));

            % The per-participant coefficients must actually vary; a degenerate
            % column would mean the conditional read-out collapsed to the
            % marginal one.
            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3, ...
                'Per-participant coefficients are constant; the read-out is not conditional.');

            % AND THE PROPERTY THAT DEFINES THIS PARAMETERIZATION: the
            % per-participant coefficient must be independent of the
            % participant's own log-MEAN effect u_p, once their own residual is
            % accounted for. Under gammascale=1 the residual is 2*mu^2/nu and u_p
            % enters directly, so this partial correlation is strongly negative
            % there. Here it must be near zero.
            rp = localPartialCorr(truth.u_p(:), T.gen_pt(:), truth.w_p(:));
            fprintf('[gamma LS subject-level recovery] partial corr(u_p, G_s | w_p) = %.3f\n', rp);
            testCase.verifyLessThan(abs(rp), 0.35, sprintf( ...
                ['Per-participant coefficients still depend on the participant''s ', ...
                'own log-mean (partial r=%.3f). Under gammascale=2 the residual ', ...
                'is decoupled from the mean, so this must be near zero; a strong ', ...
                'negative value means the log-nu residual formula leaked in.'], rp));
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function c = localTruth()
%LOCALTRUTH Planted truth, shared by all three tests so their numbers are
%   comparable. LogSig0 is the derived prior centre (psyrat_default_priors:
%   sigma_center = log(5.5) + 0.5*log(2/18)), which reproduces nu = 18 at the
%   reference mean - the SAME implied dispersion as the log-nu recovery twin's
%   beta = log(18). That is what makes the identity test a comparison of
%   parameterizations rather than of two dispersion levels.
c = struct( ...
    'Alpha',   log(5.5), ...   % log-mean intercept at z = 0
    'B',       0.30, ...       % log-mean dimension slope
    'SP',      0.35, ...       % person log-mean SD
    'SI',      0.20, ...       % trial  log-mean SD
    'LogSig0', log(5.5) + 0.5*log(2/18), ... % log residual-SD intercept at z = 0
    'BSigma',  0.25, ...       % log residual-SD dimension slope
    'SSd',     0.35, ...       % person log-sigma SD
    'Rho',     -0.30, ...      % person log-mean <-> log-sigma correlation
    'NSub',    120, ...
    'NTrl',    40);
end

function r = corr_(a, b)
%CORR_ Pearson correlation without the Statistics Toolbox (see the sibling
%   recovery tests and tests/TestBaseMatlabOnly.m).
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function [tbl, truth] = localSimulate(c)
% Persons x Trials long table generated from the LOCATION-SCALE model. The
% person (log-mean, log-sigma) block is bivariate normal with correlation rho.
% Y = gamrnd(nu/2, 2*mu/nu) has E[Y] = mu and Var[Y] = 2*mu^2/nu = sigma_p^2,
% by the definition nu = 2*(mu/sigma)^2 - so the residual SD is EXACTLY
% sigma_p, independent of mu. That independence is the generative feature the
% whole parameterization is about, and it is what makes the log-nu model
% (person-constant nu) a genuinely different model for this data.
rng(20260803);
u_t = c.SI * randn(c.NTrl, 1);
z1 = randn(c.NSub, 1); z2 = randn(c.NSub, 1);
u_p = c.SP * z1;
w_p = c.SSd * (c.Rho * z1 + sqrt(1 - c.Rho^2) * z2);

% standardized subject-level covariate (exact mean 0, SD 1 over subjects), so
% the pipeline's own z-scoring reproduces it
raw = randn(c.NSub, 1);
zp = (raw - mean(raw)) / std(raw);

sigma_p = exp(c.LogSig0 + c.BSigma * zp + w_p);

nrow = c.NSub * c.NTrl;
ids = cell(nrow,1); meas = zeros(nrow,1); dim1 = zeros(nrow,1);
r = 0;
for p = 1:c.NSub
    for t = 1:c.NTrl
        r = r + 1;
        ids{r} = sprintf('S%03d', p);
        mu = exp(c.Alpha + c.B * zp(p) + u_p(p) + u_t(t));
        nu = 2 * (mu / sigma_p(p))^2;
        meas(r) = gamrnd(nu / 2, 2 * mu / nu);
        dim1(r) = zp(p);
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});
truth = struct('u_p', u_p, 'w_p', w_p, 'zp', zp, 'sigma_p', sigma_p);
end

function [trueG, trueD] = localTrueSurface(c, zgrid, ci)
% Closed-form true G(z)/D(z) via the production converter AT THE TRUE
% PARAMETERS. Using the production converter here is deliberate and is not
% circular for this purpose: what is under test is whether the FIT recovers
% parameters that reproduce the truth, not whether the converter matches an
% independent transcription (that is the offline oracle's job).
trueG = zeros(1,numel(zgrid));
trueD = zeros(1,numel(zgrid));
for k = 1:numel(zgrid)
    z = zgrid(k);
    vc = psyrat_gamma_varcomps_ls(c.Alpha + c.B*z, c.SP, c.SI, ...
        c.LogSig0 + c.BSigma*z, c.SSd);
    [~,trueG(k)] = psyrat_rel_sing('gcoeff',2,'metric','global', ...
        'bp',sqrt(vc.sigma_p2),'wp',sqrt(vc.sigma_pi_e2+vc.sigma_i2), ...
        'i',sqrt(vc.sigma_i2),'obs',c.NTrl,'CI',ci);
    [~,trueD(k)] = psyrat_rel_sing('gcoeff',1,'metric','global', ...
        'bp',sqrt(vc.sigma_p2),'wp',sqrt(vc.sigma_pi_e2+vc.sigma_i2), ...
        'i',sqrt(vc.sigma_i2),'obs',c.NTrl,'CI',ci);
end
end

function [est, draws] = localExtract(rel)
% Posterior draws from the single stratum, in the slots gammascale=2 fills.
draws = struct();
draws.alpha0  = rel.out.mu(:,1);
draws.logsig0 = rel.out.pop_sdlog(:,1);
gro           = cell2mat(rel.out.gro_sds{1});   % draws x 2 [s_p, s_sd]
draws.s_p     = gro(:,1);
draws.s_sd    = gro(:,2);
draws.b       = cell2mat(rel.out.b{1});
draws.b_sigma = cell2mat(rel.out.b_sigma{1});
draws.s_i     = rel.out.sig_trl(:,1);
chol          = cell2mat(rel.out.chol_corrmat{1});
draws.rho     = chol(:,2,1);

est = struct('b', mean(draws.b), 'b_sigma', mean(draws.b_sigma), ...
    's_p', mean(draws.s_p), 's_sd', mean(draws.s_sd), ...
    's_i', mean(draws.s_i), 'rho', mean(draws.rho), ...
    'alpha0', mean(draws.alpha0), 'logsig0', mean(draws.logsig0));
end

function r = localPartialCorr(x, y, z)
% Partial correlation of x and y controlling for z, by residualizing both on z.
% Written out rather than calling partialcorr (Statistics Toolbox) so this
% helper stays base-MATLAB; see tests/TestBaseMatlabOnly.m.
X = [ones(numel(z),1), z(:)];
rx = x(:) - X * (X \ x(:));
ry = y(:) - X * (X \ y(:));
r = (rx' * ry) / sqrt((rx' * rx) * (ry' * ry));
end

function rel = localRunLs(testCase, dataTbl, tag)
% dynrel=2 (a Dimension 1 column) + family gamma + gammascale=2 -> analysis 11
% under the location-scale scale submodel.
rel = localRun(testCase, dataTbl, tag, {'gammascale', 2});
end

function rel = localRunLogNu(testCase, dataTbl, tag)
% The same data under the log-nu scale submodel, for the identity check.
% gammascale=1 must be EXPLICIT: log-nu stopped being the default when ruling
% 33 (2026-08-07) flipped the movable analyses to location-scale, and an
% omitted gammascale now resolves to 2 - this arm would silently refit the
% location-scale model.
rel = localRun(testCase, dataTbl, tag, {'gammascale', 1});
end

function rel = localRunLsSserr(testCase, dataTbl, tag)
% dynrel=2 + sserrvar=2 + gammascale=2 -> analysis 26, location-scale.
rel = localRun(testCase, dataTbl, tag, {'gammascale', 2, 'sserrvar', 2});
end

function rel = localRun(testCase, dataTbl, tag, extraArgs)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'dynrel', 2, ...
    extraArgs{:}, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir, tag)
p = fullfile(baseDir, tag);
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            'Could not remove recovery artifact directory ''%s''.', p);
    end
end
pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'output-*.csv', 'temp.data.R'};
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
