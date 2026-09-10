classdef TestGammaSserrLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL design under the LOCATION-SCALE parameterization
    % (analysis 6, ic_sserrvar, gammascale = 2).
    %
    % Data are simulated from the LOCATION-SCALE generative model,
    %   log mu_pt   = alpha + u_p + u_t
    %   log sigma_p = logsig0 + w_p
    %   nu_pt       = 2*(mu_pt/sigma_p)^2,   Y ~ scaled chi-square(mu_pt, nu_pt)
    % so Var(Y | p) = sigma_p^2 EXACTLY and nu varies by trial while sigma does
    % not. That asymmetry is what makes this a different model from the log-nu
    % form rather than a reparameterization of it, so simulating from it is the
    % only way the recovery can discriminate the two.
    %
    % WHAT THIS TEST DOES AND DOES NOT SHOW -- read before adding assertions.
    %
    % It was originally written to assert that a log-nu fit to the same data would
    % score high-amplitude participants as less reliable while the location-scale
    % fit would not. MEASURED 2026-08-04: both fits give essentially the same
    % per-participant coefficients (partial corr 0.1538 vs 0.1499). The log-nu
    % model has ind_nu free per person, so it absorbs the amplitude term
    % (v_p = const + 2*u_p - 2*w_p) and recovers the correct residual. At a
    % well-identified design point the parameterization does NOT change the
    % numbers.
    %
    % That is consistent with owner ruling H, which never claimed otherwise: the
    % log-nu residual "is exact for the current proportional-error model; it is
    % not a coding bug", and the difference lies in the SHRINKAGE TARGET (constant
    % CV vs constant absolute SD, which bites where person dispersion is weakly
    % identified) and in what the parameter means. Finding S16's formula-level
    % claim -- that at EQUAL ind_nu a higher-amplitude participant gets a lower
    % coefficient -- is pinned offline in TestGammaSserrLsConversion and is not in
    % question here.
    %
    % So this test now checks (1) that the location-scale scale submodel is
    % recovered, (2) that the per-participant read-out tracks the planted person
    % scale effect, and (3) the structural CONTRAST that does survive: both forms
    % recover the same per-person residual, but log-nu needs a far larger
    % person-dispersion SD to do it, because its parameter carries amplitude as
    % well as precision.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).

    % Copyright (C) 2016-2025 Peter E. Clayson
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program (gpl.txt). If not, see
    %     <http://www.gnu.org/licenses/>.

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
        function testLocationScaleSubjectLevelRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth ----
            alpha   = log(5.5);                    % log-mean grand intercept
            s_p     = 0.40;                        % person log-mean SD
            s_i     = 0.20;                        % trial  log-mean SD
            logsig0 = log(5.5) + 0.5*log(2/18);    % log residual-SD intercept
            s_sd    = 0.50;                        % person log-sigma SD
            rho     = -0.30;                       % person mean <-> log-sigma correlation

            nsub = 60; ntrl = 40;
            [dataTbl, u_p, w_p] = localSimulateGammaSserrLs( ...
                nsub, ntrl, alpha, s_p, s_i, logsig0, s_sd, rho);

            % ---- Closed-form observed-scale truth (POOLED, S17-FIX) ----
            % The numerator is a POPULATION quantity and is parameterization-
            % invariant; each participant's residual pools the mean-surface
            % person x trial term exp(2*alpha)*K (population-built - K written
            % out from the section-J identity, an independent oracle) with the
            % participant's own conditional variance. No dependence on their own
            % mean, and no trial-averaging factor on the conditional term.
            vc = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, logsig0, s_sd);
            sigma_p2_true = vc.sigma_p2;
            K = exp(s_p^2 + s_i^2) * (exp(s_p^2) - 1) * (exp(s_i^2) - 1);
            pi_mean_true = exp(2*alpha) * K;
            trueResid = pi_mean_true + exp(2*(logsig0 + w_p));             % nsub x 1
            trueG = sigma_p2_true ./ (sigma_p2_true + trueResid ./ ntrl);

            % ---- Fit ----
            rel = localRunSserrGammaLs(testCase, dataTbl, 'recov_gamma_ls_sserr');

            testCase.verifyEqual(rel.analysis, 'ic_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Recovered per-subject reliability (mirrors psyrat_relsummary) ----
            idtable = rel.out.id_matches{1};
            idtable.trls = repmat(ntrl, height(idtable), 1);
            ssrel = psyrat_ssrel( ...
                'bp',     rel.out.gro_sds{1}(:,1), ...
                'wp_pop', rel.out.pop_sdlog(:,1), ...
                'wp_ss',  rel.out.ind_sdlog{1}, ...
                'i',      rel.out.sig_trl(:,1), ...
                'gcoeff', 2, ...
                'idtable', idtable, ...
                'CI', 0.95);

            nrows = height(ssrel);
            idx = arrayfun(@(r) sscanf(char(string(ssrel.id(r))),'S%d'), (1:nrows)');
            recG   = ssrel.dep_pt;
            recLL  = ssrel.dep_ll;
            recUL  = ssrel.dep_ul;
            recRes = ssrel.ss_errvar;
            tG     = trueG(idx);
            tRes   = trueResid(idx);
            plantedU = u_p(idx);     % person log-MEAN effect
            plantedW = w_p(idx);     % person log-SIGMA effect

            % (1) the scale submodel is recovered. Since S17-FIX pop_sdlog is
            % the POOLED typical-subject residual log-SD, so its reference is
            % the pooled truth, not the bare Intercept_sigma.
            popTrue = 0.5*log(pi_mean_true + exp(2*logsig0));
            sigHat = mean(rel.out.pop_sdlog(:,1));
            sdHat  = mean(cell2mat(rel.out.gro_sds{1}(:,2)));
            rResid = localPearson(recRes, tRes);
            coverage = mean(tG >= recLL & tG <= recUL);

            % (3) THE AMPLITUDE CHECK, referenced to the TRUTH rather than to zero.
            %
            % CORRECTED 2026-08-04, third revision, and the reason is worth
            % keeping: this quantity is NOT zero on the true coefficients. D (and
            % G) is a NONLINEAR function of w_p, while partial correlation removes
            % only the LINEAR component of w_p; with u_p and w_p correlated in the
            % generating model, the leftover nonlinearity shows up in the u_p term.
            % Measured on the planted values with no estimation at all, the true
            % partial correlation is about +0.12 at n' = 48 and +0.30 at n' = 12.
            % Asserting "must be ~0" was therefore wrong regardless of what the
            % estimator does. The right reference is the truth's own value.
            rGU_partial     = localPartialCorr(plantedU, recG, plantedW);
            rGU_partial_true = localPartialCorr(plantedU, tG, plantedW);
            rGW             = localPearson(plantedW, recG);
            rGW_true        = localPearson(plantedW, tG);

            fprintf(['[LS sserr recovery] pooled pop_sdlog hat=%.3f (true %.3f)  ', ...
                's_sd hat=%.3f (true %.2f)\n  corr(resid)=%.3f  coverage=%.2f  ', ...
                'mean recG=%.3f (true %.3f)\n  partial corr(u_p, G | w_p): recovered %.4f  ', ...
                'TRUE %.4f  (the truth is NOT zero - G is nonlinear in w_p)\n  ', ...
                'corr(w_p, G): recovered %.3f  TRUE %.3f\n'], ...
                sigHat, popTrue, sdHat, s_sd, rResid, coverage, ...
                mean(recG), mean(tG), rGU_partial, rGU_partial_true, rGW, rGW_true);

            testCase.verifyLessThan(abs(sigHat - popTrue), 0.25, sprintf( ...
                'Recovered pooled population residual log-SD %.3f vs true %.3f.', sigHat, popTrue));
            testCase.verifyGreaterThan(sdHat, 0.30, sprintf( ...
                'Recovered person log-sigma SD %.3f far below the true %.2f.', sdHat, s_sd));
            testCase.verifyLessThan(sdHat, 0.75);

            % (2) subject-level differentiation
            testCase.verifyGreaterThan(rResid, 0.5, sprintf( ...
                ['Recovered per-subject residual does not track the truth ', ...
                '(Pearson r=%.3f); subject-level differentiation not recovered.'], rResid));
            testCase.verifyGreaterThanOrEqual(coverage, 0.80, sprintf( ...
                'Only %.0f%% of true per-subject coefficients inside the 95%% CIs.', ...
                100*coverage));
            testCase.verifyLessThan(abs(mean(recG) - mean(tG)), 0.05, sprintf( ...
                'Sample-mean per-subject reliability off: recovered %.3f vs true %.3f.', ...
                mean(recG), mean(tG)));

            testCase.verifyLessThan(rGW, -0.5, sprintf( ...
                ['corr(planted person scale effect, recovered coefficient) = %.3f. ', ...
                'The per-participant read-out must track the person quantity it ', ...
                'is supposed to, strongly and with the right sign.'], rGW));

            % The recovered amplitude relationship must MATCH THE TRUTH'S, not be
            % zero. Any departure is estimation error, and this is the assertion
            % that would catch a genuine re-coupling of the residual to the mean:
            % that would drive the recovered value well above the true one.
            testCase.verifyLessThan(abs(rGU_partial - rGU_partial_true), 0.20, sprintf( ...
                ['Recovered partial corr(u_p, G | w_p) = %.4f against a TRUE value ', ...
                'of %.4f. The recovered coefficient should carry the same amplitude ', ...
                'relationship the true coefficient does; a substantially larger ', ...
                'value means mean dependence has entered that is not in the truth.'], ...
                rGU_partial, rGU_partial_true));

            % (3) THE STRUCTURAL CONTRAST, and READ THE NOTE BEFORE CHANGING IT.
            %
            % MEASURED 2026-08-04, and it corrected this test rather than the code.
            % This block first asserted that a log-nu fit to the SAME data would
            % give a strongly NEGATIVE partial corr(u_p, recG | w_p) while the
            % location-scale fit gave ~0. Both fits gave essentially the SAME
            % value (0.1538 vs 0.1499). That is not a defect and not noise: the
            % log-nu model has ind_nu FREE per person, so fitted to
            % location-scale data the likelihood simply pushes
            %   v_p = const + 2*u_p - 2*w_p
            % and recovers the correct per-participant residual. At a
            % well-identified design point BOTH parameterizations produce the same
            % per-participant coefficients.
            %
            % That is consistent with owner ruling H, which says the log-nu
            % residual "is exact for the current proportional-error model; it is
            % not a coding bug", and locates the difference in the SHRINKAGE
            % TARGET (constant CV vs constant absolute SD) and in what the
            % parameter means. Do not reinstate a claim that log-nu gets
            % per-participant reliability wrong at this n - it does not.
            %
            % What IS structural, and what this now asserts, is the MECHANISM: the
            % log-nu model has to SPEND its person-dispersion parameter absorbing
            % amplitude, and the location-scale model does not. From
            % v_p = const + 2*u_p - 2*w_p with the planted s_p/s_sd/rho, the
            % implied corr(u_p, v_p) is about +0.76 and sd(v_p) about 1.46 against
            % the location-scale s_sd of 0.50. Those are large, structural, and
            % survive any amount of data, which is what the previous formulation
            % wrongly claimed of the coefficient itself.
            relNu = localRunSserrGammaLogNu(testCase, dataTbl, 'recov_gamma_lognu_contrast');
            idtableNu = relNu.out.id_matches{1};
            idtableNu.trls = repmat(ntrl, height(idtableNu), 1);
            ssrelNu = psyrat_ssrel( ...
                'bp',     relNu.out.gro_sds{1}(:,1), ...
                'wp_pop', relNu.out.pop_sdlog(:,1), ...
                'wp_ss',  relNu.out.ind_sdlog{1}, ...
                'i',      relNu.out.sig_trl(:,1), ...
                'gcoeff', 2, 'idtable', idtableNu, 'CI', 0.95);
            idxNu = arrayfun(@(r) sscanf(char(string(ssrelNu.id(r))),'S%d'), ...
                (1:height(ssrelNu))');
            rGU_partial_nu = localPartialCorr(u_p(idxNu), ssrelNu.dep_pt, w_p(idxNu));

            % NOT measurable here, recorded so it is not re-attempted: the natural
            % way to show the absorption would be corr(planted u_p, recovered
            % ind_nu). Case 6 does not store per-person ind_nu - only the SD-
            % equivalent encoding ind_sdlog - and that encoding is
            %   ind_sdlog = u_p - 0.5*v_p = w_p   (substituting v_p above),
            % i.e. it is the SAME residual log-SD under both parameterizations by
            % construction. So ind_sdlog cannot show the absorption; it is
            % precisely the quantity that agrees. What CAN be seen from the stored
            % slots is the two consequences below.
            indSdLs = mean(cell2mat(rel.out.ind_sdlog{1}), 1)';
            indSdNu = mean(cell2mat(relNu.out.ind_sdlog{1}), 1)';
            rEncoding = localPearson(indSdLs(idx), indSdNu(idxNu));
            svHatNu = mean(cell2mat(relNu.out.gro_sds{1}(:,2)));

            fprintf(['  CONTRAST on the same data:\n', ...
                '    partial corr(u_p, recG | w_p): LS = %.4f   log-nu = %.4f  ', ...
                '(EXPECTED to agree; see the note)\n', ...
                '    corr(per-person residual log-SD, LS vs log-nu) = %.3f\n', ...
                '    person dispersion SD: LS s_sd = %.3f   log-nu s_v = %.3f\n'], ...
                rGU_partial, rGU_partial_nu, rEncoding, sdHat, svHatNu);

            % (3a) THE MEASURED FINDING, now pinned rather than assumed: at this
            % design point the two parameterizations recover the SAME per-person
            % residual. If this ever stops holding, one of the two extractors has
            % changed estimand and the decision record's "not a coding bug"
            % framing needs revisiting.
            testCase.verifyGreaterThan(rEncoding, 0.90, sprintf( ...
                ['The two parameterizations recovered materially different ', ...
                'per-person residuals (corr = %.3f). At a well-identified design ', ...
                'point they are expected to agree - the log-nu model has ind_nu ', ...
                'free and absorbs the amplitude term.'], rEncoding));

            % (3b) THE STRUCTURAL COST, which is what actually differs. The log-nu
            % model carries v_p = const + 2*u_p - 2*w_p, so its person-dispersion
            % SD must be far larger than the location-scale s_sd to describe the
            % same data: about 1.46 against 0.50 for the planted values. That is
            % the "spending its dispersion parameter on amplitude" claim made
            % quantitative, and unlike the coefficient it does not wash out with n.
            testCase.verifyGreaterThan(svHatNu, 1.5 * sdHat, sprintf( ...
                ['The log-nu person-dispersion SD (%.3f) should substantially ', ...
                'exceed the location-scale one (%.3f), because it carries the ', ...
                'amplitude term as well as the precision term. If they are ', ...
                'comparable, the log-nu fit is not absorbing what the algebra ', ...
                'says it must.'], svHatNu, sdHat));
        end
    end
end

function [tbl, u_p, w_p] = localSimulateGammaSserrLs(nsub, ntrl, alpha, s_p, s_i, logsig0, s_sd, rho)
% One-facet Persons x Trials long table from the LOCATION-SCALE generative model.
% Y = gamrnd(nu/2, 2*mu/nu) has E[Y] = mu and Var[Y] = 2*mu^2/nu; with
% nu = 2*(mu/sigma)^2 that is Var[Y] = sigma^2, independent of mu. Note nu is
% computed per OBSERVATION (it moves with the trial effect) while sigma is per
% PARTICIPANT -- the asymmetry that distinguishes this model from the log-nu form.
rng(12345);
u_t = s_i * randn(ntrl, 1);                          % crossed trial main effects
z1 = randn(nsub, 1); z2 = randn(nsub, 1);
u_p = s_p  * z1;                                     % person mean effect
w_p = s_sd * (rho * z1 + sqrt(1 - rho^2) * z2);      % person log-sigma effect (corr rho)
sigma_p = exp(logsig0 + w_p);

nrow = nsub * ntrl;
ids  = cell(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    for t = 1:ntrl
        r = r + 1;
        ids{r} = sprintf('S%03d', p);
        mu = exp(alpha + u_p(p) + u_t(t));
        nu = 2 * (mu / sigma_p(p))^2;
        meas(r) = gamrnd(nu / 2, 2 * mu / nu);
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function r = localPearson(x, y)
% Pearson correlation via base MATLAB (no Statistics Toolbox).
x = x(:) - mean(x(:)); y = y(:) - mean(y(:));
r = (x' * y) / (sqrt(x' * x) * sqrt(y' * y));
end

function r = localPartialCorr(x, y, z)
% Partial correlation of x and y controlling for z, from the three pairwise
% Pearson correlations. Base MATLAB only.
rxy = localPearson(x, y);
rxz = localPearson(x, z);
ryz = localPearson(y, z);
r = (rxy - rxz*ryz) / sqrt((1 - rxz^2) * (1 - ryz^2));
end

function rel = localRunSserrGammaLs(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function rel = localRunSserrGammaLogNu(testCase, dataTbl, tag)
% The CONTRAST fit: identical data and settings, log-nu parameterization. Its
% per-participant residual contains exp(2*u_p) structurally, so the partial
% correlation this test compares against is not a matter of estimation noise.
% gammascale=1 must be EXPLICIT: this helper originally omitted it back when
% log-nu was the default, but ruling 33 (2026-08-07) flipped the movable
% analyses to location-scale, so an omitted gammascale now silently refits the
% location-scale model.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'sserrvar', 2, ...
    'gammascale', 1, ...
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

function localCleanupArtifacts(baseDir)
for p = {fullfile(baseDir, 'recov_gamma_ls_sserr'), ...
         fullfile(baseDir, 'recov_gamma_lognu_contrast')}
    if isfolder(p{1})
        try
            rmdir(p{1}, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                'Could not remove recovery artifact directory ''%s''.', p{1});
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
