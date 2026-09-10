classdef TestGammaDynrelTrtLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the LOCATION-SCALE Gamma /
    % scaled-chi-square TWO-facet (trial + occasion) DYNAMIC/conditional
    % reliability design (analysis 14, family='gamma', gammascale=2). Increment 8
    % of the gammascale rollout; the first design that leaves the static family.
    %
    % Sibling of TestGammaDynrelTrtRecovery, which covers the log-nu twin.
    % Crossed Persons x Occasions x Trials data are simulated from a log-linked
    % mean (six independent random effects) and a person-varying log RESIDUAL SD,
    % both conditioned on a standardized subject-level dimension z:
    %
    %   log mu_pot   = alpha   + b*z_p       + u_p + u_o + u_t + u_po + u_pt + u_ot
    %   log sigma_p  = logsig0 + b_sigma*z_p + w_p      (w_p corr rho with u_p)
    %   nu           = 2*(mu/sigma)^2,   Var(Y | mu, sigma) = sigma^2 exactly
    %
    % WHAT IS ASSERTED. The mean slope b, the residual-SD slope b_sigma, the
    % person-dispersion SD s_sd, and - decisively - the z-conditioned two-facet
    % coefficient surface (reltype 3) must cover the closed-form true surface at
    % z in {-1, 0, +1}.
    %
    % b_sigma IS ASSERTED, unlike increment 6/7's static residual-SD CONTRAST,
    % which those tests printed as a diagnostic instead. The reason is structural,
    % not a change of standard: there the quantity was a DERIVED difference of
    % per-cell residuals, where per-cell estimation errors compound; here b_sigma
    % is a directly sampled parameter of the fitted model. It is the interpretive
    % payload of the whole parameterization, so leaving it unasserted would leave
    % the increment's central claim untested.
    %
    % ---- THE DESIGN POINT, AND THE TRAP IT AVOIDS ----
    % Planted: b = +0.30, b_sigma = -0.20, so the implied log-nu slope is
    % b_nu = 2*b - 2*b_sigma = +1.00. The blend is therefore maximal: a log-nu fit
    % of this data would report a dispersion slope of 1.00, of which 0.60 is PURE
    % AMPLITUDE and only the remainder is precision. That is exactly the
    % conflation ruling H exists to break.
    %
    % THE DYNAMIC FORM OF INCREMENT 7's CANCELLATION TRAP, found offline before
    % any fit was run. Sweeping b_sigma at fixed b gives:
    %
    %   b_sigma   b_nu    G(-1)  G(0)  G(+1)   G swing   residual SD across z
    %    -0.20    1.00    0.827  0.891  0.917   +0.090     4.31 -> 2.94
    %     0.00    0.60    0.859  0.891  0.909   +0.050     3.53 -> 3.56
    %     0.30    0.00    0.891  0.891  0.891   -0.000     2.63 -> 4.78
    %     0.60   -0.60    0.909  0.891  0.859   -0.050     1.96 -> 6.44
    %
    % At b_sigma = b the coefficient surface is EXACTLY FLAT while the residual SD
    % nearly doubles across the grid. A recovery test planted there would report a
    % healthy G ~ 0.891 recovered at every z while being completely blind to the
    % residual-SD slope - the one quantity this parameterization exists to
    % isolate - and it is also precisely where b_nu = 0, so the log-nu twin sees no
    % dispersion slope at all. Doubly degenerate, and it looks fine. The planted
    % point is chosen at the OTHER end: the largest coefficient gradient in the
    % sweep and the largest implied blend.
    %
    % G ~ [0.83, 0.92] is deliberately away from the ceiling; increment 4's ledger
    % records G ~ 0.98 as a weak test of the coefficient surface.
    %
    % ---- THE DATA-GENERATING PROCESS IS EXACT, NOT APPROXIMATE ----
    % y = (mu/nu)*chi2(nu) is exactly Gamma(nu/2, 2*mu/nu), drawn by
    % Marsaglia-Tsang in base MATLAB. Two reasons this matters here. First, the
    % log-nu sibling can sum nu squared normals because its nu is a fixed INTEGER;
    % under location-scale nu = 2*(mu/sigma)^2 is derived, non-integer, and
    % person- AND trial-varying, so the small-nu tail is genuinely visited and a
    % Wilson-Hilferty approximation would put an unquantified bias in the very
    % ground truth this test compares against. Second, the log-nu sibling assumes
    % Statistics Toolbox gamrnd and skips without it; this test uses base MATLAB
    % only, so it survives CI's --products=MATLAB image.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            base = localCmdStanOutputDir(testCase.projectRoot());
            for tag = {'recov_gamma_dynrel_trt_ls', ...
                    'recov_gamma_dynrel_trt_ls_pairA', ...
                    'recov_gamma_dynrel_trt_ls_pairB'}
                localCleanupArtifacts(base, tag{1});
            end
        end
    end

    methods (Test)
        function testTwoFacetLocationScaleDynrelRecovered(testCase)
            T = localTruth();
            zgrid = [-1 0 1]; reltype = 3; ci = 0.95;

            dataTbl = localSimulateGammaDynrelTrtLs(T);

            % ---- Closed-form true G(z)/D(z) at each grid z, at the design
            % trial/occasion counts. Argument order matches the production
            % calculator: (alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, log_sigma, s_sd)
            % with s_po = person x occasion and s_pt = person x trial. ----
            trueG = zeros(1,numel(zgrid));
            trueD = zeros(1,numel(zgrid));
            trueRes = zeros(1,numel(zgrid));
            for k = 1:numel(zgrid)
                z = zgrid(k);
                vc = psyrat_gamma_varcomps_trt_ls(T.alpha + T.b*z, T.sP, T.sO, ...
                    T.sT, T.sPO, T.sPT, T.sOT, T.logsig0 + T.b_sigma*z, T.s_sd);
                fa = {'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2),...
                      'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2),...
                      'txo',sqrt(vc.sigma_ot2)};
                [~,trueG(k)] = psyrat_rel_trt('gcoeff',2,'reltype',reltype,...
                    'bp',sqrt(vc.sigma_p2),fa{:},'err',sqrt(vc.sigma_res2),...
                    'obs',T.ntrl,'nocc',T.nocc,'CI',ci);
                [~,trueD(k)] = psyrat_rel_trt('gcoeff',1,'reltype',reltype,...
                    'bp',sqrt(vc.sigma_p2),fa{:},'err',sqrt(vc.sigma_res2),...
                    'obs',T.ntrl,'nocc',T.nocc,'CI',ci);
                trueRes(k) = sqrt(vc.sigma_res2);
            end

            % ---- Fit ----
            rel = localRunDynrelTrtGammaLs(testCase, dataTbl, ...
                'recov_gamma_dynrel_trt_ls');
            testCase.verifyEqual(rel.analysis, 'ic_dynrel_trt');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');
            testCase.verifyEqual(rel.gammascale, 2);

            % Convergence first: a coefficient recovered from a non-converged fit
            % says nothing.
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf(['\n[dynrel_trt LS recovery] max R-hat %.4f over %d params ' ...
                '(converged = %d)\n'], maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Fit did not converge (max R-hat %.4f); recovery is not ' ...
                'interpretable, so the run is skipped rather than failed.'], maxRhat));

            % ---- Recovered draws (single stratum). These are the LOCATION-SCALE
            % slots; reading pop_lognu/b_nu here would error, which is the point of
            % the separate init branch. ----
            alpha0  = rel.out.mu(:,1);
            logsig0 = rel.out.pop_sdlog(:,1);
            gro     = cell2mat(rel.out.gro_sds{1});   % draws x 2 [s_p, s_sd]
            bdraw   = cell2mat(rel.out.b{1});
            bsdraw  = cell2mat(rel.out.b_sigma{1});
            chol    = cell2mat(rel.out.chol_corrmat{1});
            rho_ps  = chol(:,2,1);
            s_occ   = rel.out.sig_occ(:,1);
            s_trl   = rel.out.sig_trl(:,1);
            s_txp   = rel.out.sig_trlxid(:,1);
            s_oxp   = rel.out.sig_occxid(:,1);
            s_txo   = rel.out.sig_trlxocc(:,1);

            testCase.verifyFalse(isfield(rel.out, 'pop_lognu'), ...
                ['A location-scale run must not carry pop_lognu: downstream ' ...
                'consumers gate on that field existing and would report a ' ...
                'log-nu summary built from location-scale draws.']);

            bHat = mean(bdraw); bsHat = mean(bsdraw); ssdHat = mean(gro(:,2));
            fprintf(['[dynrel_trt LS recovery] b hat=%.3f (true %.2f)  ' ...
                'b_sigma hat=%.3f (true %.2f)  s_sd hat=%.3f (true %.2f)  ' ...
                'rho hat=%.3f (true %.2f)\n'], ...
                bHat, T.b, bsHat, T.b_sigma, ssdHat, T.s_sd, mean(rho_ps), T.rho);
            fprintf(['[dynrel_trt LS recovery] implied log-nu slope ' ...
                '2b-2b_sigma: recovered %.3f, true %.3f (of which %.3f is pure ' ...
                'amplitude)\n'], 2*bHat - 2*bsHat, 2*T.b - 2*T.b_sigma, 2*T.b);

            % (1) mean dimension slope recovered
            testCase.verifyLessThan(abs(bHat - T.b), 0.15, sprintf( ...
                'Recovered log-mean slope %.3f far from true %.2f.', bHat, T.b));

            % (2) DISCRIMINATING: the residual-SD slope. Sign first - a positive
            % recovered b_sigma would mean the fit put the residual on the wrong
            % side of the mean surface - then magnitude.
            testCase.verifyLessThan(bsHat, -0.05, sprintf( ...
                ['Recovered residual-SD slope %.3f is not clearly negative ' ...
                '(true %.2f). The planted design has error DECREASING with z ' ...
                'while amplitude INCREASES; getting this sign wrong means the ' ...
                'scale submodel has picked up the mean slope.'], bsHat, T.b_sigma));
            testCase.verifyLessThan(abs(bsHat - T.b_sigma), 0.25, sprintf( ...
                'Recovered residual-SD slope %.3f far from true %.2f.', ...
                bsHat, T.b_sigma));

            % (3) person residual-SD heterogeneity recovered in a band
            testCase.verifyGreaterThan(ssdHat, 0.20);
            testCase.verifyLessThan(ssdHat, 0.75);

            % (4) DECISIVE: the recovered z-conditioned two-facet G/D surface
            % covers the true surface at each grid z.
            surf = psyrat_rel_dynrel_trt_gamma_ls('alpha0',alpha0,'b',bdraw,...
                'logsig0',logsig0,'b_sigma',bsdraw,'sig_p',gro(:,1),...
                'sig_sd',gro(:,2),'sig_occ',s_occ,'sig_trl',s_trl,...
                'sig_trlxid',s_txp,'sig_occxid',s_oxp,'sig_trlxocc',s_txo,...
                'ndim',1,'z1',zgrid,'obs',T.ntrl,'nocc',T.nocc,...
                'reltype',reltype,'CI',ci);
            for k = 1:numel(zgrid)
                fprintf(['   z=%+d  G true %.3f in [%.3f, %.3f]   D true %.3f ' ...
                    'in [%.3f, %.3f]\n'], zgrid(k), trueG(k), surf.G.ll(k), ...
                    surf.G.ul(k), trueD(k), surf.D.ll(k), surf.D.ul(k));
                testCase.verifyGreaterThanOrEqual(trueG(k), surf.G.ll(k), sprintf( ...
                    'True G(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', ...
                    zgrid(k), trueG(k), surf.G.ll(k)));
                testCase.verifyLessThanOrEqual(trueG(k), surf.G.ul(k), sprintf( ...
                    'True G(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', ...
                    zgrid(k), trueG(k), surf.G.ul(k)));
                testCase.verifyGreaterThanOrEqual(trueD(k), surf.D.ll(k), sprintf( ...
                    'True D(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', ...
                    zgrid(k), trueD(k), surf.D.ll(k)));
                testCase.verifyLessThanOrEqual(trueD(k), surf.D.ul(k), sprintf( ...
                    'True D(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', ...
                    zgrid(k), trueD(k), surf.D.ul(k)));
            end

            % (5) DIAGNOSTIC, deliberately not asserted. The residual SURFACE is
            % the parameterization's interpretive payload; printing it keeps the
            % per-z recovery visible next to the coefficient, which increment 7
            % established must be read separately (its per-cell errors ran to
            % +30.5% while the contrast landed at 6.7%). No tolerance is gated
            % here because none has been derived rather than fitted to one run.
            fprintf('[dynrel_trt LS recovery] residual SD by z, true  = [%s]\n', ...
                localFmt(trueRes));
            fprintf('[dynrel_trt LS recovery] residual SD by z, recov = [%s]\n', ...
                localFmt(surf.sige_pt(:)'));
        end

        function testPairedLogNuFitAgreesOnThePersonBlock(testCase)
            % The live counterpart of the offline invariance proof: fit the SAME
            % simulated data under BOTH parameterizations and check they land in
            % the same place. Increments 4-7 each ran this; it is what produced the
            % structural log-nu cost ratio (2.9x / 2.4x / 2.3x / 2.48x so far).
            %
            % WHAT IS AND IS NOT EXACT HERE. Per ruling 28, the two complete model
            % classes are NOT exact reparameterizations of one another once trial
            % and occasion mean effects are present, so this is a CONSISTENCY check
            % with generous tolerances, not an identity test. What IS exact - and is
            % proven offline at AbsTol = 0 by
            % TestGammaDynrelTrtLsConversion/testSignalComponentsAreParameterizationInvariant
            % - is that the six observed-score SIGNAL components agree given the
            % same mean-submodel parameters. This test asks the different question
            % of whether two independent FITS recover those parameters compatibly.
            %
            % THE ASSERTION IS THE STRUCTURAL COST, not the coefficient. Substituting
            % v_p = const + 2*u_p - 2*w_p into the planted values gives
            %   Var(v_p) = 4*s_p^2 + 4*s_sd^2 - 8*rho*s_p*s_sd,
            % so the log-nu model must spend a MUCH larger person-dispersion SD to
            % describe the same data, because its single parameter carries the
            % amplitude term as well as the precision term. That is a property of
            % the parameterization rather than of any one design, and it has now
            % reproduced on four designs.
            T = localTruth();
            dataTbl = localSimulateGammaDynrelTrtLs(T);

            relLs = localRunDynrelTrtGammaLs(testCase, dataTbl, ...
                'recov_gamma_dynrel_trt_ls_pairA');
            relNu = localRunDynrelTrtGammaLogNu(testCase, dataTbl, ...
                'recov_gamma_dynrel_trt_ls_pairB');

            relLs = psyrat_checkconv(relLs);
            relNu = psyrat_checkconv(relNu);
            testCase.assumeEqual(relLs.out.conv.converged, 1, ...
                'Location-scale fit did not converge; paired check skipped.');
            testCase.assumeEqual(relNu.out.conv.converged, 1, ...
                'Log-nu fit did not converge; paired check skipped.');

            groLs = cell2mat(relLs.out.gro_sds{1});   % [s_p, s_sd]
            groNu = cell2mat(relNu.out.gro_sds{1});   % [s_p, s_v]
            s_p_ls = mean(groLs(:,1)); s_sd = mean(groLs(:,2));
            s_p_nu = mean(groNu(:,1)); s_v  = mean(groNu(:,2));
            cholLs = cell2mat(relLs.out.chol_corrmat{1});
            rho_ls = mean(cholLs(:,2,1));

            % Predicted log-nu dispersion SD from the location-scale fit's own
            % person block, via v_p = const + 2*u_p - 2*w_p.
            predSv = sqrt(4*s_p_ls^2 + 4*s_sd^2 - 8*rho_ls*s_p_ls*s_sd);

            fprintf(['\n[dynrel_trt LS paired] s_p: LS %.3f vs log-nu %.3f  |  ' ...
                'dispersion SD: s_sd %.3f vs s_v %.3f (ratio %.2fx)\n'], ...
                s_p_ls, s_p_nu, s_sd, s_v, s_v/s_sd);
            fprintf(['[dynrel_trt LS paired] predicted s_v from the LS person ' ...
                'block = %.3f, measured %.3f (%.1f%%)\n'], ...
                predSv, s_v, 100*(s_v/predSv - 1));

            % (1) The person MEAN SD is a log-mean submodel quantity, so the two
            % fits must agree on it closely - this is the live counterpart of the
            % byte-identical mean submodel in the golden diff.
            testCase.verifyLessThan(abs(s_p_ls - s_p_nu)/s_p_ls, 0.30, sprintf( ...
                ['Person log-mean SD disagrees across parameterizations ' ...
                '(%.3f vs %.3f). The mean submodel is byte-identical between ' ...
                'the two emitted models, so a large gap here means one fit is ' ...
                'reading the wrong slot.'], s_p_ls, s_p_nu));

            % (2) THE STRUCTURAL COST. Direction is the assertion; the magnitude is
            % printed. Increments 4-7 measured 2.9x / 2.4x / 2.3x / 2.48x.
            testCase.verifyGreaterThan(s_v, 1.5*s_sd, sprintf( ...
                ['The log-nu fit should need a substantially larger person ' ...
                'dispersion SD (measured %.3f vs %.3f); its parameter carries ' ...
                'the amplitude term as well as the precision term.'], s_v, s_sd));

            % (3) The algebraic prediction, as a consistency check with a generous
            % band - trial/occasion mean effects break the exact correspondence
            % (ruling 28), and increments 5/7 recorded the measured value running
            % 10-20% above prediction on the residual side.
            testCase.verifyLessThan(abs(s_v - predSv)/predSv, 0.40, sprintf( ...
                ['Measured s_v %.3f is far from the %.3f predicted by the ' ...
                'location-scale person block.'], s_v, predSv));
        end
    end
end

% =========================================================================
function T = localTruth()
% Planted truth on the log expected-score scale. See the class header for why
% b_sigma = -0.20 rather than the degenerate b_sigma = b.
T = struct();
T.alpha   = log(5.5);   T.b       = 0.30;   % log-mean intercept + dimension slope
T.logsig0 = log(3.0);   T.b_sigma = -0.20;  % log residual-SD intercept + slope
T.sP  = 0.35;   % person
T.sO  = 0.20;   % occasion
T.sT  = 0.20;   % trial
T.sPO = 0.12;   % person x occasion
T.sPT = 0.18;   % person x trial
T.sOT = 0.08;   % trial x occasion
T.s_sd = 0.40;  % person log residual-SD SD
T.rho  = -0.30; % person mean <-> log residual-SD correlation
T.nsub = 60; T.nocc = 4; T.ntrl = 12;   % 2880 rows
end

% =========================================================================
function tbl = localSimulateGammaDynrelTrtLs(T)
% Balanced crossed Persons x Occasions x Trials long table with a between-person
% standardized dimension z_p driving BOTH the log-mean and the log RESIDUAL SD.
% Six independent log-mean random effects, a person (mean, log-sigma) block
% correlated at rho, and an EXACT scaled chi-square observation with derived df.
rng(20260805);
u_p  = T.sP  * randn(T.nsub, 1);
u_o  = T.sO  * randn(T.nocc, 1);
u_t  = T.sT  * randn(T.ntrl, 1);
u_po = T.sPO * randn(T.nsub, T.nocc);
u_pt = T.sPT * randn(T.nsub, T.ntrl);
u_ot = T.sOT * randn(T.ntrl, T.nocc);

% person log residual-SD effect, correlated with u_p at rho (share the person z1)
z1 = u_p / T.sP;                     % standard normal (u_p = sP*z1)
z2 = randn(T.nsub, 1);
w_p = T.s_sd * (T.rho * z1 + sqrt(1 - T.rho^2) * z2);

% standardized subject-level covariate (exact mean 0, SD 1 over subjects)
raw = randn(T.nsub, 1);
zp = (raw - mean(raw)) / std(raw);

% each person's residual SD at their own z - NO facet terms, matching the model
sigma_p = exp(T.logsig0 + T.b_sigma * zp + w_p);

nrow = T.nsub * T.nocc * T.ntrl;
ids = cell(nrow,1); time = zeros(nrow,1);
meas = zeros(nrow,1); dim1 = zeros(nrow,1);
r = 0;
for p = 1:T.nsub
    for o = 1:T.nocc
        for t = 1:T.ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            mu = exp(T.alpha + T.b*zp(p) + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p,o) + u_pt(p,t) + u_ot(t,o));
            nu = 2 * (mu / sigma_p(p))^2;   % derived df, non-integer
            % y = (mu/nu)*chi2(nu) == Gamma(nu/2, 2*mu/nu); E[y] = mu and
            % Var(y) = 2*mu^2/nu = sigma_p^2 by construction.
            meas(r) = localGammaRand(nu/2) * (2 * mu / nu);
            dim1(r) = zp(p);
        end
    end
end
tbl = table(ids, meas, time, dim1, 'VariableNames', {'id','meas','time','dim1'});
end

% =========================================================================
function x = localGammaRand(shape)
% Exact Gamma(shape, scale = 1) draw by the Marsaglia-Tsang (2000) method, using
% only base-MATLAB randn/rand so the test runs on a CI image installed with
% --products=MATLAB. Statistics Toolbox gamrnd would be simpler and is not
% available there.
%
% The method is defined for shape >= 1; for shape < 1 it uses the standard boost
% Gamma(a) = Gamma(a+1) * U^(1/a). The rejection loop is guaranteed to terminate
% (acceptance probability is bounded well away from zero for every shape), and
% the v <= 0 guard rejects the cube-root branch where the transform is invalid.
if shape < 1
    x = localGammaRand(shape + 1) * rand^(1/shape);
    return;
end
d = shape - 1/3;
c = 1 / sqrt(9*d);
while true
    v = -1;
    while v <= 0
        z = randn;
        v = (1 + c*z)^3;
    end
    u = rand;
    if log(u) < 0.5*z^2 + d - d*v + d*log(v)
        x = d * v;
        return;
    end
end
end

% =========================================================================
function s = localFmt(v)
s = strjoin(arrayfun(@(x) sprintf('%.3f', x), v(:)', 'UniformOutput', false), ', ');
end

% =========================================================================
function rel = localRunDynrelTrtGammaLs(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'gammascale', 2, ...
    'dynrel', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

% =========================================================================
function rel = localRunDynrelTrtGammaLogNu(testCase, dataTbl, tag)
% The same run under the legacy log-nu parameterization, for the paired
% cross-parameterization check. Identical data, seed and sampler settings, so
% the only difference is the scale submodel. gammascale=1 must be EXPLICIT:
% this helper originally omitted it back when log-nu was the default, but
% ruling 33 (2026-08-07) flipped the movable analyses to location-scale, so an
% omitted gammascale now silently refits the location-scale model.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'dynrel', 2, ...
    'gammascale', 1, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
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
