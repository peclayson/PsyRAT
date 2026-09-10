classdef TestGammaDynrelRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the ONE-facet Gamma / scaled-chi-square
    % DYNAMIC/conditional (dimensional) reliability design (analysis 11,
    % family='gamma'; §20). Data are simulated from a log-linked mean and a
    % person-varying scaled-chi-square dispersion, BOTH conditioned on a
    % standardized subject-level dimension z:
    %
    %   log mu_pt = alpha + b*z_p    + u_p + u_t
    %   log nu_p  = beta  + b_nu*z_p + v_p        (v_p corr rho with u_p)
    %
    % The fit must recover the dimension slopes (b, b_nu) and the person-dispersion
    % SD (s_v), and the z-conditioned generalizability/dependability surface
    % G(z)/D(z) computed from the recovered draws (psyrat_rel_dynrel_gamma) must
    % cover the closed-form true surface at z in {-1, 0, +1}. A well-powered design
    % (many persons x many trials) keeps the dispersion side identified (§20 caveat:
    % the mean-coupled residual makes b_nu data-hungry).
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()), ...
                'recov_gamma_dynrel');
        end
    end

    methods (Test)
        function testOneFacetDynrelRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   % log-mean intercept at z = 0
            b     = 0.30;       % log-mean dimension slope
            s_p   = 0.35;       % person log-mean SD
            s_i   = 0.20;       % trial  log-mean SD
            beta  = log(18);    % log-nu intercept at z = 0 (nu_pop(0) = 18)
            b_nu  = 0.40;       % log-nu dimension slope (dispersion rises with z)
            s_v   = 0.50;       % person log-nu SD
            rho   = -0.30;      % person mean <-> log-nu correlation
            cov_pv = rho * s_p * s_v;

            nsub = 120; ntrl = 40;
            [dataTbl, zp] = localSimulateGammaDynrel(nsub, ntrl, alpha, b, s_p, s_i, ...
                beta, b_nu, s_v, rho);

            % z-grid for the recovery comparison (well inside the observed z range)
            zgrid = [-1 0 1];
            ci = 0.95;

            % ---- Closed-form true G(z)/D(z) at each grid z (the same converter the
            % production surface uses, at the true parameters) ----
            trueG = zeros(1,numel(zgrid)); trueD = zeros(1,numel(zgrid));
            for k = 1:numel(zgrid)
                z = zgrid(k);
                vc = psyrat_gamma_varcomps(alpha + b*z, s_p, s_i, ...
                    exp(beta + b_nu*z), s_v, cov_pv);
                [~,trueG(k)] = psyrat_rel_sing('gcoeff',2,'metric','global',...
                    'bp',sqrt(vc.sigma_p2),'wp',sqrt(vc.sigma_pi_e2+vc.sigma_i2),...
                    'i',sqrt(vc.sigma_i2),'obs',ntrl,'CI',ci);
                [~,trueD(k)] = psyrat_rel_sing('gcoeff',1,'metric','global',...
                    'bp',sqrt(vc.sigma_p2),'wp',sqrt(vc.sigma_pi_e2+vc.sigma_i2),...
                    'i',sqrt(vc.sigma_i2),'obs',ntrl,'CI',ci);
            end

            % ---- Fit ----
            rel = localRunDynrelGamma(testCase, dataTbl, 'recov_gamma_dynrel');
            testCase.verifyEqual(rel.analysis, 'ic_dynrel');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Recovered draws (single stratum) ----
            alpha0 = rel.out.mu(:,1);
            lognu0 = rel.out.pop_lognu(:,1);
            gro    = cell2mat(rel.out.gro_sds{1});   % draws x 2 [s_p, s_v]
            bdraw  = cell2mat(rel.out.b{1});         % draws x 1 (KDIM=1)
            bnudraw= cell2mat(rel.out.b_nu{1});      % draws x 1
            sig_i  = rel.out.sig_trl(:,1);
            chol   = cell2mat(rel.out.chol_corrmat{1});
            rho_pv = chol(:,2,1);

            bHat   = mean(bdraw);
            bnuHat = mean(bnudraw);
            svHat  = mean(gro(:,2));
            fprintf(['[dynrel recovery] b hat=%.3f (true %.2f)  b_nu hat=%.3f (true %.2f)  ',...
                's_v hat=%.3f (true %.2f)  rho hat=%.3f (true %.2f)\n'], ...
                bHat, b, bnuHat, b_nu, svHat, s_v, mean(rho_pv), rho);

            % (1) mean dimension slope recovered (tight - the mean side is well identified)
            testCase.verifyLessThan(abs(bHat - b), 0.15, sprintf( ...
                'Recovered log-mean slope %.3f far from true %.2f.', bHat, b));
            % (2) log-nu dimension slope recovered with the right sign and rough size
            testCase.verifyGreaterThan(bnuHat, 0.10, sprintf( ...
                'Recovered log-nu slope %.3f not clearly positive (true %.2f).', bnuHat, b_nu));
            testCase.verifyLessThan(abs(bnuHat - b_nu), 0.35, sprintf( ...
                'Recovered log-nu slope %.3f far from true %.2f.', bnuHat, b_nu));
            % (3) person-dispersion SD recovered (the param #1 discards)
            testCase.verifyGreaterThan(svHat, 0.30);
            testCase.verifyLessThan(svHat, 0.80);

            % (4) DECISIVE: the recovered z-conditioned G/D surface covers the true
            % surface at each grid z (folds in b, b_nu, s_p, s_v, rho jointly).
            rel_surf = psyrat_rel_dynrel_gamma('alpha0',alpha0,'b',bdraw,...
                'lognu0',lognu0,'b_nu',bnudraw,'sig_p',gro(:,1),'sig_v',gro(:,2),...
                'rho_pv',rho_pv,'sig_i',sig_i,'ndim',1,'z1',zgrid,'obs',ntrl,'CI',ci);
            for k = 1:numel(zgrid)
                fprintf('   z=%+d  G true %.3f in [%.3f, %.3f]   D true %.3f in [%.3f, %.3f]\n',...
                    zgrid(k), trueG(k), rel_surf.G.ll(k), rel_surf.G.ul(k), ...
                    trueD(k), rel_surf.D.ll(k), rel_surf.D.ul(k));
                testCase.verifyGreaterThanOrEqual(trueG(k), rel_surf.G.ll(k), sprintf( ...
                    'True G(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', zgrid(k), trueG(k), rel_surf.G.ll(k)));
                testCase.verifyLessThanOrEqual(trueG(k), rel_surf.G.ul(k), sprintf( ...
                    'True G(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', zgrid(k), trueG(k), rel_surf.G.ul(k)));
                testCase.verifyGreaterThanOrEqual(trueD(k), rel_surf.D.ll(k), sprintf( ...
                    'True D(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', zgrid(k), trueD(k), rel_surf.D.ll(k)));
                testCase.verifyLessThanOrEqual(trueD(k), rel_surf.D.ul(k), sprintf( ...
                    'True D(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', zgrid(k), trueD(k), rel_surf.D.ul(k)));
            end
            % (5) monotonicity: a positive b_nu (dispersion rises with z, nu rises,
            % 2*mu^2/nu residual relatively lower) yields G increasing across z.
            testCase.verifyTrue(rel_surf.G.pt(end) > rel_surf.G.pt(1), ...
                'Recovered G(z) should rise across z for a positive log-nu slope.');
        end
    end
end

function [tbl, zp] = localSimulateGammaDynrel(nsub, ntrl, alpha, b, s_p, s_i, ...
    beta, b_nu, s_v, rho)
% One-facet Persons x Trials long table with a between-person standardized
% dimension z_p driving BOTH the log-mean and the log-nu. The person (mean,
% log-nu) block is bivariate normal with correlation rho; the continuous-df
% scaled chi-square Y = gamrnd(nu_p/2, 2*mu/nu_p) has E[Y]=mu, Var[Y]=2*mu^2/nu_p.
% dim1 is stored raw = z_p (already standardized), so the pipeline's subject-level
% z-scoring reproduces z_p (mean 0, SD ~1); the tiny residual re-standardization is
% harmless (the truth is compared on the recovered scale, at matched z values).
rng(20260710);
u_t = s_i * randn(ntrl, 1);
z1 = randn(nsub, 1); z2 = randn(nsub, 1);
u_p = s_p * z1;
v_p = s_v * (rho * z1 + sqrt(1 - rho^2) * z2);
% standardized subject-level covariate (exact mean 0, SD 1 over subjects)
raw = randn(nsub, 1);
zp = (raw - mean(raw)) / std(raw);
nu_p = exp(beta + b_nu * zp + v_p);

nrow = nsub * ntrl;
ids = cell(nrow,1); meas = zeros(nrow,1); dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    for t = 1:ntrl
        r = r + 1;
        ids{r} = sprintf('S%03d', p);
        mu = exp(alpha + b * zp(p) + u_p(p) + u_t(t));
        meas(r) = gamrnd(nu_p(p) / 2, 2 * mu / nu_p(p));
        dim1(r) = zp(p);
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id', 'meas', 'dim1'});
end

function rel = localRunDynrelGamma(testCase, dataTbl, tag)
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
