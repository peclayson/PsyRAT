classdef TestGammaDynrelTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the TWO-facet (trial + occasion) Gamma /
    % scaled-chi-square DYNAMIC/conditional reliability design (analysis 14,
    % family='gamma'; §20). Crossed Persons x Occasions x Trials data are simulated
    % from a log-linked mean (six independent random effects) and a person-varying
    % scaled-chi-square dispersion, BOTH conditioned on a standardized subject-level
    % dimension z:
    %
    %   log mu_pot = alpha + b*z_p    + u_p + u_o + u_t + u_po + u_pt + u_ot
    %   log nu_p   = beta  + b_nu*z_p + v_p      (v_p corr rho with u_p)
    %
    % The fit must recover the dimension slopes (b, b_nu) and the person-dispersion
    % SD (s_v), and the z-conditioned two-facet coefficient surface (reltype 3, both
    % facets random) computed from the recovered draws (psyrat_rel_dynrel_trt_gamma)
    % must cover the closed-form true surface at z in {-1, 0, +1}.
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
                'recov_gamma_dynrel_trt');
        end
    end

    methods (Test)
        function testTwoFacetDynrelRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   b    = 0.30;   % log-mean intercept + dimension slope
            sP  = 0.35;         % person
            sO  = 0.20;         % occasion
            sT  = 0.20;         % trial
            sPO = 0.12;         % person x occasion
            sPT = 0.18;         % person x trial
            sOT = 0.08;         % trial x occasion
            beta = log(18);     b_nu = 0.30;   % log-nu intercept + dimension slope
            s_v  = 0.40;        rho  = -0.30;  % person dispersion SD + mean<->log-nu corr
            cov_pv = rho * sP * s_v;

            nsub = 60; nocc = 4; ntrl = 12;
            zgrid = [-1 0 1]; reltype = 3; ci = 0.95;

            dataTbl = localSimulateGammaDynrelTrt(nsub, nocc, ntrl, alpha, b, ...
                sP, sO, sT, sPO, sPT, sOT, beta, b_nu, s_v, rho);

            % ---- Closed-form true G(z)/D(z) at each grid z (reltype 3, at the design
            % trial/occasion counts). Argument order matches psyrat_gamma_extract_trt:
            % (alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu, s_v, cov_pv). ----
            trueG = zeros(1,numel(zgrid)); trueD = zeros(1,numel(zgrid));
            for k = 1:numel(zgrid)
                z = zgrid(k);
                vc = psyrat_gamma_varcomps_trt(alpha + b*z, sP, sO, sT, sPO, sPT, sOT, ...
                    exp(beta + b_nu*z), s_v, cov_pv);
                fa = {'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2),...
                      'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2),...
                      'txo',sqrt(vc.sigma_ot2)};
                [~,trueG(k)] = psyrat_rel_trt('gcoeff',2,'reltype',reltype,...
                    'bp',sqrt(vc.sigma_p2),fa{:},'err',sqrt(vc.sigma_res2),...
                    'obs',ntrl,'nocc',nocc,'CI',ci);
                [~,trueD(k)] = psyrat_rel_trt('gcoeff',1,'reltype',reltype,...
                    'bp',sqrt(vc.sigma_p2),fa{:},'err',sqrt(vc.sigma_res2),...
                    'obs',ntrl,'nocc',nocc,'CI',ci);
            end

            % ---- Fit ----
            rel = localRunDynrelTrtGamma(testCase, dataTbl, 'recov_gamma_dynrel_trt');
            testCase.verifyEqual(rel.analysis, 'ic_dynrel_trt');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Recovered draws (single stratum) ----
            alpha0 = rel.out.mu(:,1);
            lognu0 = rel.out.pop_lognu(:,1);
            gro    = cell2mat(rel.out.gro_sds{1});   % draws x 2 [s_p, s_v]
            bdraw  = cell2mat(rel.out.b{1});
            bnudraw= cell2mat(rel.out.b_nu{1});
            chol   = cell2mat(rel.out.chol_corrmat{1});
            rho_pv = chol(:,2,1);
            s_occ  = rel.out.sig_occ(:,1);
            s_trl  = rel.out.sig_trl(:,1);
            s_txp  = rel.out.sig_trlxid(:,1);
            s_oxp  = rel.out.sig_occxid(:,1);
            s_txo  = rel.out.sig_trlxocc(:,1);

            bHat = mean(bdraw); bnuHat = mean(bnudraw); svHat = mean(gro(:,2));
            fprintf(['[dynrel_trt recovery] b hat=%.3f (true %.2f)  b_nu hat=%.3f (true %.2f)  ',...
                's_v hat=%.3f (true %.2f)  rho hat=%.3f (true %.2f)\n'], ...
                bHat, b, bnuHat, b_nu, svHat, s_v, mean(rho_pv), rho);

            % (1) mean dimension slope recovered
            testCase.verifyLessThan(abs(bHat - b), 0.15, sprintf( ...
                'Recovered log-mean slope %.3f far from true %.2f.', bHat, b));
            % (2) log-nu dimension slope: right sign and rough magnitude
            testCase.verifyGreaterThan(bnuHat, 0.05, sprintf( ...
                'Recovered log-nu slope %.3f not clearly positive (true %.2f).', bnuHat, b_nu));
            testCase.verifyLessThan(abs(bnuHat - b_nu), 0.35);
            % (3) person-dispersion SD recovered
            testCase.verifyGreaterThan(svHat, 0.20);
            testCase.verifyLessThan(svHat, 0.75);

            % (4) DECISIVE: the recovered z-conditioned two-facet G/D surface covers
            % the true surface at each grid z.
            rel_surf = psyrat_rel_dynrel_trt_gamma('alpha0',alpha0,'b',bdraw,...
                'lognu0',lognu0,'b_nu',bnudraw,'sig_p',gro(:,1),'sig_v',gro(:,2),...
                'rho_pv',rho_pv,'sig_occ',s_occ,'sig_trl',s_trl,'sig_trlxid',s_txp,...
                'sig_occxid',s_oxp,'sig_trlxocc',s_txo,'ndim',1,'z1',zgrid,...
                'obs',ntrl,'nocc',nocc,'reltype',reltype,'CI',ci);
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
        end
    end
end

function tbl = localSimulateGammaDynrelTrt(nsub, nocc, ntrl, alpha, b, ...
    sP, sO, sT, sPO, sPT, sOT, beta, b_nu, s_v, rho)
% Balanced crossed Persons x Occasions x Trials long table with a between-person
% standardized dimension z_p driving BOTH the log-mean and the log-nu. Six
% independent log-mean random effects (person, occasion, trial + three two-way
% interactions), a person (mean, log-nu) block correlated at rho, and a
% continuous-df scaled chi-square observation Y = gamrnd(nu_p/2, 2*mu/nu_p).
rng(20260710);
u_p  = sP  * randn(nsub, 1);
u_o  = sO  * randn(nocc, 1);
u_t  = sT  * randn(ntrl, 1);
u_po = sPO * randn(nsub, nocc);
u_pt = sPT * randn(nsub, ntrl);
u_ot = sOT * randn(ntrl, nocc);
% person log-nu effect, correlated with u_p at rho (share the person z1)
z1 = u_p / sP;                       % standard normal (u_p = sP*z1)
z2 = randn(nsub, 1);
v_p = s_v * (rho * z1 + sqrt(1 - rho^2) * z2);
% standardized subject-level covariate (exact mean 0, SD 1 over subjects)
raw = randn(nsub, 1);
zp = (raw - mean(raw)) / std(raw);
nu_p = exp(beta + b_nu * zp + v_p);

nrow = nsub * nocc * ntrl;
ids = cell(nrow,1); time = zeros(nrow,1); meas = zeros(nrow,1); dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            mu = exp(alpha + b*zp(p) + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p,o) + u_pt(p,t) + u_ot(t,o));
            meas(r) = gamrnd(nu_p(p) / 2, 2 * mu / nu_p(p));
            dim1(r) = zp(p);
        end
    end
end
tbl = table(ids, meas, time, dim1, 'VariableNames', {'id','meas','time','dim1'});
end

function rel = localRunDynrelTrtGamma(testCase, dataTbl, tag)
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
