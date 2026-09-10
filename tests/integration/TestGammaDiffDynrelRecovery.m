classdef TestGammaDiffDynrelRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the one-facet GROUP-LEVEL DYNAMIC
    % DIFFERENCE-score design under the Gamma / scaled-chi-square family
    % (analysis 12, family='gamma', non-concurrent). Two component cells are
    % simulated jointly from log-linked means and per-event, person-varying
    % scaled-chi-square dispersions, BOTH conditioned on a standardized
    % subject-level dimension z:
    %
    %   log mu_pet = b[e] + b_dim[e]*z_p    + u_p[p,e] + u_t[t,e]
    %   log nu_pe  = b_nu[e] + b_nu_dim[e]*z_p + v_p[p,e]
    %
    % with a 4-D person block over (mean_1, mean_2, lognu_1, lognu_2) carrying the
    % cross-event mean correlation and the within-event mean<->log-nu couplings,
    % and a 2-D trial block carrying the cross-event trial correlation.
    %
    % The trial facet is simulated on PsyRAT's own definition: t indexes the
    % within-person, within-event trial POSITION (psyrat_withinrunindex), so
    % position t of event 1 pairs with position t of event 2. See S3 in
    % PSYRAT_AUDIT_FINDINGS.md - the owner-supplied reference instead uses global
    % trial identity, an open convention question. Simulating on the toolbox's
    % definition keeps the generative model matched to the estimand under test.
    %
    % DELIBERATELY UNBALANCED per-event trial counts (25 vs 60), because
    % difference designs routinely are (the reference ERP design is 49 error vs
    % 341 correct trials) and the per-event n' path is a specific thing to prove.
    %
    % The DECISIVE assertion is (5): the z-conditioned G(z)/D(z) surface computed
    % from the recovered draws must cover the closed-form true surface at every
    % grid z. That folds every parameter in jointly, so it is a stronger check
    % than any single marginal.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent),
    % and R-hat-gated via assume so a non-converged fit SKIPS rather than
    % false-passing (the convention adopted for the copula recoveries).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()), ...
                'recov_gamma_diff_dynrel');
        end
    end

    methods (Test)
        function testGroupLevelDiffDynrelRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score / log-nu scales) ----
            b       = [log(5.5)  log(4.5)];  % per-event log-mean intercepts at z=0
            b_dim   = [0.30      0.10];      % per-event log-mean dimension slopes
            b_nu    = [log(18)   log(20)];   % per-event log-nu intercepts at z=0
            b_nudim = [0.35      0.30];      % per-event log-nu dimension slopes
            sd_id   = [0.35 0.32 0.45 0.40]; % [mean_1 mean_2 lognu_1 lognu_2]
            sd_trl  = [0.20 0.18];
            rho_p   = 0.60;   % cross-event person-MEAN correlation (drives the difference)
            rho_pv1 = -0.30;  % within-event mean<->log-nu, event 1 (estimand #2)
            rho_pv2 = -0.25;  % within-event mean<->log-nu, event 2
            rho_v   = 0.50;   % cross-event person log-nu correlation (not used downstream)
            rho_i   = 0.30;   % cross-event TRIAL correlation

            nsub = 120; n1 = 25; n2 = 60;   % unbalanced on purpose
            dataTbl = localSimulateGammaDiffDynrel(nsub, n1, n2, b, b_dim, ...
                b_nu, b_nudim, sd_id, sd_trl, rho_p, rho_pv1, rho_pv2, rho_v, rho_i);

            zgrid = [-1 0 1];
            ci = 0.95;
            obs = [n1 n2];

            % ---- Closed-form true G(z)/D(z), via the same converters the
            % production surface uses, evaluated at the TRUE parameters ----
            [trueG, trueD] = localTrueSurface(zgrid, obs, ci, b, b_dim, b_nu, ...
                b_nudim, sd_id, sd_trl, rho_p, rho_pv1, rho_pv2, rho_i);

            % ---- Fit ----
            rel = localRunDiffDynrelGamma(testCase, dataTbl, 'recov_gamma_diff_dynrel');
            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Gate on convergence via assume: a non-converged fit is filtered
            % (Incomplete), never a false PASS. Same R-hat>=1.1 + n_eff criteria as
            % the production psyrat_checkconv path.
            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma diff dynrel recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma dynamic ', ...
                'difference fit non-converged (max R-hat=%.4f >= 1.1 / n_eff gate); ', ...
                'recovery not validated in this environment.'], maxRhat));

            % ---- Recovered draws (single stratum) ----
            bb     = cell2mat(rel.out.b{1});         % draws x 2
            bnu    = cell2mat(rel.out.b_nu{1});      % draws x 2
            bd     = cell2mat(rel.out.b_dim{1});     % draws x 2 (KDIM=2)
            bnd    = cell2mat(rel.out.b_nu_dim{1});  % draws x 2
            sdid   = cell2mat(rel.out.sd_id{1});     % draws x 4
            sdtrl  = cell2mat(rel.out.sd_trl{1});    % draws x 2
            corid  = cell2mat(rel.out.cor_id{1});    % draws x 4 x 4
            cori   = cell2mat(rel.out.cor_i{1});     % draws x 1

            fprintf(['   b_dim hat=[%.3f %.3f] (true [%.2f %.2f])  ', ...
                'b_nu_dim hat=[%.3f %.3f] (true [%.2f %.2f])  cor_p hat=%.3f (true %.2f)\n'], ...
                mean(bd(:,1)), mean(bd(:,2)), b_dim(1), b_dim(2), ...
                mean(bnd(:,1)), mean(bnd(:,2)), b_nudim(1), b_nudim(2), ...
                mean(corid(:,1,2)), rho_p);

            % (1) per-event n' reached REL from the unbalanced design
            testCase.verifyEqual(rel.out.ntrials_ev(:,1)', [n1 n2], ...
                'Per-event trial counts were not recovered from the unbalanced design.');

            % (2) log-MEAN dimension slopes: well identified, so assert tightly
            testCase.verifyLessThan(abs(mean(bd(:,1)) - b_dim(1)), 0.12, sprintf( ...
                'Recovered event-1 log-mean slope %.3f far from true %.2f.', ...
                mean(bd(:,1)), b_dim(1)));
            testCase.verifyLessThan(abs(mean(bd(:,2)) - b_dim(2)), 0.12, sprintf( ...
                'Recovered event-2 log-mean slope %.3f far from true %.2f.', ...
                mean(bd(:,2)), b_dim(2)));

            % (3) log-NU dimension slopes: the dispersion side is data-hungry under
            % the mean-coupled residual (SCIENTIFIC_FORMULA_AUDIT.md section 20), and
            % event 1 carries only n1 trials, so these shrink toward 0. Assert the
            % SIGN and a generous magnitude only - a tight bound here would encode a
            % known identification limit as a regression trip-wire.
            testCase.verifyGreaterThan(mean(bnd(:,1)), 0.08, sprintf( ...
                'Recovered event-1 log-nu slope %.3f not clearly positive (true %.2f).', ...
                mean(bnd(:,1)), b_nudim(1)));
            testCase.verifyGreaterThan(mean(bnd(:,2)), 0.08, sprintf( ...
                'Recovered event-2 log-nu slope %.3f not clearly positive (true %.2f).', ...
                mean(bnd(:,2)), b_nudim(2)));
            testCase.verifyLessThan(abs(mean(bnd(:,1)) - b_nudim(1)), 0.30);
            testCase.verifyLessThan(abs(mean(bnd(:,2)) - b_nudim(2)), 0.30);

            % (4) the cross-event person-MEAN correlation drives the whole
            % difference numerator, so it must be recovered well
            testCase.verifyLessThan(abs(mean(corid(:,1,2)) - rho_p), 0.15, sprintf( ...
                'Recovered cross-event person-mean correlation %.3f far from true %.2f.', ...
                mean(corid(:,1,2)), rho_p));

            % (5) DECISIVE: the recovered z-conditioned surface covers the true
            % surface at each grid z (folds in every parameter jointly).
            rs = psyrat_rel_diffdynrel_gamma('b',bb,'b_nu',bnu,'b_dim',bd,...
                'b_nu_dim',bnd,'sd_id',sdid,'sd_trl',sdtrl,'cor_id',corid,...
                'cor_i',cori,'ndim',1,'z1',zgrid,'obs',obs,'CI',ci);
            for k = 1:numel(zgrid)
                fprintf('   z=%+d  G true %.3f in [%.3f, %.3f]   D true %.3f in [%.3f, %.3f]\n', ...
                    zgrid(k), trueG(k), rs.G.ll(k), rs.G.ul(k), ...
                    trueD(k), rs.D.ll(k), rs.D.ul(k));
                testCase.verifyGreaterThanOrEqual(trueG(k), rs.G.ll(k), sprintf( ...
                    'True G(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', ...
                    zgrid(k), trueG(k), rs.G.ll(k)));
                testCase.verifyLessThanOrEqual(trueG(k), rs.G.ul(k), sprintf( ...
                    'True G(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', ...
                    zgrid(k), trueG(k), rs.G.ul(k)));
                testCase.verifyGreaterThanOrEqual(trueD(k), rs.D.ll(k), sprintf( ...
                    'True D(z=%+d)=%.3f below recovered 95%% CI lower %.3f.', ...
                    zgrid(k), trueD(k), rs.D.ll(k)));
                testCase.verifyLessThanOrEqual(trueD(k), rs.D.ul(k), sprintf( ...
                    'True D(z=%+d)=%.3f above recovered 95%% CI upper %.3f.', ...
                    zgrid(k), trueD(k), rs.D.ul(k)));
            end

            % (6) positive log-nu slopes raise nu with z, lowering the mean-coupled
            % residual, so the recovered G(z) must increase across z.
            testCase.verifyTrue(rs.G.pt(end) > rs.G.pt(1), ...
                'Recovered G(z) should rise across z for positive log-nu slopes.');
        end
    end
end

function tbl = localSimulateGammaDiffDynrel(nsub, n1, n2, b, b_dim, b_nu, ...
    b_nudim, sd_id, sd_trl, rho_p, rho_pv1, rho_pv2, rho_v, rho_i)
% Two-event Persons x Trials long table. The person block is 4-D over
% (mean_1, mean_2, lognu_1, lognu_2); the trial block is 2-D over events and is
% indexed by within-person, within-event trial POSITION, matching PsyRAT's
% psyrat_withinrunindex trial facet. Y = gamrnd(nu/2, 2*mu/nu) has E[Y]=mu and
% Var[Y]=2*mu^2/nu, i.e. the scaled chi-square the model assumes.
% dim1 is stored raw = z_p (already standardized), so the pipeline's subject-level
% z-scoring reproduces z_p.
rng(20260720);

R = eye(4);
R(1,2) = rho_p;   R(2,1) = rho_p;
R(1,3) = rho_pv1; R(3,1) = rho_pv1;
R(2,4) = rho_pv2; R(4,2) = rho_pv2;
R(3,4) = rho_v;   R(4,3) = rho_v;
Up  = (chol(R,'lower') * randn(4,nsub))';    % nsub x 4
u_p = Up(:,1:2) .* sd_id(1:2);
v_p = Up(:,3:4) .* sd_id(3:4);

NTRL = max(n1,n2);
Ut   = (chol([1 rho_i; rho_i 1],'lower') * randn(2,NTRL))';
u_t  = Ut .* sd_trl;

raw = randn(nsub,1);
zp  = (raw - mean(raw)) / std(raw);

nper = [n1 n2];
enames = {'e1','e2'};
rows = nsub * (n1 + n2);
ids = cell(rows,1); evc = cell(rows,1);
meas = zeros(rows,1); dim1 = zeros(rows,1);
r = 0;
for p = 1:nsub
    for e = 1:2
        for t = 1:nper(e)
            r = r + 1;
            ids{r} = sprintf('S%03d', p);
            evc{r} = enames{e};
            mu = exp(b(e)    + b_dim(e)*zp(p)   + u_p(p,e) + u_t(t,e));
            nu = exp(b_nu(e) + b_nudim(e)*zp(p) + v_p(p,e));
            meas(r) = gamrnd(nu/2, 2*mu/nu);
            dim1(r) = zp(p);
        end
    end
end
tbl = table(ids, evc, meas, dim1, ...
    'VariableNames', {'id','event','meas','dim1'});
end

function [trueG, trueD] = localTrueSurface(zgrid, obs, ci, b, b_dim, b_nu, ...
    b_nudim, sd_id, sd_trl, rho_p, rho_pv1, rho_pv2, rho_i)
% Closed-form true G(z)/D(z): the production converters evaluated at the TRUE
% parameters, assembled exactly as psyrat_rel_diffdynrel_gamma does per z.
sp = sd_id(1:2); sv = sd_id(3:4); si = sd_trl;
cpv = [rho_pv1*sp(1)*sv(1), rho_pv2*sp(2)*sv(2)];
trueG = zeros(1,numel(zgrid));
trueD = zeros(1,numel(zgrid));
for k = 1:numel(zgrid)
    z = zgrid(k);
    aA = b(1) + b_dim(1)*z;
    aB = b(2) + b_dim(2)*z;
    vA = psyrat_gamma_varcomps(aA, sp(1), si(1), exp(b_nu(1)+b_nudim(1)*z), sv(1), cpv(1));
    vB = psyrat_gamma_varcomps(aB, sp(2), si(2), exp(b_nu(2)+b_nudim(2)*z), sv(2), cpv(2));
    cc = psyrat_gamma_crosscov(aA, sp(1), si(1), aB, sp(2), si(2), rho_p, rho_i);
    idv = zeros(1,2,2); trv = zeros(1,2,2);
    idv(1,1,1) = vA.sigma_p2; idv(1,2,2) = vB.sigma_p2;
    idv(1,1,2) = cc.cov_p_obs; idv(1,2,1) = cc.cov_p_obs;
    trv(1,1,1) = vA.sigma_i2; trv(1,2,2) = vB.sigma_i2;
    trv(1,1,2) = cc.cov_i_obs; trv(1,2,1) = cc.cov_i_obs;
    er = [0.5*log(vA.sigma_pi_e2), 0.5*log(vB.sigma_pi_e2)];
    g = psyrat_diffrel('bp',idv,'bt',trv,'er_var',er,'obs',obs,'CI',ci,...
        'est','gen','wp_cov',0);
    d = psyrat_diffrel('bp',idv,'bt',trv,'er_var',er,'obs',obs,'CI',ci,...
        'est','dep','wp_cov',0);
    trueG(k) = g.pt; trueD(k) = d.pt;
end
end

function rel = localRunDiffDynrelGamma(testCase, dataTbl, tag)
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
    'diffest', 2, ...
    'chains', 2, ...
    'warmup', 1000, ...
    'sampling', 1000, ...
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
        % best effort
    end
end
end
