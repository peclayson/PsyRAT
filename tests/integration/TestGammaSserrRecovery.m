classdef TestGammaSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL design (sserr, analysis 6). The gamma subject-level model IS
    % the one-facet gamma model: the person-varying log-dispersion log nu_p =
    % beta + v_p already gives each participant their own scaled-chi-square
    % residual, so per-subject reliability comes from the CONDITIONAL residual
    %   sigma_pi_e2(p) = 2*exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p))
    % (psyrat_gamma_extract_sserr), fed through the Gaussian subject-level path
    % (psyrat_ssrel). This test simulates a person-VARYING dispersion (so the
    % per-subject residuals genuinely differ) and confirms:
    %   (1) the fit recovers the person log-nu SD s_v (the mechanism that makes
    %       the design subject-level);
    %   (2) the recovered per-subject residual DIFFERENTIATES participants (it
    %       tracks the true per-subject residual across the sample); and
    %   (3) the true per-subject generalizability coefficients fall inside the
    %       recovered 95% credible intervals (coverage), with the sample-mean
    %       coefficient recovered.
    %
    % Ground truth is the closed form validated offline (no CmdStan) in
    % TestGammaSserrConversion. CmdStan-gated (skips Incomplete, not Failed, when
    % CmdStan/MatlabStan absent).

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
        function testSubjectLevelReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   % log-mean grand intercept
            s_p   = 0.40;       % person log-mean SD
            s_i   = 0.20;       % trial  log-mean SD
            beta  = log(15);    % population log-nu intercept
            s_v   = 0.60;       % person log-nu SD (dispersion heterogeneity)
            rho   = -0.30;      % person mean <-> log-nu correlation
            nu_pop = exp(beta);
            cov_pv = rho * s_p * s_v;

            nsub = 60; ntrl = 40;
            [dataTbl, u_p, v_p] = localSimulateGammaSserr( ...
                nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho);

            % ---- Closed-form observed-scale truth ----
            % Group person universe variance (numerator; estimand-invariant) and
            % each participant's own CONDITIONAL residual (observed scale).
            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, nu_pop, s_v, cov_pv);
            sigma_p2_true = vc.sigma_p2;
            trueResid = 2 .* exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p));  % nsub x 1
            trueG = sigma_p2_true ./ (sigma_p2_true + trueResid ./ ntrl);    % generalizability

            % ---- Fit ----
            rel = localRunSserrGamma(testCase, dataTbl, 'recov_gamma_sserr');

            testCase.verifyEqual(rel.analysis, 'ic_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
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

            % align the recovered per-subject quantities to the truth by parsing
            % the original id string ('S###' -> subject index)
            nrows = height(ssrel);
            recResid = zeros(nrows,1); tG = zeros(nrows,1);
            recG = ssrel.dep_pt; recLL = ssrel.dep_ll; recUL = ssrel.dep_ul;
            for r = 1:nrows
                idx = sscanf(char(string(ssrel.id(r))), 'S%d');
                recResid(r) = ssrel.ss_errvar(r);
                tG(r) = trueG(idx);
            end
            trueResidAligned = arrayfun(@(r) trueResid(sscanf(char(string(ssrel.id(r))),'S%d')), (1:nrows)');

            % (1) s_v recovered (person log-nu SD is gro_sds column 2, provenance)
            svHat = mean(cell2mat(rel.out.gro_sds{1}(:,2)));
            fprintf(['[sserr recovery] s_v hat=%.3f (true %.2f)  corr(resid)=%.3f  ', ...
                'coverage=%.2f  mean recG=%.3f (true %.3f)\n'], ...
                svHat, s_v, localPearson(recResid, trueResidAligned), ...
                mean(tG >= recLL & tG <= recUL), mean(recG), mean(tG));
            testCase.verifyGreaterThan(svHat, 0.40, sprintf( ...
                'Recovered person log-nu SD %.3f far below the true %.2f.', svHat, s_v));
            testCase.verifyLessThan(svHat, 0.85);

            % (2) subject-level DIFFERENTIATION: recovered per-subject residual
            % tracks the true per-subject residual across participants.
            rResid = localPearson(recResid, trueResidAligned);
            testCase.verifyGreaterThan(rResid, 0.5, sprintf( ...
                ['Recovered per-subject residual does not track the truth ', ...
                '(Pearson r=%.3f); subject-level differentiation not recovered.'], rResid));

            % (3) COVERAGE: most true per-subject coefficients fall inside their
            % recovered 95%% CIs, and the sample-mean coefficient is recovered.
            coverage = mean(tG >= recLL & tG <= recUL);
            testCase.verifyGreaterThanOrEqual(coverage, 0.80, sprintf( ...
                'Only %.0f%% of true per-subject coefficients inside the 95%% CIs.', ...
                100*coverage));
            testCase.verifyLessThan(abs(mean(recG) - mean(tG)), 0.05, sprintf( ...
                'Sample-mean per-subject reliability off: recovered %.3f vs true %.3f.', ...
                mean(recG), mean(tG)));
        end
    end
end

function [tbl, u_p, v_p] = localSimulateGammaSserr(nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho)
% One-facet Persons x Trials long table from a log-linked mean (person + crossed
% trial main effects) and a person-VARYING scaled-chi-square dispersion. Returns
% the per-subject (u_p, v_p) so the caller can form the CONDITIONAL per-subject
% residual truth. Y = gamrnd(nu_p/2, 2*mu/nu_p) has E[Y]=mu, Var[Y]=2*mu^2/nu_p.
rng(12345);
u_t = s_i * randn(ntrl, 1);                          % crossed trial main effects
z1 = randn(nsub, 1); z2 = randn(nsub, 1);
u_p = s_p * z1;                                      % person mean effect
v_p = s_v * (rho * z1 + sqrt(1 - rho^2) * z2);       % person log-nu effect (corr rho)
nu_p = exp(beta + v_p);

nrow = nsub * ntrl;
ids  = cell(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    for t = 1:ntrl
        r = r + 1;
        ids{r}  = sprintf('S%03d', p);
        mu = exp(alpha + u_p(p) + u_t(t));
        meas(r) = gamrnd(nu_p(p) / 2, 2 * mu / nu_p(p));
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function r = localPearson(x, y)
% Pearson correlation via base MATLAB (no Statistics Toolbox).
x = x(:) - mean(x(:)); y = y(:) - mean(y(:));
r = (x' * y) / (sqrt(x' * x) * sqrt(y' * y));
end

function rel = localRunSserrGamma(testCase, dataTbl, tag)
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
p = fullfile(baseDir, 'recov_gamma_sserr');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            'Could not remove recovery artifact directory ''%s''.', p);
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
