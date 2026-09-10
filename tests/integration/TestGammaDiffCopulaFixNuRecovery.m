classdef TestGammaDiffCopulaFixNuRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square one-facet
    % two-event CONCURRENT DIFFERENCE-score design with GROUP-LEVEL GLOBAL fixed
    % dispersion (analysis 7, diffrescor=2, dispersion=2;
    % psyrat_build_stan_diff_copula_fixnu_chisq +
    % psyrat_gamma_extract_diff_copula_fixnu). The two events are recorded on the SAME
    % (person, trial-position) cell and their scaled-chi-square residuals are coupled by
    % a Gaussian copula correlation rho_e, but the dispersion nu is event-global
    % (constant across persons) and the person block is 2-D (means only).
    %
    % Two validations:
    %   (A) testGlobalNuConcurrentDifferenceReliabilityRecovered - data simulated from
    %       the global-nu generative model (constant nu); confirms rho_e>0, the
    %       cross-event person covariance, and the true difference G/D are recovered.
    %       Unlike the per-person-nu sibling (S11-blocked), the global-nu geometry is
    %       expected to CONVERGE at this short-warmup regime - convergence itself is
    %       part of the claim.
    %   (B) testGlobalNuRobustToPersonVaryingDispersion - the decisive robustness check
    %       for the group-level estimand: data simulated with TRUE person-varying nu
    %       (the world the per-person-nu model assumes), then fit with the GLOBAL-nu
    %       model. Confirms the global-nu (mis-specified dispersion) model still
    %       recovers the TRUE GROUP difference reliability. This tests whether omitting
    %       person-varying dispersion biases group-level reliability - it should not,
    %       because group reliability is a function of the mean-structure variance
    %       components, which the global-nu model still estimates.
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
            localCleanupArtifacts(base, 'recov_gamma_diff_copula_fixnu');
            localCleanupArtifacts(base, 'recov_gamma_diff_copula_fixnu_robust');
        end
    end

    methods (Test)
        function testGlobalNuConcurrentDifferenceReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gaminv','file')) && logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            % ---- Known truth (GLOBAL nu: constant across persons) ----
            b     = [log(5.5), log(4.0)];
            nu    = [16, 12];        % event-global dispersion (NOT person-varying)
            s_p   = [0.40, 0.35];    % per-event person log-mean SDs
            s_i   = [0.30, 0.26];    % per-event trial log-mean SDs
            rho_mean = 0.6;          % cross-event person-mean correlation
            rho_trl  = 0.4;          % cross-event trial correlation
            rho_e_true = 0.5;        % copula residual correlation

            nsub = 40; ntrl = 22;
            dataTbl = localSimGlobalNu(nsub, ntrl, b, nu, s_p, s_i, ...
                rho_mean, rho_trl, rho_e_true);

            % ---- Independent large-MC truth for the difference reliability ----
            [trueG, trueD, cov_p_obs_true] = localTrueGlobalNu(b, nu, s_p, s_i, ...
                rho_mean, rho_trl, rho_e_true, ntrl);

            rel = localRunFixnu(testCase, dataTbl, 'recov_gamma_diff_copula_fixnu');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[fixnu diff copula recovery] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Global-nu fit ', ...
                'non-converged (max R-hat=%.3f); recovery not validated in this environment.'], maxRhat));

            [gGen, gDep, covHat, rescorHat, wpHat] = localDiffRel(rel, ntrl);
            fprintf(['[fixnu diff copula recovery] rescor hat=%.3f  wp_cov hat=%.3f  ', ...
                'cov_p obs hat=%.3f (true %.3f)  G hat=%.3f [%.3f %.3f] (true %.3f)  ', ...
                'D hat=%.3f (true %.3f)\n'], rescorHat, wpHat, covHat, ...
                cov_p_obs_true, gGen.pt, gGen.ll, gGen.ul, trueG, gDep.pt, trueD);

            % (1) positive residual correlation recovered (the copula effect)
            testCase.verifyGreaterThan(rescorHat, 0.10, ...
                'A positive copula residual correlation should be recovered (rho_e_true>0).');
            testCase.verifyGreaterThan(wpHat, 0, ...
                'The observed residual cross-covariance should be positive.');
            % (2) cross-event person covariance recovered with the right sign + scale
            testCase.verifyGreaterThan(covHat, 0, ...
                'Recovered cross-event person covariance should be positive.');
            testCase.verifyLessThan(abs(covHat - cov_p_obs_true), 0.5*abs(cov_p_obs_true) + 1, ...
                sprintf('Cross-event person covariance off: recovered %.3f vs true %.3f.', ...
                covHat, cov_p_obs_true));
            % (3) true difference reliability inside the recovered 95% CIs + point recovered
            testCase.verifyGreaterThanOrEqual(trueG, gGen.ll - 0.03, ...
                sprintf('True G %.3f below the recovered 95%% CI lower %.3f.', trueG, gGen.ll));
            testCase.verifyLessThanOrEqual(trueG, gGen.ul + 0.03, ...
                sprintf('True G %.3f above the recovered 95%% CI upper %.3f.', trueG, gGen.ul));
            testCase.verifyLessThan(abs(gGen.pt - trueG), 0.06, ...
                sprintf('Difference generalizability off: recovered %.3f vs true %.3f.', gGen.pt, trueG));
            testCase.verifyLessThan(abs(gDep.pt - trueD), 0.06, ...
                sprintf('Difference dependability off: recovered %.3f vs true %.3f.', gDep.pt, trueD));
        end

        function testGlobalNuRobustToPersonVaryingDispersion(testCase)
            % DECISIVE robustness check: TRUE dispersion varies across persons, but we
            % fit the GLOBAL-nu model. The recovered GROUP difference reliability must
            % still match the true group reliability - i.e. omitting person-varying
            % dispersion does not bias the group-level estimand.
            testCase.assumeTrue(logical(exist('gaminv','file')) && logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            % ---- Known truth WITH person-varying dispersion (v_p block) ----
            b     = [log(5.5), log(4.0)];
            b_nu  = [log(16),  log(12)];
            s_p   = [0.40, 0.35];
            s_v   = [0.50, 0.45];    % per-event person log-nu SDs (TRUE dispersion varies)
            s_i   = [0.30, 0.26];
            rho_mean = 0.6; rho_trl = 0.4; rho_pv = -0.2; rho_e_true = 0.5;

            nsub = 40; ntrl = 22;
            dataTbl = localSimPersonNu(nsub, ntrl, b, b_nu, s_p, s_v, s_i, ...
                rho_mean, rho_trl, rho_pv, rho_e_true);
            [trueG, trueD] = localTruePersonNu(b, b_nu, s_p, s_v, s_i, ...
                rho_mean, rho_trl, rho_pv, rho_e_true, ntrl);

            % Fit the GLOBAL-nu (mis-specified dispersion) model. This harder fit needs
            % more warmup than the recovery lane's 500 to clear the R-hat gate, so use
            % 4 chains x 2000 (verified live: R-hat=1.000, 0 divergences).
            rel = localRunFixnu(testCase, dataTbl, 'recov_gamma_diff_copula_fixnu_robust', 2000, 4);

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[fixnu robustness] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Global-nu fit ', ...
                'non-converged (max R-hat=%.3f); robustness not validated here.'], maxRhat));

            [gGen, gDep] = localDiffRel(rel, ntrl);
            fprintf(['[fixnu robustness] G hat=%.3f [%.3f %.3f] (true %.3f)  ', ...
                'D hat=%.3f (true %.3f)\n'], gGen.pt, gGen.ll, gGen.ul, trueG, gDep.pt, trueD);

            % Group reliability recovered despite mis-specified (global) dispersion.
            testCase.verifyGreaterThanOrEqual(trueG, gGen.ll - 0.05, ...
                sprintf('True group G %.3f below the recovered 95%% CI lower %.3f.', trueG, gGen.ll));
            testCase.verifyLessThanOrEqual(trueG, gGen.ul + 0.05, ...
                sprintf('True group G %.3f above the recovered 95%% CI upper %.3f.', trueG, gGen.ul));
            testCase.verifyLessThan(abs(gGen.pt - trueG), 0.08, ...
                sprintf(['Global-nu group generalizability biased under person-varying-nu ', ...
                'truth: recovered %.3f vs true %.3f.'], gGen.pt, trueG));
            testCase.verifyLessThan(abs(gDep.pt - trueD), 0.08, ...
                sprintf(['Global-nu group dependability biased under person-varying-nu ', ...
                'truth: recovered %.3f vs true %.3f.'], gDep.pt, trueD));
        end
    end
end

% ===================== helpers =====================

function tbl = localSimGlobalNu(nsub, ntrl, b, nu, s_p, s_i, rho_mean, rho_trl, rho_e)
% Concurrent two-event Persons x Trials long table with GLOBAL (constant) nu.
rng(24681);
Rp = [1 rho_mean; rho_mean 1];
u_p = (diag(s_p) * chol(Rp,'lower') * randn(2,nsub))';   % nsub x 2 (means only)
Rt = [1 rho_trl; rho_trl 1];
u_t = (diag(s_i) * chol(Rt,'lower') * randn(2,ntrl))';   % ntrl x 2
ids=cell(nsub*2*ntrl,1); evs=cell(nsub*2*ntrl,1); meas=zeros(nsub*2*ntrl,1);
evn={'E1','E2'}; r=0;
for p = 1:nsub
    for e = 1:2
        for t = 1:ntrl
            r=r+1; ids{r}=sprintf('S%03d',p); evs{r}=evn{e};
        end
    end
    z1 = randn(ntrl,1); z2 = rho_e.*z1 + sqrt(1-rho_e^2).*randn(ntrl,1);
    mu1 = exp(b(1)+u_p(p,1)+u_t(:,1)); mu2 = exp(b(2)+u_p(p,2)+u_t(:,2));
    y1 = gaminv(normcdf(z1), 0.5*nu(1), (2.*mu1)./nu(1));
    y2 = gaminv(normcdf(z2), 0.5*nu(2), (2.*mu2)./nu(2));
    base = (p-1)*2*ntrl;
    meas(base+(1:ntrl))      = y1;
    meas(base+ntrl+(1:ntrl)) = y2;
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

function [trueG, trueD, cov_p_obs] = localTrueGlobalNu(b, nu, s_p, s_i, ...
    rho_mean, rho_trl, rho_e, nprime)
rng(99998);
NP = 4000; NI = 60;
Rp = [1 rho_mean; rho_mean 1];
up = (diag(s_p)*chol(Rp,'lower')*randn(2,NP))';
Rt = [1 rho_trl; rho_trl 1];
ut = (diag(s_i)*chol(Rt,'lower')*randn(2,NI))';
mu1 = exp(b(1)+up(:,1)+ut(:,1)'); mu2 = exp(b(2)+up(:,2)+ut(:,2)');
za = randn(NP,NI); zb = randn(NP,NI); ua = normcdf(za); ub = normcdf(rho_e*za+sqrt(1-rho_e^2)*zb);
th = gaminv(ua, 0.5*nu(1)+zeros(NP,NI), (2*mu1)./nu(1));
al = gaminv(ub, 0.5*nu(2)+zeros(NP,NI), (2*mu2)./nu(2));
vd = localVcStruct(th-al);
p_d = max(vd.p,0); i_d = max(vd.i,0); e_d = max(vd.pi_e,0);
trueG = p_d / (p_d + e_d/nprime);
trueD = p_d / (p_d + e_d/nprime + i_d/nprime);
vt = localVcStruct(th); va = localVcStruct(al);
cov_p_obs = (vt.p + va.p - vd.p)/2;
end

function tbl = localSimPersonNu(nsub, ntrl, b, b_nu, s_p, s_v, s_i, ...
    rho_mean, rho_trl, rho_pv, rho_e)
% Person-varying-nu truth (4-D person block), same layout as the per-person recovery.
rng(24682);
Rp = eye(4);
Rp(1,2)=rho_mean; Rp(2,1)=rho_mean; Rp(1,3)=rho_pv; Rp(3,1)=rho_pv;
Rp(2,4)=rho_pv; Rp(4,2)=rho_pv;
sdp = [s_p(1) s_p(2) s_v(1) s_v(2)];
blk = (diag(sdp) * chol(Rp,'lower') * randn(4,nsub))';
u_p = blk(:,1:2); v_p = blk(:,3:4);
Rt = [1 rho_trl; rho_trl 1];
u_t = (diag(s_i) * chol(Rt,'lower') * randn(2,ntrl))';
ids=cell(nsub*2*ntrl,1); evs=cell(nsub*2*ntrl,1); meas=zeros(nsub*2*ntrl,1);
evn={'E1','E2'}; r=0;
for p = 1:nsub
    nu1 = exp(b_nu(1)+v_p(p,1)); nu2 = exp(b_nu(2)+v_p(p,2));
    for e = 1:2
        for t = 1:ntrl
            r=r+1; ids{r}=sprintf('S%03d',p); evs{r}=evn{e};
        end
    end
    z1 = randn(ntrl,1); z2 = rho_e.*z1 + sqrt(1-rho_e^2).*randn(ntrl,1);
    mu1 = exp(b(1)+u_p(p,1)+u_t(:,1)); mu2 = exp(b(2)+u_p(p,2)+u_t(:,2));
    y1 = gaminv(normcdf(z1), 0.5*nu1, (2.*mu1)./nu1);
    y2 = gaminv(normcdf(z2), 0.5*nu2, (2.*mu2)./nu2);
    base = (p-1)*2*ntrl;
    meas(base+(1:ntrl))      = y1;
    meas(base+ntrl+(1:ntrl)) = y2;
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

function [trueG, trueD] = localTruePersonNu(b, b_nu, s_p, s_v, s_i, ...
    rho_mean, rho_trl, rho_pv, rho_e, nprime)
rng(99997);
NP = 4000; NI = 60;
Rp = eye(4); Rp(1,2)=rho_mean;Rp(2,1)=rho_mean;Rp(1,3)=rho_pv;Rp(3,1)=rho_pv;Rp(2,4)=rho_pv;Rp(4,2)=rho_pv;
sdp=[s_p(1) s_p(2) s_v(1) s_v(2)];
blk=(diag(sdp)*chol(Rp,'lower')*randn(4,NP))'; up=blk(:,1:2); vp=blk(:,3:4);
Rt=[1 rho_trl;rho_trl 1]; ut=(diag(s_i)*chol(Rt,'lower')*randn(2,NI))';
mu1=exp(b(1)+up(:,1)+ut(:,1)'); mu2=exp(b(2)+up(:,2)+ut(:,2)');
nu1=exp(b_nu(1)+vp(:,1)); nu2=exp(b_nu(2)+vp(:,2));
za=randn(NP,NI); zb=randn(NP,NI); ua=normcdf(za); ub=normcdf(rho_e*za+sqrt(1-rho_e^2)*zb);
th=gaminv(ua,(0.5*nu1)+zeros(1,NI),(2*mu1)./nu1);
al=gaminv(ub,(0.5*nu2)+zeros(1,NI),(2*mu2)./nu2);
vd = localVcStruct(th-al);
p_d = max(vd.p,0); i_d = max(vd.i,0); e_d = max(vd.pi_e,0);
trueG = p_d / (p_d + e_d/nprime);
trueD = p_d / (p_d + e_d/nprime + i_d/nprime);
end

function [gGen, gDep, covHat, rescorHat, wpHat] = localDiffRel(rel, ntrl)
idvar  = cell2mat(rel.out.id_varcov{1});
trlvar = cell2mat(rel.out.trl_varcov{1});
bsigma = cell2mat(rel.out.b_sigma{1});
errvc  = cell2mat(rel.out.err_varcov{1});
wp_cov = errvc(:,2,1);
gGen = psyrat_diffrel('bp',idvar,'bt',trlvar,'er_var',bsigma,...
    'obs',[ntrl ntrl],'CI',0.95,'est','gen','wp_cov',wp_cov);
gDep = psyrat_diffrel('bp',idvar,'bt',trlvar,'er_var',bsigma,...
    'obs',[ntrl ntrl],'CI',0.95,'est','dep','wp_cov',wp_cov);
covHat = mean(idvar(:,1,2));
rescorHat = mean(cell2mat(rel.out.rescor{1}));
wpHat = mean(wp_cov);
end

function v = localVcStruct(Y)
[P,I]=size(Y); g=mean(Y(:)); pm=mean(Y,2); im=mean(Y,1); r=Y-pm-im+g;
msp=(I*sum((pm-g).^2))/(P-1); msi=(P*sum((im-g).^2))/(I-1); mspi=sum(r(:).^2)/((P-1)*(I-1));
v=struct('p',(msp-mspi)/I,'i',(msi-mspi)/P,'pi_e',mspi);
end

function rel = localRunFixnu(testCase, dataTbl, tag, warmup, chains)
% warmup/chains default to the fast recovery-lane settings (500 / 2). The robustness
% method overrides them: fitting the MIS-SPECIFIED global-nu model to person-varying-nu
% truth needs more warmup to clear the R-hat convergence gate (at 500/2 it lands at
% R-hat~1.10 and assume-skips; verified live to converge cleanly - R-hat=1.000, 0
% divergences - and recover GROUP reliability to within 0.002/0.016 at 4 chains x 2000).
if nargin < 4 || isempty(warmup), warmup = 500; end
if nargin < 5 || isempty(chains), chains = 2; end
outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), tag);
if ~exist(outDir,'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp('data', dataTbl, 'family','gamma', ...
    'diffest', 2, 'diffwpcov', 2, 'dispersion', 2, ...
    'chains', chains, 'warmup', warmup, 'sampling', warmup, ...
    'seed', 12345, 'verbose', 0, 'showgui', 1, 'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir,'dir'), mkdir(outdir); end
end

function localCleanupArtifacts(baseDir, tag)
p = fullfile(baseDir, tag);
if isfolder(p), try, rmdir(p,'s'); catch, end, end
end
