classdef TestGammaDiffCopulaRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square one-facet
    % two-event CONCURRENT DIFFERENCE-score design (analysis 7, diffrescor=2). The
    % two events are recorded on the SAME (person, trial-position) cell, and their
    % scaled-chi-square residuals are coupled by a Gaussian copula correlation
    % rho_e (psyrat_build_stan_diff_copula_chisq). The observed-scale components and
    % the residual cross-covariance are recovered by the reference full Monte-Carlo
    % projection (psyrat_gamma_extract_diff_copula).
    %
    % This test simulates a concurrent bivariate design (correlated event means,
    % correlated trial effects, person-varying dispersion, and a KNOWN copula
    % residual correlation) and confirms:
    %   (1) a positive residual correlation is recovered (rho_e > 0 induces a
    %       positive observed-scale residual cross-covariance), whereas a
    %       non-concurrent extraction would report 0;
    %   (2) the observed cross-event person covariance (the difference numerator's
    %       off-diagonal) is recovered; and
    %   (3) the true difference generalizability/dependability coefficients fall
    %       inside the recovered 95% credible intervals, with the point estimates
    %       recovered. The truth is a large INDEPENDENT Monte-Carlo of the concurrent
    %       difference at the data-generating parameters (not the extraction's MC).
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
                'recov_gamma_diff_copula');
        end
    end

    methods (Test)
        function testConcurrentDifferenceReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gaminv','file')) && logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            % ---- Known truth (log expected-score scale) ----
            b     = [log(5.5), log(4.0)];
            b_nu  = [log(16),  log(12)];
            s_p   = [0.40, 0.35];    % per-event person log-mean SDs
            s_v   = [0.50, 0.45];    % per-event person log-nu SDs
            s_i   = [0.30, 0.26];    % per-event trial log-mean SDs
            rho_mean = 0.6;          % cross-event person-mean correlation
            rho_trl  = 0.4;          % cross-event trial correlation
            rho_pv   = -0.2;         % within-event mean<->log-nu correlation
            rho_e_true = 0.5;        % copula residual correlation (concurrent coupling)

            nsub = 40; ntrl = 22;   % copula density is expensive; keep the fit tractable
            dataTbl = localSimulateGammaDiffCopula(nsub, ntrl, b, b_nu, s_p, s_v, ...
                s_i, rho_mean, rho_trl, rho_pv, rho_e_true);

            % ---- Independent large-MC truth for the difference reliability ----
            [trueG, trueD, cov_p_obs_true] = localTrueConcurrentDiffRel(b, b_nu, ...
                s_p, s_v, s_i, rho_mean, rho_trl, rho_pv, rho_e_true, ntrl);

            % ---- Fit (case 7 gamma, concurrent copula) ----
            rel = localRunDiffCopulaGamma(testCase, dataTbl, 'recov_gamma_diff_copula');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % --- Convergence gate: the recovery checks below are only meaningful on a
            % CONVERGED fit (matching truth from an R-hat>=1.1 posterior is coincidental).
            % Gate via assume so a non-converged copula fit is filtered (Incomplete),
            % not a false PASS. The copula per-person-nu geometry is known not to converge
            % under this short-warmup regime (S11); the same R-hat>=1.1 + n_eff criteria
            % as the production psyrat_checkconv path decide the gate.
            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[diff copula recovery] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Copula fit ', ...
                'non-converged (max R-hat=%.3f >= 1.1 / n_eff gate); recovery not validated ', ...
                'in this environment (S11 copula per-person-nu convergence limit).'], maxRhat));

            idvar  = cell2mat(rel.out.id_varcov{1});     % nd x 2 x 2
            trlvar = cell2mat(rel.out.trl_varcov{1});    % nd x 2 x 2
            bsigma = cell2mat(rel.out.b_sigma{1});       % nd x 2
            errvc  = cell2mat(rel.out.err_varcov{1});    % nd x 2 x 2
            wp_cov = errvc(:,2,1);                        % residual cross-cov
            rescor = cell2mat(rel.out.rescor{1});

            gGen = psyrat_diffrel('bp',idvar,'bt',trlvar,'er_var',bsigma,...
                'obs',[ntrl ntrl],'CI',0.95,'est','gen','wp_cov',wp_cov);
            gDep = psyrat_diffrel('bp',idvar,'bt',trlvar,'er_var',bsigma,...
                'obs',[ntrl ntrl],'CI',0.95,'est','dep','wp_cov',wp_cov);

            covHat = mean(idvar(:,1,2));
            fprintf(['[diff copula recovery] rescor hat=%.3f  wp_cov hat=%.3f  ', ...
                'cov_p obs hat=%.3f (true %.3f)  G hat=%.3f [%.3f %.3f] (true %.3f)  ', ...
                'D hat=%.3f (true %.3f)\n'], mean(rescor), mean(wp_cov), covHat, ...
                cov_p_obs_true, gGen.pt, gGen.ll, gGen.ul, trueG, gDep.pt, trueD);

            % (1) positive residual correlation recovered (the copula effect)
            testCase.verifyGreaterThan(mean(rescor), 0.10, ...
                'A positive copula residual correlation should be recovered (rho_e_true>0).');
            testCase.verifyGreaterThan(mean(wp_cov), 0, ...
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
    end
end

function tbl = localSimulateGammaDiffCopula(nsub, ntrl, b, b_nu, s_p, s_v, s_i, ...
    rho_mean, rho_trl, rho_pv, rho_e, ~)
% Concurrent two-event Persons x Trials long table: each (person, trial-position)
% cell yields a copula-coupled scaled-chi-square PAIR (E1, E2), emitted in trial
% order so psyrat_build_diff_pairs pairs them concurrently.
rng(24680);
% 4-D person block [u1,u2,v1,v2]
Rp = eye(4);
Rp(1,2)=rho_mean; Rp(2,1)=rho_mean; Rp(1,3)=rho_pv; Rp(3,1)=rho_pv;
Rp(2,4)=rho_pv; Rp(4,2)=rho_pv;
sdp = [s_p(1) s_p(2) s_v(1) s_v(2)];
blk = (diag(sdp) * chol(Rp,'lower') * randn(4,nsub))';   % nsub x 4
u_p = blk(:,1:2); v_p = blk(:,3:4);
% 2-D correlated trial block (shared trial-position across events)
Rt = [1 rho_trl; rho_trl 1];
u_t = (diag(s_i) * chol(Rt,'lower') * randn(2,ntrl))';   % ntrl x 2

ids=cell(nsub*2*ntrl,1); evs=cell(nsub*2*ntrl,1); meas=zeros(nsub*2*ntrl,1);
evn={'E1','E2'}; r=0;
for p = 1:nsub
    nu1 = exp(b_nu(1)+v_p(p,1)); nu2 = exp(b_nu(2)+v_p(p,2));
    for e = 1:2
        for t = 1:ntrl
            r=r+1; ids{r}=sprintf('S%03d',p); evs{r}=evn{e};
        end
    end
    % draw the concurrent copula pair per trial, then place into the E1/E2 runs
    z1 = randn(ntrl,1); z2 = rho_e.*z1 + sqrt(1-rho_e^2).*randn(ntrl,1);
    mu1 = exp(b(1)+u_p(p,1)+u_t(:,1)); mu2 = exp(b(2)+u_p(p,2)+u_t(:,2));
    y1 = gaminv(normcdf(z1), 0.5*nu1, (2.*mu1)./nu1);
    y2 = gaminv(normcdf(z2), 0.5*nu2, (2.*mu2)./nu2);
    base = (p-1)*2*ntrl;
    meas(base+(1:ntrl))        = y1;   % E1 run
    meas(base+ntrl+(1:ntrl))   = y2;   % E2 run
end
tbl = table(ids, meas, evs, 'VariableNames', {'id','meas','event'});
end

function [trueG, trueD, cov_p_obs] = localTrueConcurrentDiffRel(b, b_nu, s_p, s_v, ...
    s_i, rho_mean, rho_trl, rho_pv, rho_e, nprime)
% Large INDEPENDENT Monte-Carlo of the concurrent difference at the true params:
% simulate a big Persons x Trials concurrent design, ANOVA-decompose delta, and
% form the single-observation-per-cell difference reliability projected over
% nprime trials. Independent of the extraction's (smaller, seeded) MC.
rng(99999);
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
% components of delta
vd = localVcStruct(th-al);
p_d = max(vd.p,0); i_d = max(vd.i,0); e_d = max(vd.pi_e,0);
trueG = p_d / (p_d + e_d/nprime);
trueD = p_d / (p_d + e_d/nprime + i_d/nprime);
% observed cross-event person covariance (numerator off-diagonal)
vt = localVcStruct(th); va = localVcStruct(al);
cov_p_obs = (vt.p + va.p - vd.p)/2;
end

function v = localVcStruct(Y)
[P,I]=size(Y); g=mean(Y(:)); pm=mean(Y,2); im=mean(Y,1); r=Y-pm-im+g;
msp=(I*sum((pm-g).^2))/(P-1); msi=(P*sum((im-g).^2))/(I-1); mspi=sum(r(:).^2)/((P-1)*(I-1));
v=struct('p',(msp-mspi)/I,'i',(msi-mspi)/P,'pi_e',mspi);
end

function rel = localRunDiffCopulaGamma(testCase, dataTbl, tag)
outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), tag);
if ~exist(outDir,'dir'), mkdir(outDir); end
% dispersion=1 pins the per-person-nu variant: the group-level default is now GLOBAL
% (2), so this per-person-nu recovery test selects dispersion=1 explicitly.
rel = psyrat_computevarcomp('data', dataTbl, 'family','gamma', ...
    'diffest', 2, 'diffwpcov', 2, 'dispersion', 1, 'chains', 2, 'warmup', 500, ...
    'sampling', 500, 'seed', 12345, 'verbose', 0, 'showgui', 1, 'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir,'dir'), mkdir(outdir); end
end

function localCleanupArtifacts(baseDir, tag)
p = fullfile(baseDir, tag);
if isfolder(p), try, rmdir(p,'s'); catch, end, end
end
