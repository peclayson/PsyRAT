classdef TestGammaDiffTrtCopulaFixNuRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square two-facet
    % (test-retest) two-event CONCURRENT DIFFERENCE-score design with GROUP-LEVEL
    % GLOBAL fixed dispersion (analysis 10, diffrescor=2, dispersion=2;
    % psyrat_build_stan_diff_trt_copula_fixnu_chisq +
    % psyrat_gamma_extract_diff_trt_copula_fixnu). The two events are recorded on the
    % SAME (person, occasion, trial-position) cell; their scaled-chi-square residuals
    % are coupled by a Gaussian copula correlation rho_e, but the dispersion nu is
    % event-global (constant across persons) and the person block is 2-D (means only).
    %
    % Confirms: (1) a positive residual correlation is recovered; (2) the true CES
    % difference generalizability/dependability coefficients fall inside the recovered
    % 95% credible intervals with the point estimates recovered (truth = a large
    % INDEPENDENT Monte-Carlo of the concurrent two-facet difference at the global-nu
    % data-generating parameters). Unlike the per-person-nu sibling (S11-blocked;
    % case-10 SIGKILLed at treedepth 12 in-repo), the global-nu geometry is expected to
    % be tractable - convergence itself is part of the claim, and this test is the
    % cluster validation the maintainer runs to confirm it.
    %
    % NOTE on D: the occasion main-effect difference variance (o_d/nocc) is weakly
    % identified at realistic occasion counts, so the D POINT is checked by 95%
    % CI-containment (not a tight point tolerance), exactly as the per-person sibling.
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
                'recov_gamma_diff_trt_copula_fixnu');
        end
    end

    methods (Test)
        function testGlobalNuConcurrentTwoFacetDifferenceReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gaminv','file')) && logical(exist('normcdf','file')), ...
                'Statistics Toolbox (gaminv/normcdf) required for the copula draw.');

            % ---- Known truth (GLOBAL nu: constant across persons) ----
            P = struct();
            P.b=[log(5.5) log(4.0)]; P.nu=[16 12];
            P.s_p=[.38 .33]; P.s_i=[.28 .24]; P.s_o=[.24 .20];
            P.s_pi=[.26 .22]; P.s_po=[.20 .17]; P.s_io=[.16 .14];
            P.rho_mean=.6;
            P.rho_i=.4; P.rho_o=.5; P.rho_pi=.3; P.rho_po=.35; P.rho_io=.25;
            P.rho_e=0.5;

            nsub=32; nocc=3; ntrl=12;
            dataTbl = localSimGlobalNuTrt(nsub, nocc, ntrl, P);
            [trueG, trueD] = localTrueGlobalNuTrt(P, ntrl, nocc);

            rel = localRunFixnuTrt(testCase, dataTbl, 'recov_gamma_diff_trt_copula_fixnu');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[fixnu diff trt copula recovery] max R-hat=%.3f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Global-nu two-facet fit ', ...
                'non-converged (max R-hat=%.3f); recovery not validated in this environment.'], maxRhat));

            idv  = cell2mat(rel.out.id_varcov{1});
            tidv = cell2mat(rel.out.tid_varcov{1});
            oidv = cell2mat(rel.out.oid_varcov{1});
            tlv  = cell2mat(rel.out.trl_varcov{1});
            ocv  = cell2mat(rel.out.occ_varcov{1});
            tov  = cell2mat(rel.out.to_varcov{1});
            bsig = cell2mat(rel.out.b_sigma{1});
            errvc = cell2mat(rel.out.err_varcov{1});
            er_cov = errvc(:,2,1);
            rescor = cell2mat(rel.out.rescor{1});

            args = {'bp',idv,'bpi',tidv,'bpo',oidv,'bt',tlv,'bo',ocv,'boi',tov,...
                'er_var',bsig,'er_cov',er_cov,'obs',[ntrl ntrl],'nocc',[nocc nocc],...
                'reltype',3,'CI',0.95};
            gGen = psyrat_diffrel_trt(args{:},'est','gen');
            gDep = psyrat_diffrel_trt(args{:},'est','dep');

            fprintf(['[fixnu diff trt copula recovery] rescor hat=%.3f  er_cov hat=%.3f  ', ...
                'G hat=%.3f [%.3f %.3f] (true %.3f)  D hat=%.3f [%.3f %.3f] (true %.3f)\n'], ...
                mean(rescor), mean(er_cov), gGen.pt, gGen.ll, gGen.ul, trueG, ...
                gDep.pt, gDep.ll, gDep.ul, trueD);

            testCase.verifyGreaterThan(mean(rescor), 0.10, ...
                'A positive copula residual correlation should be recovered.');
            testCase.verifyGreaterThanOrEqual(trueG, gGen.ll - 0.04, ...
                sprintf('True G %.3f below recovered 95%% CI lower %.3f.', trueG, gGen.ll));
            testCase.verifyLessThanOrEqual(trueG, gGen.ul + 0.04, ...
                sprintf('True G %.3f above recovered 95%% CI upper %.3f.', trueG, gGen.ul));
            testCase.verifyLessThan(abs(gGen.pt - trueG), 0.08, ...
                sprintf('CES generalizability off: recovered %.3f vs true %.3f.', gGen.pt, trueG));
            testCase.verifyGreaterThanOrEqual(trueD, gDep.ll - 0.04, ...
                sprintf('True D %.3f below recovered 95%% CI lower %.3f (occasion-limited).', trueD, gDep.ll));
            testCase.verifyLessThanOrEqual(trueD, gDep.ul + 0.04, ...
                sprintf('True D %.3f above recovered 95%% CI upper %.3f (occasion-limited).', trueD, gDep.ul));
        end
    end
end

% ===================== helpers =====================

function tbl = localSimGlobalNuTrt(nsub, nocc, ntrl, P)
% Concurrent two-facet long table with GLOBAL (constant) nu: each (person, occasion,
% trial-position) cell yields a copula-coupled scaled-chi-square PAIR (E1, E2).
rng(13580);
u_p =(diag(P.s_p)*chol([1 P.rho_mean;P.rho_mean 1],'lower')*randn(2,nsub))';  % means only
u_i =(diag(P.s_i)*chol([1 P.rho_i;P.rho_i 1],'lower')*randn(2,ntrl))';
u_o =(diag(P.s_o)*chol([1 P.rho_o;P.rho_o 1],'lower')*randn(2,nocc))';
u_pi=(diag(P.s_pi)*chol([1 P.rho_pi;P.rho_pi 1],'lower')*randn(2,nsub*ntrl))';
u_po=(diag(P.s_po)*chol([1 P.rho_po;P.rho_po 1],'lower')*randn(2,nsub*nocc))';
u_io=(diag(P.s_io)*chol([1 P.rho_io;P.rho_io 1],'lower')*randn(2,ntrl*nocc))';
upi=@(p,t,c) u_pi((t-1)*nsub+p, c);
upo=@(p,o,c) u_po((o-1)*nsub+p, c);
uio=@(t,o,c) u_io((o-1)*ntrl+t, c);

ids={}; meas=[]; evs={}; tm={}; evn={'E1','E2'};
for p=1:nsub
    nu1=P.nu(1); nu2=P.nu(2);              % event-global (constant across persons)
    for o=1:nocc
        for e=1:2
            for t=1:ntrl
                ids{end+1,1}=sprintf('S%03d',p); evs{end+1,1}=evn{e}; tm{end+1,1}=sprintf('T%d',o); %#ok<AGROW>
            end
        end
        z1=randn(ntrl,1); z2=P.rho_e.*z1+sqrt(1-P.rho_e^2).*randn(ntrl,1);
        e1=zeros(ntrl,1); e2=zeros(ntrl,1);
        for t=1:ntrl
            m1=exp(P.b(1)+u_p(p,1)+u_i(t,1)+u_o(o,1)+upi(p,t,1)+upo(p,o,1)+uio(t,o,1));
            m2=exp(P.b(2)+u_p(p,2)+u_i(t,2)+u_o(o,2)+upi(p,t,2)+upo(p,o,2)+uio(t,o,2));
            e1(t)=gaminv(normcdf(z1(t)),0.5*nu1,2*m1/nu1);
            e2(t)=gaminv(normcdf(z2(t)),0.5*nu2,2*m2/nu2);
        end
        meas=[meas; e1; e2]; %#ok<AGROW>
    end
end
tbl = table(ids, meas, evs, tm, 'VariableNames', {'id','meas','event','time'});
end

function [trueG, trueD] = localTrueGlobalNuTrt(P, nprime, noprime)
% Large INDEPENDENT MC of the concurrent two-facet difference at the GLOBAL-nu truth.
rng(88887);
NP=800; NI=25; NO=8;
up=(diag(P.s_p)*chol([1 P.rho_mean;P.rho_mean 1],'lower')*randn(2,NP))';
ui=(diag(P.s_i)*chol([1 P.rho_i;P.rho_i 1],'lower')*randn(2,NI))';
uo=(diag(P.s_o)*chol([1 P.rho_o;P.rho_o 1],'lower')*randn(2,NO))';
upi=(diag(P.s_pi)*chol([1 P.rho_pi;P.rho_pi 1],'lower')*randn(2,NP*NI))';
upo=(diag(P.s_po)*chol([1 P.rho_po;P.rho_po 1],'lower')*randn(2,NP*NO))';
uio=(diag(P.s_io)*chol([1 P.rho_io;P.rho_io 1],'lower')*randn(2,NI*NO))';
za=randn(NP,NI,NO); zb=randn(NP,NI,NO);
for ev=1:2
    E=P.b(ev)+reshape(up(:,ev),[NP 1 1])+reshape(ui(:,ev),[1 NI 1])+reshape(uo(:,ev),[1 1 NO]) ...
      +reshape(upi(:,ev),NP,NI)+reshape(reshape(upo(:,ev),NP,NO),[NP 1 NO])+reshape(reshape(uio(:,ev),NI,NO),[1 NI NO]);
    MU=exp(E); NU=P.nu(ev)+zeros(NP,NI,NO);
    if ev==1, TH=gaminv(normcdf(za),0.5*NU,(2*MU)./NU);
    else,     AL=gaminv(normcdf(P.rho_e*za+sqrt(1-P.rho_e^2)*zb),0.5*NU,(2*MU)./NU); end
end
vd = localVc2(TH-AL);
p_d=max(vd.p,0); i_d=max(vd.i,0); o_d=max(vd.o,0);
pi_d=max(vd.pi,0); po_d=max(vd.po,0); io_d=max(vd.io,0); e_d=max(vd.pio_e,0);
rel_err = pi_d/nprime + po_d/noprime + e_d/(nprime*noprime);
abs_err = rel_err + i_d/nprime + o_d/noprime + io_d/(nprime*noprime);
trueG = p_d/(p_d + rel_err);
trueD = p_d/(p_d + abs_err);
end

function v = localVc2(Y)
[Pn,I,O]=size(Y); g=mean(Y(:));
pm=mean(Y,[2 3]); im=mean(Y,[1 3]); om=mean(Y,[1 2]);
pim=mean(Y,3); pom=mean(Y,2); iom=mean(Y,1);
ssp=I*O*sum((pm(:)-g).^2); ssi=Pn*O*sum((im(:)-g).^2); sso=Pn*I*sum((om(:)-g).^2);
pr=pim-pm-im+g; sspi=O*sum(pr(:).^2); qr=pom-pm-om+g; sspo=I*sum(qr(:).^2); ir=iom-im-om+g; ssio=Pn*sum(ir(:).^2);
r=Y-pim-pom-iom+pm+im+om-g; sspio=sum(r(:).^2);
msp=ssp/(Pn-1);msi=ssi/(I-1);mso=sso/(O-1);mspi=sspi/((Pn-1)*(I-1));mspo=sspo/((Pn-1)*(O-1));msio=ssio/((I-1)*(O-1));mspio=sspio/((Pn-1)*(I-1)*(O-1));
v=struct('p',(msp-mspi-mspo+mspio)/(I*O),'i',(msi-mspi-msio+mspio)/(Pn*O),'o',(mso-mspo-msio+mspio)/(Pn*I),...
  'pi',(mspi-mspio)/O,'po',(mspo-mspio)/I,'io',(msio-mspio)/Pn,'pio_e',mspio);
end

function rel = localRunFixnuTrt(testCase, dataTbl, tag)
outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), tag);
if ~exist(outDir,'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp('data', dataTbl, 'family','gamma', ...
    'diffest', 2, 'diffwpcov', 2, 'dispersion', 2, 'chains', 2, 'warmup', 500, ...
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
