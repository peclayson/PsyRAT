function out = psyrat_rel_dynrel_trt_gamma(varargin)
%Dynamic/conditional trial + occasion TWO-facet reliability for the scaled-chi-
%square / Gamma(log link) family, as a function of the standardized dimension
%predictor(s) (gamma analog of psyrat_rel_dynrel_trt; analysis 14, family =
%'gamma').
%
% out = psyrat_rel_dynrel_trt_gamma('alpha0',alpha0,'b',b,'lognu0',lognu0,...
%   'b_nu',b_nu,'sig_p',sig_p,'sig_v',sig_v,'rho_pv',rho_pv,...
%   'sig_occ',sig_occ,'sig_trl',sig_trl,'sig_trlxid',sig_trlxid,...
%   'sig_occxid',sig_occxid,'sig_trlxocc',sig_trlxocc,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',ntrl,'nocc',nocc,...
%   'reltype',reltype,'CI',.95)
%
%Unlike the Gaussian dynrel_trt (psyrat_rel_dynrel_trt), the gamma family has no
%log-residual scale submodel: the observed-score residual is the MEAN-COUPLED
%scaled-chi-square dispersion 2*mu^2/nu, so the log-mean and log-nu submodels are
%both dimension-conditioned,
%
%   log mu(z) = alpha0 + b . X(z)          log nu(z) = lognu0 + b_nu . X(z)
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two. At each z the log-scale crossed-design parameters are mapped to
%OBSERVED-score-scale variance components with the same lognormal-moment
%converter used by the static test-retest gamma path (psyrat_gamma_varcomps_trt),
%and the two-facet coefficient is the validated test-retest psyrat_rel_trt
%evaluated on those components. psyrat_rel_trt's three reltypes give:
%
%   reltype 1 (coefficient of equivalence / internal consistency: occasion fixed)
%   reltype 2 (coefficient of stability: trial fixed)
%   reltype 3 (coefficient of trial equivalence and stability: both random)
%
%Generalizability (gcoeff 2) excludes, and dependability (gcoeff 1) includes, the
%facet main effects, exactly as the static test-retest gamma path. NOTE that ALL
%components vary with z here (everything is on the multiplicative/log scale).
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's dispersion heterogeneity (estimand #2), the same convention the
%static gamma family uses: the person log-nu SD (sig_v) and its correlation with
%the person mean (rho_pv) are propagated into the converter as (s_v, cov_pv).
%There is no separate typical-person / population-average toggle.
%
%Inputs (name-value; all posterior-draw inputs are column vectors of length
%ndraws, or draws x KDIM for the slope blocks)
% alpha0     - log expected-score grand intercept draws (Intercept).
% b          - log-mean dimension-slope draws (draws x KDIM).
% lognu0     - log degrees-of-freedom grand intercept draws (Intercept_nu).
% b_nu       - log-nu dimension-slope draws (draws x KDIM).
% sig_p      - person log-mean SD draws (gro_sds(:,1) = s_p).
% sig_v      - person log-nu SD draws (gro_sds(:,2) = s_v).
% rho_pv     - person mean <-> log-nu correlation draws (chol_corrmat(:,2,1)).
% sig_occ    - occasion main-effect log-mean SD draws.
% sig_trl    - trial main-effect log-mean SD draws (s_i).
% sig_trlxid - trial x person log-mean SD draws (s_pt).
% sig_occxid - occasion x person log-mean SD draws (s_po).
% sig_trlxocc- trial x occasion log-mean SD draws (s_ot).
% ndim       - number of dimensions: 1 or 2.
% z1         - vector of standardized values for dimension 1 (the grid axis).
% obs        - n' (trial count) at which to evaluate the coefficients (scalar).
% CI         - credible-interval width in (0,1], e.g., .95.
%
%Optional Inputs
% z2       - vector of standardized values for dimension 2 (required if ndim=2).
% nocc     - n' (occasion count); default 1.
% reltype  - 1, 2, or 3 (see above); default 3 (both facets random).
%
%Output (struct) matches psyrat_rel_dynrel_trt: ndim/obs/nocc/reltype/ci, z1[/z2],
% G/D/ICCg/ICCd (pt/ll/ul surfaces), and sige_pt (posterior-mean observed-scale
% residual SD sigma_pto,e(z) over the grid).

% Copyright (C) 2016-2026 Peter E. Clayson
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

if mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_rel_dynrel_trt_gamma for more information about inputs'));
end

alpha0      = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b           = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
lognu0      = psyrat_req(varargin,'lognu0','Please input lognu0 (Intercept_nu draws).');
b_nu        = psyrat_req(varargin,'b_nu','Please input b_nu (log-nu dimension slopes).');
sig_p       = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_v       = psyrat_req(varargin,'sig_v','Please input sig_v (person log-nu SD draws).');
rho_pv      = psyrat_req(varargin,'rho_pv','Please input rho_pv (mean<->log-nu corr draws).');
sig_occ     = psyrat_req(varargin,'sig_occ','Please input sig_occ (occasion log-mean SD draws).');
sig_trl     = psyrat_req(varargin,'sig_trl','Please input sig_trl (trial log-mean SD draws).');
sig_trlxid  = psyrat_req(varargin,'sig_trlxid','Please input sig_trlxid (trial x person SD draws).');
sig_occxid  = psyrat_req(varargin,'sig_occxid','Please input sig_occxid (occasion x person SD draws).');
sig_trlxocc = psyrat_req(varargin,'sig_trlxocc','Please input sig_trlxocc (trial x occasion SD draws).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1 (dimension-1 grid axis).');
obs         = psyrat_req(varargin,'obs','Please input obs (n'' trial count).');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

z2      = psyrat_opt(varargin,'z2',[]);
nocc    = psyrat_opt(varargin,'nocc',1);
reltype = psyrat_opt(varargin,'reltype',3);

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end
if ~isscalar(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a positive scalar.\n');
end
if ~isscalar(nocc) || nocc < 1
    error('varargin:nocc','WARNING: nocc must be a positive scalar.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1, 2, or 3.\n');
end

%column vectors of draws
alpha0      = alpha0(:);
lognu0      = lognu0(:);
sig_p       = sig_p(:);
sig_v       = sig_v(:);
rho_pv      = rho_pv(:);
sig_occ     = sig_occ(:);
sig_trl     = sig_trl(:);
sig_trlxid  = sig_trlxid(:);
sig_occxid  = sig_occxid(:);
sig_trlxocc = sig_trlxocc(:);

%log-scale person mean<->log-nu covariance (pass the covariance directly to the
%converter; estimand #2), identical to psyrat_gamma_extract_trt.
cov_pv = rho_pv .* sig_p .* sig_v;

%b/b_nu must be draws x KDIM (KDIM = 1 for one dimension; 3 for two).
kdim = size(b,2);
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if kdim ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        kdim,ndim,expectedk);
end
if size(b_nu,2) ~= expectedk
    error('varargin:bnu',...
        'WARNING: b_nu has %d columns but ndim=%d expects %d.\n',...
        size(b_nu,2),ndim,expectedk);
end

z1 = z1(:)';
nz1 = numel(z1);
if ndim == 2
    z2 = z2(:)';
    nz2 = numel(z2);
else
    nz2 = 1;
end

%preallocate the coefficient surfaces (rows index z1, columns index z2)
gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);
D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);
ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);

for a = 1:nz1
    for c = 1:nz2
        %dimension design row matching the Stan model design
        if ndim == 1
            xrow = z1(a);
        else
            xrow = [z1(a) z2(c) z1(a)*z2(c)];
        end

        %z-conditioned log-mean intercept and population nu at this z
        alpha_z = alpha0 + b * xrow(:);
        nu_z    = exp(lognu0 + b_nu * xrow(:));

        %map the log-scale crossed-design parameters to OBSERVED-score-scale
        %variance components at this z (estimand #2). The argument order follows
        %psyrat_gamma_extract_trt: (alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu, ...)
        %with s_po = sig_occxid, s_pt = sig_trlxid, s_ot = sig_trlxocc.
        vc = psyrat_gamma_varcomps_trt(alpha_z, sig_p, sig_occ, sig_trl, ...
            sig_occxid, sig_trlxid, sig_trlxocc, nu_z, sig_v, cov_pv);

        %SD-equivalents feeding psyrat_rel_trt exactly as psyrat_gamma_extract_trt
        %(bp<-sigma_p, bo<-sigma_o, bt<-sigma_t, txp<-sigma_pt, oxp<-sigma_po,
        %txo<-sigma_ot, err<-sigma_res). psyrat_rel_trt squares its SD inputs.
        bp_z  = sqrt(vc.sigma_p2);
        bo_z  = sqrt(vc.sigma_o2);
        bt_z  = sqrt(vc.sigma_t2);
        txp_z = sqrt(vc.sigma_pt2);
        oxp_z = sqrt(vc.sigma_po2);
        txo_z = sqrt(vc.sigma_ot2);
        err_z = sqrt(vc.sigma_res2);

        sige_pt(a,c) = mean(err_z);

        facetargs = {'bo',bo_z,'bt',bt_z,'txp',txp_z,'oxp',oxp_z,'txo',txo_z};

        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',2,'reltype',reltype,'bp',bp_z,facetargs{:},...
            'err',err_z,'obs',obs,'nocc',nocc,'CI',ciperc);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',1,'reltype',reltype,'bp',bp_z,facetargs{:},...
            'err',err_z,'obs',obs,'nocc',nocc,'CI',ciperc);
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',2,'reltype',reltype,'bp',bp_z,facetargs{:},...
            'err',err_z,'obs',1,'nocc',1,'CI',ciperc);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',1,'reltype',reltype,'bp',bp_z,facetargs{:},...
            'err',err_z,'obs',1,'nocc',1,'CI',ciperc);
    end
end

out = struct();
out.ndim = ndim;
out.obs = obs;
out.nocc = nocc;
out.reltype = reltype;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2
    out.z2 = z2;
end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.G    = gridops.squeeze(G,ndim);
out.D    = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
