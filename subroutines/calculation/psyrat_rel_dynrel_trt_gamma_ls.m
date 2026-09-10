function out = psyrat_rel_dynrel_trt_gamma_ls(varargin)
%Dynamic/conditional trial + occasion TWO-facet reliability for the LOCATION-
%SCALE scaled-chi-square / Gamma(log link) family, as a function of the
%standardized dimension predictor(s) (analysis 14, family = 'gamma',
%gammascale = 2).
%
% out = psyrat_rel_dynrel_trt_gamma_ls('alpha0',alpha0,'b',b,...
%   'logsig0',logsig0,'b_sigma',b_sigma,'sig_p',sig_p,'sig_sd',sig_sd,...
%   'sig_occ',sig_occ,'sig_trl',sig_trl,'sig_trlxid',sig_trlxid,...
%   'sig_occxid',sig_occxid,'sig_trlxocc',sig_trlxocc,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',ntrl,'nocc',nocc,...
%   'reltype',reltype,'CI',.95)
%
%This is the SIBLING of psyrat_rel_dynrel_trt_gamma, not an overload - the same
%convention psyrat_rel_dynrel_gamma_ls already follows for the one-facet design.
%Both describe the same scaled-chi-square likelihood; they differ in which
%quantity carries the scale submodel, and they are DIFFERENT ESTIMANDS whose
%numbers are not comparable.
%
%   log mu(z)    = alpha0  + b       . X(z)   (mean at dimension location z)
%   log sigma(z) = logsig0 + b_sigma . X(z)   (residual SD at z)
%   nu           = 2*(mu/sigma)^2,  Var(Y | mu,sigma) = sigma^2
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two dimensions. At each z the log-scale crossed-design parameters are mapped
%to OBSERVED-score-scale variance components by psyrat_gamma_varcomps_trt_ls, and
%the two-facet coefficient is the validated test-retest psyrat_rel_trt evaluated
%on those components. psyrat_rel_trt's three reltypes give:
%
%   reltype 1 (coefficient of equivalence / internal consistency: occasion fixed)
%   reltype 2 (coefficient of stability: trial fixed)
%   reltype 3 (coefficient of trial equivalence and stability: both random)
%
%Generalizability (gcoeff 2) excludes, and dependability (gcoeff 1) includes, the
%facet main effects, exactly as the log-nu sibling and the static test-retest
%gamma path.
%
%WHY THIS PARAMETERIZATION EXISTS. Under the log-nu sibling the conditional
%variance is 2*mu(z)^2/nu(z), which is MEAN-COUPLED, so a slope on log-nu blends
%a mean-level effect with a dispersion effect and cannot be read as a statement
%about precision alone. Here the residual is decoupled from the mean surface, so
%b_sigma is a pure ABSOLUTE residual-SD slope. Since
%log sigma(z) = log mu(z) + 0.5*log(2) - 0.5*log nu(z), the two slope sets are
%related by b_nu = 2*b - 2*b_sigma; the log-nu model is a coherent constant-CV
%(proportional-error) model, not an invalid one. Owner ruling H; see
%documentation/gamma_scale_submodel_decision.md.
%
%THE MEAN SLOPE DOES NOT CANCEL HERE, AND THAT IS THE POINT. Under the log-nu
%sibling every observed-scale component - including the residual - carries a
%common exp(2*alpha(z)) factor, so a mean-level slope scales numerator and
%denominator alike and divides out of the coefficient; only b_nu moves it. That
%cancellation is documented for gammascale = 1 and does NOT hold here:
%psyrat_gamma_varcomps_trt_ls's e_cond_var is exp(2*log_sigma + 2*s_sd^2), which
%contains no alpha, while all six signal components scale with exp(2*alpha(z)).
%A POSITIVE mean slope therefore raises reliability with z - a fixed absolute
%residual SD is a smaller share of a larger expected score. Both b and b_sigma
%move the coefficient under gammascale = 2, and any user-facing text carried over
%from the log-nu path ("only the dispersion slope moves it") is false here.
%
%WHAT THE FIVE CROSSED FACETS DO AND DO NOT DO. Occasion, trial, and the three
%two-way interactions live on the MEAN submodel in both parameterizations, and
%the emitted Stan model's mean half is byte-identical across the two. They
%therefore enter the six signal components identically and are not part of what
%this parameterization changes. There are deliberately no scale-submodel
%counterparts: the residual is the highest-order p x t x o term, confounded with
%pure error at one replicate per cell, exactly as in the log-nu sibling and the
%Gaussian psyrat_build_stan_dynrel_trt.
%
%PRECISELY WHAT "DECOUPLED FROM THE MEAN SURFACE" MEANS HERE, because the phrase
%is easy to over-read. In the fitted Stan model nu is DERIVED as 2*(mu/sigma)^2
%after mu has absorbed all five facets, and the density is evaluated at (mu, nu).
%So the shape parameter is facet- and mean-dependent by construction. What is
%decoupled is the conditional VARIANCE:
%
%   Var(Y | mu, sigma) = 2*mu^2/nu = sigma^2      exactly, no mu, no facets
%
%and the conditional variance is the ONLY moment this file consumes -- e_cond_var
%feeds sigma_res2 and nothing else reads the conditional distribution. That is why
%the decoupling claim is sound for every quantity reported here. It does NOT
%extend to higher moments: conditional skewness is sqrt(8/nu) = 2*sigma/mu, so a
%facet or a mean slope does move the SHAPE of the observation distribution. Do not
%restate this as "the facets/mean do not affect the residual distribution" -- the
%correct statement is that they do not affect its VARIANCE.
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's residual-SD heterogeneity, E[sigma_p^2] over w_p ~ N(0, s_sd^2),
%so the components sum to the observed-score variance by the law of total
%variance. This is the same marginal convention the log-nu sibling and the static
%gamma family use, and there is no typical-person / population-average toggle.
%
%NOTE there is deliberately NO rho / cov_pv input, exactly as in the one-facet
%location-scale sibling, and the log-nu version's cov_pv computation has no
%counterpart here. The person mean <-> log-sigma correlation IS estimated (it is
%the LKJ block in the Stan model) and IS reported in the variance-component
%table, but the residual never touches mu, so that correlation cannot enter the
%observed-scale conversion. psyrat_gamma_varcomps_trt_ls accepts no such argument
%by design; adding one "for symmetry" would be wrong.
%
%THIS IS AN EXTRAPOLATION, NOT AN IMPLEMENTATION OF RAST & CLAYSON. That paper's
%location-scale model is Gaussian throughout (family = gaussian(); no link
%function, no non-normal outcome anywhere in it). The log link, the observed-scale
%signal transformation, and the derived shape nu = 2*(mu/sigma)^2 mean the
%coefficient produced here is NOT behaviorally identical to the paper's. Cite the
%paper for the idea; do not cite it for this model's behavior.
%
%Inputs (name-value; all posterior-draw inputs are column vectors of length
%ndraws, or draws x KDIM for the slope blocks)
% alpha0     - log expected-score grand intercept draws (Intercept). z = 0 is the
%              OVERALL-sample mean of each dimension (standardized once across the
%              whole sample; see psyrat_zscore_subjectlevel).
% b          - log-mean dimension-slope draws (draws x KDIM).
% logsig0    - log residual-SD grand intercept draws (Intercept_sigma).
% b_sigma    - log residual-SD dimension-slope draws (draws x KDIM).
% sig_p      - person log-mean SD draws (gro_sds(:,1) = s_p).
% sig_sd     - person log-sigma SD draws (gro_sds(:,2) = s_sd). NOTE this is the
%              SAME SLOT that holds the person log-nu SD (s_v) under
%              gammascale = 1 and a DIFFERENT quantity; passing one where the
%              other belongs produces numbers, not an error.
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
%Output (struct) matches psyrat_rel_dynrel_trt_gamma: ndim/obs/nocc/reltype/ci,
% z1[/z2], G/D/ICCg/ICCd (pt/ll/ul surfaces), and sige_pt (posterior-mean
% observed-scale residual SD sigma_pto,e(z) over the grid).
%
%See also PSYRAT_REL_DYNREL_TRT_GAMMA, PSYRAT_GAMMA_VARCOMPS_TRT_LS,
%PSYRAT_SSREL_DYNREL_TRT_GAMMA_LS, PSYRAT_REL_DYNREL_GAMMA_LS, PSYRAT_REL_TRT.

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
        'See help psyrat_rel_dynrel_trt_gamma_ls for more information about inputs'));
end

alpha0      = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b           = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
logsig0     = psyrat_req(varargin,'logsig0','Please input logsig0 (Intercept_sigma draws).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (log residual-SD dimension slopes).');
sig_p       = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_sd      = psyrat_req(varargin,'sig_sd','Please input sig_sd (person log-sigma SD draws).');
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
logsig0     = logsig0(:);
sig_p       = sig_p(:);
sig_sd      = sig_sd(:);
sig_occ     = sig_occ(:);
sig_trl     = sig_trl(:);
sig_trlxid  = sig_trlxid(:);
sig_occxid  = sig_occxid(:);
sig_trlxocc = sig_trlxocc(:);

%b/b_sigma must be draws x KDIM (KDIM = 1 for one dimension; 3 for two). The
%expected column-count check below catches a transposed/mis-shaped input.
kdim = size(b,2);
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if kdim ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        kdim,ndim,expectedk);
end
if size(b_sigma,2) ~= expectedk
    error('varargin:bsigma',...
        'WARNING: b_sigma has %d columns but ndim=%d expects %d.\n',...
        size(b_sigma,2),ndim,expectedk);
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

        %z-conditioned log-mean intercept and population log residual SD at z.
        %Note both submodels are evaluated here: unlike the log-nu sibling, the
        %mean slope does not cancel downstream (see the header).
        alpha_z  = alpha0  + b       * xrow(:);
        logsig_z = logsig0 + b_sigma * xrow(:);

        %map the log-scale crossed-design parameters to OBSERVED-score-scale
        %variance components at this z, marginal over the person residual-SD
        %population (s_sd). The argument order follows psyrat_gamma_varcomps_trt:
        %(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, ...) with s_po = sig_occxid,
        %s_pt = sig_trlxid, s_ot = sig_trlxocc. Only the trailing dispersion
        %arguments differ from the log-nu call: (nu, s_v, cov_pv) becomes
        %(log_sigma, s_sd). No cov argument exists: the residual never touches mu.
        %
        %*** SLOT 8 CHANGES SCALE, NOT JUST NAME. DO NOT ADD exp() HERE. ***
        %The log-nu sibling passes nu_z = EXP(lognu0 + b_nu*xrow), because slot 8
        %of psyrat_gamma_varcomps_trt is nu, a NATURAL-scale precision. Slot 8 of
        %psyrat_gamma_varcomps_trt_ls is log_sigma, a LOG-scale quantity, so
        %logsig_z is passed UNEXPONENTIATED - the converter exponentiates it
        %internally (grep e_cond_var = exp(2.*log_sigma ...)). Wrapping this in
        %exp() "for symmetry with the log-nu call" would produce finite, plausible,
        %silently WRONG numbers rather than an error. The two calls look parallel
        %and are not.
        vc = psyrat_gamma_varcomps_trt_ls(alpha_z, sig_p, sig_occ, sig_trl, ...
            sig_occxid, sig_trlxid, sig_trlxocc, logsig_z, sig_sd);

        %SD-equivalents feeding psyrat_rel_trt exactly as the log-nu sibling
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
