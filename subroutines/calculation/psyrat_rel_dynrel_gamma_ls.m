function out = psyrat_rel_dynrel_gamma_ls(varargin)
%Dynamic/conditional ONE-facet reliability for the LOCATION-SCALE scaled-chi-
%square / Gamma(log link) family, as a function of the standardized dimension
%predictor(s) (analysis 11, family = 'gamma', gammascale = 2).
%
% out = psyrat_rel_dynrel_gamma_ls('alpha0',alpha0,'b',b,'logsig0',logsig0,...
%   'b_sigma',b_sigma,'sig_p',sig_p,'sig_sd',sig_sd,'sig_i',sig_i,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',nprime,'CI',.95)
%
%This is the SIBLING of psyrat_rel_dynrel_gamma, not an overload - the same
%convention psyrat_gamma_varcomps_ls already follows. Both describe the same
%scaled-chi-square likelihood; they differ in which quantity carries the scale
%submodel, and they are DIFFERENT ESTIMANDS whose numbers are not comparable.
%
%   log mu(z)    = alpha0  + b       . X(z)      (mean at dimension location z)
%   log sigma(z) = logsig0 + b_sigma . X(z)      (residual SD at z)
%   nu           = 2*(mu/sigma)^2,  Var(Y | mu,sigma) = sigma^2
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two dimensions. At each z the log-scale parameters are mapped to
%OBSERVED-score-scale variance components by psyrat_gamma_varcomps_ls, and
%reliability is the validated one-facet coefficient (psyrat_rel_sing) evaluated
%on those components:
%
%   G(z) = sigma_p^2(z) / (sigma_p^2(z) + sigma_pi,e^2(z) / n')                (Rocha Table 2)
%   D(z) = sigma_p^2(z) / (sigma_p^2(z) + (sigma_pi,e^2(z) + sigma_i^2(z)) / n')
%
%so generalizability excludes, and dependability includes, the trial main
%effect, exactly as the log-nu sibling and the static one-facet gamma path.
%
%WHY THIS PARAMETERIZATION EXISTS. Under the log-nu sibling the conditional
%variance is 2*mu(z)^2/nu(z), which is MEAN-COUPLED, so a slope on log-nu blends
%a mean-level effect with a dispersion effect and cannot be read as a statement
%about precision alone. Here the residual is decoupled from the mean surface, so
%b_sigma is a pure residual-SD slope. Owner ruling 2026-08-03; see
%documentation/gamma_scale_submodel_decision.md.
%
%THE MEAN SLOPE DOES NOT CANCEL HERE, AND THAT IS THE POINT.
%Under the log-nu sibling EVERY observed-scale term - including the residual -
%carries a common exp(2*alpha(z)) factor, so a mean-level slope b scales
%numerator and denominator alike and divides out of the coefficient; only b_nu
%moves it. That cancellation is documented for gammascale = 1 and verified there
%to ~1e-15. It does NOT hold here. psyrat_gamma_varcomps_ls's e_cond_var is
%exp(2*log_sigma + 2*s_sd^2), which contains no alpha, so writing the relative
%residual as sigma_pi,e^2(z) = exp(2*alpha(z))*(V - A - B) + E gives
%
%   G(z) = A / ( A + (V - A - B)/n' + E/(n' * exp(2*alpha(z))) )
%
%with A, B, V and E free of alpha. The last term shrinks as alpha(z) grows, so a
%POSITIVE mean slope raises reliability with z: a fixed absolute residual SD is
%a smaller share of a larger expected score. Both b and b_sigma therefore move
%the coefficient under gammascale = 2, and any user-facing text carried over
%from the log-nu path ("only the dispersion slope moves it") is false here.
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's residual-SD heterogeneity, E[sigma_p^2] over w_p ~ N(0, s_sd^2),
%so the components sum to the observed-score variance by the law of total
%variance. This is the same marginal convention the log-nu sibling and the
%static gamma family use, and there is no typical-person / population-average
%toggle (the marginal argument of psyrat_rel_dynrel).
%
%NOTE there is deliberately NO rho / cov_pv input, and the log-nu sibling's
%single cov_pv computation outside the grid loop has no counterpart here. The
%person mean <-> log-sigma correlation IS estimated (it is the LKJ block in the
%Stan model) and IS reported in the variance-component table, but the residual
%never touches mu, so that correlation cannot enter the observed-scale
%conversion. psyrat_gamma_varcomps_ls accepts no such argument by design; adding
%one "for symmetry" would be wrong.
%
%Inputs (name-value; all posterior-draw inputs are column vectors of length
%ndraws, or draws x KDIM for the slope blocks)
% alpha0  - log expected-score grand intercept draws (Intercept). z = 0 is the
%           OVERALL-sample mean of each dimension (standardized once across the
%           whole sample; see psyrat_zscore_subjectlevel).
% b       - log-mean dimension-slope draws (draws x KDIM).
% logsig0 - log residual-SD grand intercept draws (Intercept_sigma).
% b_sigma - log residual-SD dimension-slope draws (draws x KDIM).
% sig_p   - person log-mean SD draws (gro_sds(:,1) = s_p).
% sig_sd  - person log-sigma SD draws (gro_sds(:,2) = s_sd). NOTE this is the
%           SAME SLOT that holds the person log-nu SD (s_v) under gammascale = 1
%           and a DIFFERENT quantity; passing one where the other belongs
%           produces numbers, not an error.
% sig_i   - trial main-effect SD draws on the log mean (sig_trl = s_i).
% ndim    - number of dimensions: 1 or 2.
% z1      - vector of standardized values for dimension 1 (the grid axis).
% obs     - n' (trial count) at which to evaluate the coefficients (scalar).
% CI      - credible-interval width in (0,1], e.g., .95.
%
%Optional Inputs
% z2      - vector of standardized values for dimension 2 (required if ndim=2).
%
%Output (struct)
% out.ndim, out.obs, out.ci
% out.z1            - dimension-1 grid axis (1 x nz1)
% out.z2            - dimension-2 grid axis (1 x nz2; only when ndim=2)
% out.G, out.D      - structs with fields pt/ll/ul (global coefficients over the
%                     grid). ndim=1 -> 1 x nz1; ndim=2 -> nz1 x nz2 (rows index
%                     z1, columns index z2).
% out.ICCg, out.ICCd- relative/absolute single-observation ICC over the grid.
% out.sige_pt       - posterior-mean OBSERVED-scale residual SD sigma_pi,e(z) over
%                     the grid (for reporting the residual surface).
%
%See also PSYRAT_REL_DYNREL_GAMMA, PSYRAT_GAMMA_VARCOMPS_LS,
%PSYRAT_SSREL_DYNREL_GAMMA_LS, PSYRAT_REL_SING.

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
        'See help psyrat_rel_dynrel_gamma_ls for more information about inputs'));
end

alpha0  = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b       = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
logsig0 = psyrat_req(varargin,'logsig0','Please input logsig0 (Intercept_sigma draws).');
b_sigma = psyrat_req(varargin,'b_sigma','Please input b_sigma (log residual-SD dimension slopes).');
sig_p   = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_sd  = psyrat_req(varargin,'sig_sd','Please input sig_sd (person log-sigma SD draws).');
sig_i   = psyrat_req(varargin,'sig_i','Please input sig_i (trial log-mean SD draws).');
ndim    = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1      = psyrat_req(varargin,'z1','Please input z1 (dimension-1 grid axis).');
obs     = psyrat_req(varargin,'obs','Please input obs (n'' trial count).');
ciperc  = psyrat_req(varargin,'CI','Please input CI.');

z2      = psyrat_opt(varargin,'z2',[]);

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end
if ~isscalar(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a positive scalar.\n');
end

%column vectors of draws
alpha0  = alpha0(:);
logsig0 = logsig0(:);
sig_p   = sig_p(:);
sig_sd  = sig_sd(:);
sig_i   = sig_i(:);

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

        %map the log-scale parameters to OBSERVED-score-scale variance components
        %at this z, marginal over the person residual-SD population (s_sd). No
        %cov argument exists: the residual never touches mu.
        vc = psyrat_gamma_varcomps_ls(alpha_z, sig_p, sig_i, logsig_z, sig_sd);

        %SD-equivalents in the production one-facet convention: bp = sigma_p(z);
        %total-within wp = sqrt(sigma_pi,e^2(z) + sigma_i^2(z)); i = sigma_i(z).
        %psyrat_rel_sing subtracts i^2 to recover the relative residual for G and
        %keeps it for D, matching the log-nu sibling and Rocha Table 2.
        bp_z = sqrt(vc.sigma_p2);
        i_z  = sqrt(vc.sigma_i2);
        wp_z = sqrt(vc.sigma_pi_e2 + vc.sigma_i2);

        sige_pt(a,c) = mean(sqrt(vc.sigma_pi_e2));

        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',2,'metric','global','bp',bp_z,'wp',wp_z,'i',i_z,...
            'obs',obs,'CI',ciperc);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',1,'metric','global','bp',bp_z,'wp',wp_z,'i',i_z,...
            'obs',obs,'CI',ciperc);
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',2,'metric','icc','bp',bp_z,'wp',wp_z,'i',i_z,'CI',ciperc);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',1,'metric','icc','bp',bp_z,'wp',wp_z,'i',i_z,'CI',ciperc);
    end
end

out = struct();
out.ndim = ndim;
out.obs = obs;
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
