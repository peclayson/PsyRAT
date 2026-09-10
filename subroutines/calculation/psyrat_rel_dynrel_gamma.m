function out = psyrat_rel_dynrel_gamma(varargin)
%Dynamic/conditional ONE-facet reliability for the scaled-chi-square / Gamma(log
%link) family, as a function of the standardized dimension predictor(s) (gamma
%analog of psyrat_rel_dynrel; analysis 11, family = 'gamma').
%
% out = psyrat_rel_dynrel_gamma('alpha0',alpha0,'b',b,'lognu0',lognu0,...
%   'b_nu',b_nu,'sig_p',sig_p,'sig_v',sig_v,'rho_pv',rho_pv,'sig_i',sig_i,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',nprime,'CI',.95)
%
%Unlike the Gaussian dynrel (psyrat_rel_dynrel), the gamma family has no
%log-residual scale submodel. The observed-score residual is the MEAN-COUPLED
%scaled-chi-square dispersion Var(Y|mu,nu) = 2*mu^2/nu, so BOTH the mean and the
%dispersion are dimension-conditioned through the log-linked submodels
%
%   log mu(z) = alpha0 + b . X(z)          (mean at dimension location z)
%   log nu(z) = lognu0 + b_nu . X(z)        (dispersion at z)
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two dimensions. At each z the log-scale parameters are mapped to
%OBSERVED-score-scale variance components with the same lognormal-moment
%converter used by the static one-facet gamma path (psyrat_gamma_varcomps), and
%reliability is the validated one-facet coefficient (psyrat_rel_sing) evaluated on
%those components:
%
%   G(z) = sigma_p^2(z) / (sigma_p^2(z) + sigma_pi,e^2(z) / n')                (Rocha Table 2)
%   D(z) = sigma_p^2(z) / (sigma_p^2(z) + (sigma_pi,e^2(z) + sigma_i^2(z)) / n')
%
%so generalizability excludes, and dependability includes, the trial main
%effect, exactly as the static one-facet gamma path. NOTE that ALL of
%sigma_p^2(z), sigma_i^2(z), and sigma_pi,e^2(z) vary with z here (everything is
%on the multiplicative/log scale), unlike the Gaussian dynrel where only the
%residual varies with z.
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's dispersion heterogeneity (estimand #2), the same convention the
%static gamma family uses: the person log-nu SD (sig_v) and its correlation with
%the person mean (rho_pv) are propagated into the converter as (s_v, cov_pv). This
%is the single gamma reliability estimand; there is no separate typical-person /
%population-average toggle (the marginal argument of psyrat_rel_dynrel), because
%the coupled dispersion is integrated by the law of total variance.
%
%Inputs (name-value; all posterior-draw inputs are column vectors of length
%ndraws, or draws x KDIM for the slope blocks)
% alpha0 - log expected-score grand intercept draws (Intercept). z = 0 is the
%          OVERALL-sample mean of each dimension (standardized once across the
%          whole sample; see psyrat_zscore_subjectlevel).
% b      - log-mean dimension-slope draws (draws x KDIM).
% lognu0 - log degrees-of-freedom grand intercept draws (Intercept_nu).
% b_nu   - log-nu dimension-slope draws (draws x KDIM).
% sig_p  - person log-mean SD draws (gro_sds(:,1) = s_p).
% sig_v  - person log-nu SD draws (gro_sds(:,2) = s_v).
% rho_pv - person mean <-> log-nu correlation draws (chol_corrmat(:,2,1)).
% sig_i  - trial main-effect SD draws on the log mean (sig_trl = s_i).
% ndim   - number of dimensions: 1 or 2.
% z1     - vector of standardized values for dimension 1 (the grid axis).
% obs    - n' (trial count) at which to evaluate the coefficients (scalar).
% CI     - credible-interval width in (0,1], e.g., .95.
%
%Optional Inputs
% z2     - vector of standardized values for dimension 2 (required if ndim=2).
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
        'See help psyrat_rel_dynrel_gamma for more information about inputs'));
end

alpha0 = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b      = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
lognu0 = psyrat_req(varargin,'lognu0','Please input lognu0 (Intercept_nu draws).');
b_nu   = psyrat_req(varargin,'b_nu','Please input b_nu (log-nu dimension slopes).');
sig_p  = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_v  = psyrat_req(varargin,'sig_v','Please input sig_v (person log-nu SD draws).');
rho_pv = psyrat_req(varargin,'rho_pv','Please input rho_pv (mean<->log-nu corr draws).');
sig_i  = psyrat_req(varargin,'sig_i','Please input sig_i (trial log-mean SD draws).');
ndim   = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1     = psyrat_req(varargin,'z1','Please input z1 (dimension-1 grid axis).');
obs    = psyrat_req(varargin,'obs','Please input obs (n'' trial count).');
ciperc = psyrat_req(varargin,'CI','Please input CI.');

z2     = psyrat_opt(varargin,'z2',[]);

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
alpha0 = alpha0(:);
lognu0 = lognu0(:);
sig_p  = sig_p(:);
sig_v  = sig_v(:);
rho_pv = rho_pv(:);
sig_i  = sig_i(:);

%log-scale person mean<->log-nu covariance on the log scale (pass the covariance
%directly to the converter; estimand #2). This is the identical (s_v, cov_pv)
%coupling used by psyrat_gamma_extract_single for the static one-facet path.
cov_pv = rho_pv .* sig_p .* sig_v;

%b/b_nu must be draws x KDIM (KDIM = 1 for one dimension; 3 for two). The
%expected column-count check below catches a transposed/mis-shaped input.
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

        %map the log-scale parameters to OBSERVED-score-scale variance components
        %at this z (estimand #2, marginal over person dispersion via s_v/cov_pv).
        vc = psyrat_gamma_varcomps(alpha_z, sig_p, sig_i, nu_z, sig_v, cov_pv);

        %SD-equivalents in the production one-facet convention: bp = sigma_p(z);
        %total-within wp = sqrt(sigma_pi,e^2(z) + sigma_i^2(z)); i = sigma_i(z).
        %psyrat_rel_sing subtracts i^2 to recover the relative residual for G and
        %keeps it for D, matching the static one-facet gamma / Rocha Table 2 path.
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
