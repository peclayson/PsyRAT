function out = psyrat_rel_dodiffdynrel(varargin)
%Dynamic/conditional reliability SURFACE of the four-cell NONCONCURRENT
%difference-of-differences for the GAUSSIAN family, as a function of the
%standardized dimension predictor(s) (analysis 28, family = 'gaussian').
%The Gaussian sibling of PSYRAT_REL_DODIFFDYNREL_GAMMA_LS (no family marker
%in the name, the toolbox convention). Derivation:
%SCIENTIFIC_FORMULA_AUDIT.md section 26.
%
% out = psyrat_rel_dodiffdynrel('b',b,'b_sigma',b_sigma,...
%   'b_dim',b_dim,'b_sigma_dim',b_sigma_dim,'sd_id',sd_id,'cor_id',cor_id,...
%   'sd_trl',sd_trl,'cor_trl',cor_trl,'ndim',ndim,'z1',z1axis,...
%   'obs',nprime,'CI',.95)
%
%THE ESTIMAND is the gamma sibling's, unchanged: the TYPICAL-PERSON read-out
%(every person effect - location AND scale - at zero, so the residual at z is
%the population sigma_q(z) = exp(b_sigma_q + b_sigma_dim_q . X(z)) exactly;
%owner-approved 2026-08-12, carried to the Gaussian arm 2026-08-20). Label any
%display of this surface as the typical-person reference curve; the
%per-participant table is the primary output.
%
%What the identity link changes (section 26): the signal component matrices
%are the model covariance blocks directly and are CONSTANT IN Z - the surface
%varies with z only through the residual diagonal - and the residual channel
%is exactly diagonal, so the relative error reduces to
%sum_q c_q^2*sigma_q(z)^2/n_q. The composite kernel, the harmonic count rule,
%the coefficient formulas, and the nonconcurrent boundary are identical to the
%gamma sibling.
%
%Inputs (name-value; cells in dodmap constituent order)
% b, b_sigma      - [ndraws x 4] per-cell intercepts: b on the NATURAL scale
%                   (identity link - no offsets exist for this family),
%                   b_sigma the ABSOLUTE log residual SD
% b_dim, b_sigma_dim - [ndraws x 4 x KDIM]
% sd_id           - [ndraws x 8], locations on the natural scale; cor_id -
%                   [ndraws x 8 x 8]
% sd_trl          - [ndraws x 4]; cor_trl - [ndraws x 4 x 4]
% ndim            - 1 or 2
% z1 (, z2)       - grid axes of standardized dimension values
% obs             - [1 x 4] per-cell trial counts n'
% CI              - credible interval width in (0,1)
% cvec            - optional signed contrast (default [1 -1 -1 1])
%
%Output struct out - identical shape to the gamma sibling:
% z1 (, z2), G, D, ICCg, ICCd (ll/pt/ul each), sige_pt, invalid_frac,
% estimand = 'typical person (all person effects at zero)'
%
%See also PSYRAT_SSREL_DODIFFDYNREL, PSYRAT_DOD_COMPOSITE,
%PSYRAT_DOD_COMPONENTS, PSYRAT_REL_DODIFFDYNREL_GAMMA_LS.

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
%

for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
    if any(strcmpi(coupled{1},varargin(1:2:end)))
        error('psyrat_rel_dodiffdynrel:noresidualcoupling',... %Error code and associated error
            ['''%s'' is not an input: all six observation-level residual '...
            'covariances are fixed to zero by the nonconcurrent estimand.'],...
            coupled{1});
    end
end

b           = psyrat_req(varargin,'b','Please input b (per-cell natural-scale mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-cell absolute log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (per-cell mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (per-cell log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (joint person SDs, draws x 8).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (joint person correlations, draws x 8 x 8).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-cell trial SDs, draws x 4).');
cor_trl     = psyrat_req(varargin,'cor_trl','Please input cor_trl (cross-cell trial correlations).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1 (grid axis).');
obs         = psyrat_req(varargin,'obs','Please input obs (per-cell trial counts n'').');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

ind = find(strcmpi('cvec',varargin),1);
if ~isempty(ind)
    cvec = varargin{ind+1};
else
    cvec = [1 -1 -1 1];
end

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2
    z2 = psyrat_req(varargin,'z2','Please input z2 (second grid axis).');
end
ndraws = size(b,1);
kdim = 2*ndim - 1;
if size(b_dim,2) ~= 4 || size(b_dim,3) ~= kdim || ...
        size(b_sigma_dim,2) ~= 4 || size(b_sigma_dim,3) ~= kdim
    error('varargin:b_dim',... %Error code and associated error
        'WARNING: b_dim and b_sigma_dim must be [ndraws x 4 x %d] for ndim = %d.\n',...
        kdim,ndim);
end
if numel(obs) ~= 4 || any(obs <= 0) || any(isnan(obs))
    error('varargin:obs',... %Error code and associated error
        'WARNING: obs must be four positive per-cell trial counts.\n');
end
obs = obs(:)';

sd_p = sd_id(:,1:4);
cor_p = cor_id(:,1:4,1:4);

z1 = z1(:)'; nz1 = numel(z1);
if ndim == 2; z2 = z2(:)'; nz2 = numel(z2); else; nz2 = 1; end

gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);  D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);  ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);
invalid_frac = nan(nz1,nz2);
ciedge = (1-ciperc)/2;

for a = 1:nz1
    for c = 1:nz2
        if ndim == 1
            zt = z1(a);
        else
            zt = [z1(a); z2(c); z1(a)*z2(c)];
        end
        zt = zt(:);

        alpha = zeros(ndraws,4);
        lsig = zeros(ndraws,4);
        for q = 1:4
            alpha(:,q) = b(:,q) + reshape(b_dim(:,q,:),ndraws,kdim) * zt;
            lsig(:,q) = b_sigma(:,q) + reshape(b_sigma_dim(:,q,:),ndraws,kdim) * zt;
        end

        %typical person: the population sigma(z), no person effects
        mats = psyrat_dod_components('facet','trial',...
            'alpha',alpha,'log_sigma',lsig,...
            'sd_p',sd_p,'cor_p',cor_p,'sd_trl',sd_trl,'cor_trl',cor_trl);

        u    = psyrat_dod_composite('varcov',mats.p,'cvec',cvec,'rule','unscaled');
        relu = psyrat_dod_composite('varcov',mats.pi_e,'cvec',cvec,'rule','unscaled');
        iu   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','unscaled');
        relh = psyrat_dod_composite('varcov',mats.pi_e,'cvec',cvec,'rule','harmonic','n',obs);
        ih   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','harmonic','n',obs);

        Gd = local_ratio(u.value, u.value + relh.value);
        Dd = local_ratio(u.value, u.value + relh.value + ih.value);
        ig = local_ratio(u.value, u.value + relu.value);
        id_ = local_ratio(u.value, u.value + relu.value + iu.value);

        tol = 1e-10;
        badcoef = @(x) ~isfinite(x) | x < -tol | x > 1 + tol;
        inval = u.invalid_negative | relh.invalid_negative | ...
            ih.invalid_negative | badcoef(Gd) | badcoef(Dd);

        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = local_sum(Gd,ciedge);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = local_sum(Dd,ciedge);
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = local_sum(ig,ciedge);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = local_sum(id_,ciedge);
        sige_pt(a,c) = mean(sqrt(max(relu.value,0)),'omitnan');
        invalid_frac(a,c) = mean(inval);
    end
end

out = struct();
out.ndim = ndim; out.obs = obs; out.ci = ciperc;
out.z1 = z1;
if ndim == 2; out.z2 = z2; end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.invalid_frac = gridops.squeeze1d(invalid_frac,ndim);
out.G = gridops.squeeze(G,ndim);
out.D = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);
out.estimand = 'typical person (all person effects at zero)';

end

function r = local_ratio(num,den)
r = num ./ den;
r(~isfinite(num) | ~isfinite(den) | num < 0 | den <= 0) = NaN;
end

function [ll,pt,ul] = local_sum(draws,ciedge)
q = quantile(draws,[ciedge 1-ciedge]);
pt = mean(draws,'omitnan');
ll = q(1); ul = q(2);
end
