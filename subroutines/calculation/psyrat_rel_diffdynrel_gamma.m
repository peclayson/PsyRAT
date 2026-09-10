function out = psyrat_rel_diffdynrel_gamma(varargin)
%Dynamic/conditional reliability of a two-event DIFFERENCE score for the
%scaled-chi-square / Gamma(log link) family, as a function of the standardized
%dimension predictor(s) (gamma analog of psyrat_rel_diffdynrel; analysis 12,
%family = 'gamma', group-level, NON-concurrent).
%
% out = psyrat_rel_diffdynrel_gamma('b',b,'b_nu',b_nu,'b_dim',b_dim,...
%   'b_nu_dim',b_nu_dim,'sd_id',sd_id,'sd_trl',sd_trl,'cor_id',cor_id,...
%   'cor_i',cor_i,'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',nprime,'CI',.95)
%
%Unlike the Gaussian dynamic difference (psyrat_rel_diffdynrel), the gamma family
%has no log-residual scale submodel. The observed-score residual is the
%MEAN-COUPLED scaled-chi-square dispersion Var(Y|mu,nu) = 2*mu^2/nu, so BOTH the
%per-event mean and the per-event dispersion are dimension-conditioned through the
%log-linked submodels
%
%   log mu_k(z) = b_k    + b_dim_k    . X(z)      (event k mean at z)
%   log nu_k(z) = b_nu_k + b_nu_dim_k . X(z)      (event k dispersion at z)
%
%where X(z) is the dimension-term vector: [z1] for one dimension; [z1 z2 z1*z2]
%for two. The event-crossed slope blocks follow the case-12 Xdim convention -
%ODD columns are event 1, EVEN columns event 2.
%
%KEY DIFFERENCE FROM THE GAUSSIAN PATH. psyrat_rel_diffdynrel holds the person and
%trial (co)variance matrices FIXED across z and varies only the residual, because
%under an identity link the mean slopes do not move the variance components. Here
%every observed-scale component carries exp(2*alpha(z) + V), so sigma^2_p(z),
%sigma^2_i(z), the cross-event covariances AND the residual all vary with z. This
%function therefore REBUILDS id_varcov(z) / trl_varcov(z) at every grid point from
%the raw log-scale draws, using the same converters as the static gamma difference
%(psyrat_gamma_varcomps for the per-event diagonals, psyrat_gamma_crosscov for the
%cross-event off-diagonals), and then reuses the validated psyrat_diffrel:
%
%   sigma^2_p,Delta(z) = sigma^2_p,A(z) + sigma^2_p,B(z) - 2*cov_p,AB(z)
%   G_Delta(z) = sigma^2_p,Delta(z) / (sigma^2_p,Delta(z) + relative error(z))
%   D_Delta(z) = G_Delta(z) with the trial contrast variance added to the error
%
%exactly the Rocha Table 6 one-facet difference shape, z-conditioned. Residual
%cross-covariance is 0 (non-concurrent), matching psyrat_gamma_crosscov and the
%static case-7 gamma difference path.
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's dispersion heterogeneity (estimand #2), the same convention the rest
%of the gamma family uses: each event's person log-nu SD and its correlation with
%that event's person mean are propagated into the converter as (s_v, cov_pv).
%Those are log-scale quantities and therefore z-INDEPENDENT, so they are formed
%once outside the grid loop.
%
%Inputs (name-value; all posterior-draw inputs are draws x N matrices)
% b        - per-event log-mean intercepts at z = 0 (draws x 2). z = 0 is the
%            OVERALL-sample mean of each dimension (standardized once across the
%            whole sample; see psyrat_zscore_subjectlevel).
% b_nu     - per-event log-nu intercepts at z = 0 (draws x 2).
% b_dim    - event-crossed log-mean dimension slopes (draws x KDIM; KDIM = 2 for
%            one dimension, 6 for two). Odd columns event 1, even columns event 2.
% b_nu_dim - event-crossed log-nu dimension slopes (draws x KDIM), same layout.
% sd_id    - person log-scale SDs (draws x 4): [mean_1, mean_2, lognu_1, lognu_2].
% sd_trl   - per-event trial main-effect log-mean SDs (draws x 2).
% cor_id   - person correlation matrix draws (draws x 4 x 4). Supplies the
%            cross-event person-mean correlation cor_id(:,1,2) and the WITHIN-event
%            mean<->log-nu correlations cor_id(:,1,3) / cor_id(:,2,4).
% cor_i    - cross-event trial correlation on the log scale (draws x 1).
% ndim     - number of dimensions: 1 or 2.
% z1       - grid axis for dimension 1.
% obs      - n' per event. Scalar (applied to both events) or [n1 n2]. Difference
%            designs are routinely unbalanced, so prefer the two-element form.
% CI       - credible-interval width in (0,1].
%
%Optional Inputs
% z2       - grid axis for dimension 2 (required if ndim = 2).
%
%Output (struct)
% out.ndim, out.obs, out.ci
% out.z1 (, out.z2)   - grid axes
% out.G, out.D        - difference-score surfaces (ll/pt/ul); ndim=1 -> 1 x nz1;
%                       ndim=2 -> nz1 x nz2 (rows index z1, columns index z2).
% out.ICCg, out.ICCd  - relative/absolute single-observation difference ICC.
% out.sige_pt         - posterior-mean residual CONTRAST SD sqrt(sigma^2_pi,e,Delta)
%                       at each grid point (for reporting the residual surface).
%
% See also PSYRAT_REL_DIFFDYNREL, PSYRAT_REL_DYNREL_GAMMA, PSYRAT_GAMMA_VARCOMPS,
% PSYRAT_GAMMA_CROSSCOV, PSYRAT_DIFFREL.

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

if mod(length(varargin),2)
    error('varargin:incomplete',...
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_rel_diffdynrel_gamma for more information about inputs'));
end

b        = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_nu     = psyrat_req(varargin,'b_nu','Please input b_nu (per-event log-nu intercepts).');
b_dim    = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_nu_dim = psyrat_req(varargin,'b_nu_dim','Please input b_nu_dim (log-nu dimension slopes).');
sd_id    = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl   = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
cor_id   = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_i    = psyrat_req(varargin,'cor_i','Please input cor_i (cross-event trial correlation).');
ndim     = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1       = psyrat_req(varargin,'z1','Please input z1.');
obs      = psyrat_req(varargin,'obs','Please input obs (n'').');
ciperc   = psyrat_req(varargin,'CI','Please input CI.');
z2       = psyrat_opt(varargin,'z2',[]);

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end

%n' per event. A scalar is broadcast to both events (balanced design); a
%two-element vector carries genuinely unbalanced per-event replication straight
%through to psyrat_diffrel, which already scales each event's variance by its own
%n'.
if isscalar(obs)
    obs = [obs obs];
end
%~isfinite is tested explicitly: NaN < 1 is FALSE in MATLAB, so a NaN n' (e.g. a
%median over an empty per-event trial count) would pass a bare `any(obs < 1)` and
%then divide every error term by NaN, yielding an all-NaN surface and an empty
%plot with no diagnostic.
if numel(obs) ~= 2 || any(~isfinite(obs)) || any(obs < 1)
    error('varargin:obs',...
        'WARNING: obs (n'') must be a finite positive scalar or a 2-element [n1 n2].\n');
end

%expected number of event-crossed slope columns: 2 events x (1 or 3 dim terms)
expectedk = (ndim == 1) * 2 + (ndim == 2) * 6;
if size(b_dim,2) ~= expectedk
    error('varargin:bdim',...
        'WARNING: b_dim has %d columns but ndim=%d expects %d.\n',...
        size(b_dim,2),ndim,expectedk);
end
if size(b_nu_dim,2) ~= expectedk
    error('varargin:bnudim',...
        'WARNING: b_nu_dim has %d columns but ndim=%d expects %d.\n',...
        size(b_nu_dim,2),ndim,expectedk);
end
if size(sd_id,2) ~= 4
    error('varargin:sdid',...
        ['WARNING: sd_id must have 4 columns '...
        '[mean_1, mean_2, lognu_1, lognu_2]; got %d.\n'],size(sd_id,2));
end

%event-crossed slopes: odd columns -> event 1, even columns -> event 2, matching
%the case-12 Xdim construction [meas1.*z1, meas2.*z1, meas1.*z2, meas2.*z2, ...].
bdim1 = b_dim(:,1:2:end);      %draws x (KDIM/2)
bdim2 = b_dim(:,2:2:end);
bnud1 = b_nu_dim(:,1:2:end);
bnud2 = b_nu_dim(:,2:2:end);

%per-event person mean SDs, trial SDs, and person log-nu SDs (all log scale, all
%z-independent).
sp_a = sd_id(:,1); sp_b = sd_id(:,2);
sv_a = sd_id(:,3); sv_b = sd_id(:,4);
si_a = sd_trl(:,1); si_b = sd_trl(:,2);

%WITHIN-event person mean<->log-nu covariances carry estimand #2 per event; the
%CROSS-event person-mean correlation drives the difference numerator. Both are
%log-scale, hence z-independent, so form them once. (Cross-event log-nu
%correlations do not enter a cell's own residual - see psyrat_gamma_varcomps.)
cov_pv_a = cor_id(:,1,3) .* sp_a .* sv_a;
cov_pv_b = cor_id(:,2,4) .* sp_b .* sv_b;
cor_p    = cor_id(:,1,2);
cor_i    = cor_i(:);

z1 = z1(:)'; nz1 = numel(z1);
if ndim == 2; z2 = z2(:)'; nz2 = numel(z2); else; nz2 = 1; end

gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);  D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);  ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);

ndraws = size(b,1);
id_varcov  = zeros(ndraws,2,2);
trl_varcov = zeros(ndraws,2,2);

%One warning per surface, not one per grid point: the kernel warns on every
%call and this loop calls it nz1*nz2 times. Suppress here, aggregate, report
%once. The gamma difference variant deliberately passes UNEQUAL per-event trial
%counts at every point, so without this a single surface would emit hundreds of
%identical blocks.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for a = 1:nz1
    for c = 1:nz2
        %dimension-term vector matching the case-12 design
        if ndim == 1
            zterms = z1(a);
        else
            zterms = [z1(a); z2(c); z1(a)*z2(c)];
        end
        zt = zterms(:);

        %z-conditioned per-event log-mean intercepts and population nu
        alpha_a = b(:,1)    + bdim1 * zt;
        alpha_b = b(:,2)    + bdim2 * zt;
        nu_a    = exp(b_nu(:,1) + bnud1 * zt);
        nu_b    = exp(b_nu(:,2) + bnud2 * zt);

        %OBSERVED-score-scale components at this z (estimand #2 per event) and the
        %cross-event observed-scale covariances. Both converters are the ones the
        %STATIC gamma difference (case 7) uses; only the intercepts move with z.
        vc_a = psyrat_gamma_varcomps(alpha_a, sp_a, si_a, nu_a, sv_a, cov_pv_a);
        vc_b = psyrat_gamma_varcomps(alpha_b, sp_b, si_b, nu_b, sv_b, cov_pv_b);
        cc   = psyrat_gamma_crosscov(alpha_a, sp_a, si_a, ...
                                     alpha_b, sp_b, si_b, cor_p, cor_i);

        %assemble the 2x2 person and trial (co)variance arrays psyrat_diffrel
        %expects, exactly as psyrat_gamma_extract_diff does for the static case.
        id_varcov(:,1,1) = vc_a.sigma_p2;  id_varcov(:,2,2) = vc_b.sigma_p2;
        id_varcov(:,1,2) = cc.cov_p_obs;   id_varcov(:,2,1) = cc.cov_p_obs;
        trl_varcov(:,1,1) = vc_a.sigma_i2; trl_varcov(:,2,2) = vc_b.sigma_i2;
        trl_varcov(:,1,2) = cc.cov_i_obs;  trl_varcov(:,2,1) = cc.cov_i_obs;

        %psyrat_diffrel reads the residual as a log-scale SD and squares it
        %(exp(er_var).^2 = sigma_pi_e2), so encode as 0.5*log(sigma_pi_e2) -
        %the identical encoding used by psyrat_gamma_extract_diff.
        er_var_z = [0.5*log(vc_a.sigma_pi_e2(:)), 0.5*log(vc_b.sigma_pi_e2(:))];

        %residual cross-covariance is 0 (non-concurrent), matching
        %psyrat_gamma_crosscov's cov_e_obs and the static case-7 gamma path.
        g = psyrat_diffrel('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
            'obs',obs,'CI',ciperc,'est','gen','wp_cov',0);
        d = psyrat_diffrel('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
            'obs',obs,'CI',ciperc,'est','dep','wp_cov',0);

        G.ll(a,c) = g.ll; G.pt(a,c) = g.pt; G.ul(a,c) = g.ul;
        D.ll(a,c) = d.ll; D.pt(a,c) = d.pt; D.ul(a,c) = d.ul;
        ICCg.ll(a,c) = g.icc_ll; ICCg.pt(a,c) = g.icc_pt; ICCg.ul(a,c) = g.icc_ul;
        ICCd.ll(a,c) = d.icc_ll; ICCd.pt(a,c) = d.icc_pt; ICCd.ul(a,c) = d.icc_ul;

        if isfield(g,'admissibility')
            adm_all = [adm_all, g.admissibility, d.admissibility]; %#ok<AGROW>
        end

        %residual CONTRAST SD at this z (non-concurrent, so the cross term is 0)
        sige_pt(a,c) = mean(sqrt(vc_a.sigma_pi_e2 + vc_b.sigma_pi_e2));
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

out = struct();
out.ndim = ndim; out.obs = obs; out.ci = ciperc;
out.z1 = z1;
if ndim == 2; out.z2 = z2; end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.G = gridops.squeeze(G,ndim);
out.D = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
