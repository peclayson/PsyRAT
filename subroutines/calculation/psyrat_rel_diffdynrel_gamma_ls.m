function out = psyrat_rel_diffdynrel_gamma_ls(varargin)
%Dynamic/conditional reliability of a two-event DIFFERENCE score for the
%LOCATION-SCALE scaled-chi-square / Gamma(log link) family, as a function of the
%standardized dimension predictor(s) (analysis 12, family = 'gamma',
%gammascale = 2, group-level, NON-concurrent).
%
% out = psyrat_rel_diffdynrel_gamma_ls('b',b,'b_sigma',b_sigma,'b_dim',b_dim,...
%   'b_sigma_dim',b_sigma_dim,'sd_id',sd_id,'sd_trl',sd_trl,'cor_id',cor_id,...
%   'cor_i',cor_i,'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',nprime,'CI',.95)
%
%This is the location-scale sibling of PSYRAT_REL_DIFFDYNREL_GAMMA. Both describe
%the same scaled-chi-square likelihood; they differ in which quantity carries the
%dispersion. Here the per-event residual is modeled DIRECTLY as a log residual SD
%with its own dimension slopes,
%
%   log mu_k(z)    = b_k       + b_dim_k       . X(z)   (event k mean at z)
%   log sigma_k(z) = b_sigma_k + b_sigma_dim_k . X(z)   (event k residual SD at z)
%   nu_k           = 2*(mu_k/sigma_k)^2   =>  Var(Y | mu, sigma) = sigma_k^2
%
%where X(z) is the dimension-term vector: [z1] for one dimension; [z1 z2 z1*z2]
%for two. The event-crossed slope blocks follow the case-12 Xdim convention -
%ODD columns are event 1, EVEN columns event 2.
%
%WHY THIS EXISTS. Under the log-nu sibling the conditional variance is
%2*mu^2/nu, which is MEAN-COUPLED, so a dimension slope on log-nu blends a
%mean-level effect with a dispersion effect and cannot be read as a statement
%about precision alone. Here the CONDITIONAL variance is exactly sigma_k(z)^2 and
%no longer moves with the amplitude, so b_sigma_dim is a pure residual-SD slope
%per event. The log-nu sibling is NOT invalid - it is a coherent constant-CV
%model. See documentation/gamma_scale_submodel_decision.md sections B and I.
%
%BE PRECISE ABOUT WHICH RESIDUAL IS DECOUPLED. It is tempting to write "the
%residual is decoupled from the mean surface entirely". That is WRONG and was
%caught by a test written from it. Only vc.e_cond_var is decoupled. The reported
%component vc.sigma_pi_e2 is recovered by SUBTRACTION, so it confounds the person
%x trial interaction of the lognormal MEAN surface with the conditional variance,
%and every mean-surface term scales as exp(2*alpha):
%
%   sigma_pi_e2(z) = exp(2*alpha(z))*K + exp(2*log_sigma(z) + 2*s_sd^2)
%
%with K free of alpha (K = 0 iff s_p = 0 or s_i = 0). So out.sige_pt can vary
%with z through b_dim even at b_sigma_dim = 0. Pinned by
%TestGammaDiffDynrelLsConversion.
%
%WHAT IS AND IS NOT SHARED WITH THE LOG-NU SIBLING. The MEAN submodel is
%identical, so the per-event mean-surface components (sigma_p2, sigma_i2) and
%BOTH cross-event covariances are computed by exactly the same code. Only the
%residual differs. Concretely, per z this function rebuilds id_varcov(z) /
%trl_varcov(z) using psyrat_gamma_varcomps_LS for the per-event diagonals and the
%UNCHANGED psyrat_gamma_crosscov for the off-diagonals, then reuses the validated
%psyrat_diffrel:
%
%   sigma^2_p,Delta(z) = sigma^2_p,A(z) + sigma^2_p,B(z) - 2*cov_p,AB(z)
%   G_Delta(z) = sigma^2_p,Delta(z) / (sigma^2_p,Delta(z) + relative error(z))
%   D_Delta(z) = G_Delta(z) with the trial contrast variance added to the error
%
%WHY PSYRAT_GAMMA_CROSSCOV IS REUSED UNCHANGED, derived rather than assumed. All
%eight of its arguments (per-event alpha, s_p, s_i plus cor_p/cor_i) are
%MEAN-submodel quantities; no dispersion quantity enters it under either
%parameterization. The reason that is legitimate rather than coincidental: the
%cross-event covariance is a covariance of two conditional MEANS, and
%E[Y_k | u, w] = mu_k whatever the dispersion parameterization, because the scale
%parameter governs the second moment of Y given mu, not the first. See
%gamma_scale_submodel_decision.md section I.3 for the full argument.
%
%NO cov_pv COUNTERPART, AND ADDING ONE WOULD BE WRONG. The log-nu sibling threads
%the within-event person mean<->log-nu covariance into its converter because its
%residual carries exp(-2*cov_pv). Here the residual never touches mu, so the
%person mean<->log-sigma correlation cannot affect E[sigma_p^2]. The same argument
%kills the CROSS-event mean<->scale correlations cor_id(1,4) / cor_id(2,3): they
%enter no component of this design. All are still ESTIMATED (the LKJ block) and
%extracted into REL.out.cor_id; they simply do not enter the conversion. Of them,
%only the two WITHIN-event correlations cor_id(1,3) / cor_id(2,4) are surfaced to
%the user, by local_gamma_diff_varcomp_table_ls. The cross-event pair (1,4)/(2,3)
%and the cross-event scale correlation (3,4) are estimated but not reported
%anywhere - do not describe them as reported.
%
%ESTIMAND. The observation-variance term is MARGINALIZED over the person
%population's residual-SD heterogeneity, per event, via the s_sd argument of
%psyrat_gamma_varcomps_ls (E[sigma_p^2] = exp(2*log_sigma(z) + 2*s_sd^2)). NOTE
%log_sigma_k(z) IS z-dependent - it carries b_sigma_dim - so it is rebuilt inside
%the grid loop; only s_sd is z-free and formed once outside it. (The log-nu
%sibling's corresponding sentence says its (s_v, cov_pv) pair is z-independent,
%which is true THERE; do not carry that wording across.)
%
%NOTE this marginalization has no counterpart in the owner's reference bundle,
%which plugs in each participant's realized scale effect instead; live CmdStan
%recovery is the only check on it.
%
%Inputs (name-value; all posterior-draw inputs are draws x N matrices)
% b           - per-event log-mean intercepts at z = 0 (draws x 2). z = 0 is the
%               OVERALL-sample mean of each dimension (standardized once across
%               the whole sample; see psyrat_zscore_subjectlevel).
% b_sigma     - per-event log residual-SD intercepts at z = 0 (draws x 2).
% b_dim       - event-crossed log-mean dimension slopes (draws x KDIM; KDIM = 2
%               for one dimension, 6 for two). Odd columns event 1, even event 2.
% b_sigma_dim - event-crossed log residual-SD dimension slopes (draws x KDIM),
%               same layout.
% sd_id       - person log-scale SDs (draws x 4):
%               [mean_1, mean_2, logsigma_1, logsigma_2].
% sd_trl      - per-event trial main-effect log-mean SDs (draws x 2).
% cor_id      - person correlation matrix draws (draws x 4 x 4). Supplies the
%               cross-event person-mean correlation cor_id(:,1,2). The
%               mean<->log-sigma blocks are NOT consumed here (see above).
% cor_i       - cross-event trial correlation on the log scale (draws x 1).
% ndim        - number of dimensions: 1 or 2.
% z1          - grid axis for dimension 1.
% obs         - n' per event. Scalar (applied to both events) or [n1 n2].
%               Difference designs are routinely unbalanced, so prefer the
%               two-element form.
% CI          - credible-interval width in (0,1].
%
%Optional Inputs
% z2          - grid axis for dimension 2 (required if ndim = 2).
%
%Output (struct)
% out.ndim, out.obs, out.ci
% out.z1 (, out.z2)   - grid axes
% out.G, out.D        - difference-score surfaces (ll/pt/ul); ndim=1 -> 1 x nz1;
%                       ndim=2 -> nz1 x nz2 (rows index z1, columns index z2).
% out.ICCg, out.ICCd  - relative/absolute single-observation difference ICC.
% out.sige_pt         - posterior-mean residual CONTRAST SD
%                       sqrt(sigma^2_pi,e,Delta) at each grid point.
%
% See also PSYRAT_REL_DIFFDYNREL_GAMMA, PSYRAT_GAMMA_VARCOMPS_LS,
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
        'See help psyrat_rel_diffdynrel_gamma_ls for more information about inputs'));
end

b           = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_i       = psyrat_req(varargin,'cor_i','Please input cor_i (cross-event trial correlation).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1.');
obs         = psyrat_req(varargin,'obs','Please input obs (n'').');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');
z2          = psyrat_opt(varargin,'z2',[]);
%CONCURRENT (analysis 19) only. Empty is the NON-concurrent default and leaves
%every existing caller byte-identical: cases 12 and 13 never pass it. When
%supplied it must carry one copula correlation per draw, on the SAME retained
%draw set as every other argument - the correlation exists only for the thinned
%stage-1 draws, so the caller subsamples before calling rather than recycling.
rho         = psyrat_opt(varargin,'rho',[]);
if ~isempty(rho)
    rho = rho(:);
    if numel(rho) ~= size(b,1)
        error('varargin:rho',...
            ['WARNING: rho has %d entries but the marginal draws have %d rows. '...
            'Subsample every draw array to REL.out.copula.sel before calling.\n'],...
            numel(rho),size(b,1));
    end
    %Same rejection as the two-facet twin: an out-of-range or non-finite rho
    %would otherwise die mid-grid inside the quadrature kernel with a less
    %actionable identifier.
    if any(~isfinite(rho)) || any(abs(rho) >= 1)
        error('varargin:rho',...
            'WARNING: every copula correlation must be finite with magnitude below one.\n');
    end
end

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
if size(b_sigma_dim,2) ~= expectedk
    error('varargin:bsigmadim',...
        'WARNING: b_sigma_dim has %d columns but ndim=%d expects %d.\n',...
        size(b_sigma_dim,2),ndim,expectedk);
end
if size(sd_id,2) ~= 4
    error('varargin:sdid',...
        ['WARNING: sd_id must have 4 columns '...
        '[mean_1, mean_2, logsigma_1, logsigma_2]; got %d.\n'],size(sd_id,2));
end

%event-crossed slopes: odd columns -> event 1, even columns -> event 2, matching
%the case-12 Xdim construction [meas1.*z1, meas2.*z1, meas1.*z2, meas2.*z2, ...].
bdim1 = b_dim(:,1:2:end);      %draws x (KDIM/2)
bdim2 = b_dim(:,2:2:end);
bsd1  = b_sigma_dim(:,1:2:end);
bsd2  = b_sigma_dim(:,2:2:end);

%per-event person mean SDs, trial SDs, and person log-sigma SDs (all log scale,
%all z-independent).
sp_a = sd_id(:,1); sp_b = sd_id(:,2);
ssd_a = sd_id(:,3); ssd_b = sd_id(:,4);
si_a = sd_trl(:,1); si_b = sd_trl(:,2);

%The CROSS-event person-mean correlation drives the difference numerator. Unlike
%the log-nu sibling there is no within-event mean<->scale covariance to form:
%psyrat_gamma_varcomps_ls has no cov_pv argument because the residual never
%touches mu. Every cor_id entry involving slots 3-4 is therefore deliberately
%unread here - (1,3), (2,4), (1,4), (2,3) and (3,4) alike.
cor_p = cor_id(:,1,2);
cor_i = cor_i(:);

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

        %z-conditioned per-event log-mean intercepts and log residual SDs. Each
        %event is evaluated at its OWN z-specific mean and scale predictors
        %(ruling 13); asymmetric slopes on either submodel legitimately move the
        %coefficient and must not be suppressed (ruling 15).
        alpha_a = b(:,1) + bdim1 * zt;
        alpha_b = b(:,2) + bdim2 * zt;
        lsig_a  = b_sigma(:,1) + bsd1 * zt;
        lsig_b  = b_sigma(:,2) + bsd2 * zt;

        %OBSERVED-score-scale components at this z and the cross-event
        %observed-scale covariances. NOTE the asymmetry between the two calls:
        %the converter is the LOCATION-SCALE one, while the cross-event helper is
        %the SAME one the log-nu path uses, because its whole argument list is
        %mean-submodel (see the header).
        vc_a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, si_a, lsig_a, ssd_a);
        vc_b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, si_b, lsig_b, ssd_b);
        cc   = psyrat_gamma_crosscov(alpha_a, sp_a, si_a, ...
                                     alpha_b, sp_b, si_b, cor_p, cor_i);

        %assemble the 2x2 person and trial (co)variance arrays psyrat_diffrel
        %expects, exactly as psyrat_gamma_extract_diff_ls does for the static case.
        id_varcov(:,1,1) = vc_a.sigma_p2;  id_varcov(:,2,2) = vc_b.sigma_p2;
        id_varcov(:,1,2) = cc.cov_p_obs;   id_varcov(:,2,1) = cc.cov_p_obs;
        trl_varcov(:,1,1) = vc_a.sigma_i2; trl_varcov(:,2,2) = vc_b.sigma_i2;
        trl_varcov(:,1,2) = cc.cov_i_obs;  trl_varcov(:,2,1) = cc.cov_i_obs;

        %psyrat_diffrel reads the residual as a log-scale SD and squares it
        %(exp(er_var).^2 = sigma_pi_e2), so encode as 0.5*log(sigma_pi_e2) -
        %the identical encoding used by psyrat_gamma_extract_diff_ls.
        er_var_z = [0.5*log(vc_a.sigma_pi_e2(:)), 0.5*log(vc_b.sigma_pi_e2(:))];

        %Residual cross-covariance is 0 (non-concurrent, ruling 16): the two
        %events come from different physical trials, so the conditional draws are
        %independent given mu and sigma. This matches psyrat_gamma_crosscov's
        %cov_e_obs and the static case-7 gamma path.
        %
        %wp_cov = 0 also OMITS the model-implied nonlinear interaction cross
        %covariance EY_A*EY_B*expm1(c_p)*expm1(c_i), which the owner's reference
        %retains. Owner decision 2026-08-05: documented as a measured
        %approximation, not fixed here - measured at 4e-6 to 6e-5 of the residual
        %variance beside it. See gamma_scale_submodel_decision.md section I.4.
        %
        %CONCURRENT OVERRIDE (analysis 19). When a copula correlation is supplied
        %the two events ARE measured on the same trials, so neither omission
        %above applies: the residual covariance exists and the interaction cross
        %term is pooled into it, exactly as the per-participant path does. The
        %typical person is u_p = 0 on both submodels, so the expected means that
        %set the quadrature point carry the population intercept at this z and
        %marginalize over the trial facet only.
        if isempty(rho)
            wp_cov_z = 0;
        else
            E_a = exp(alpha_a + 0.5.*si_a.^2);
            E_b = exp(alpha_b + 0.5.*si_b.^2);
            nu_a = 2 .* (E_a ./ exp(lsig_a)).^2;
            nu_b = 2 .* (E_b ./ exp(lsig_b)).^2;
            cov_e = zeros(ndraws,1);
            for q = 1:ndraws
                qq = psyrat_gamma_copula_rescov(nu_a(q), nu_b(q), rho(q));
                cov_e(q) = E_a(q) .* E_b(q) .* qq.value;
            end
            ccr = psyrat_gamma_crosscov_rescor(alpha_a, sp_a, si_a, ...
                alpha_b, sp_b, si_b, cor_p, cor_i, cov_e);
            wp_cov_z = ccr.cov_pi_e_obs;
        end

        g = psyrat_diffrel('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
            'obs',obs,'CI',ciperc,'est','gen','wp_cov',wp_cov_z);
        d = psyrat_diffrel('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
            'obs',obs,'CI',ciperc,'est','dep','wp_cov',wp_cov_z);

        G.ll(a,c) = g.ll; G.pt(a,c) = g.pt; G.ul(a,c) = g.ul;
        D.ll(a,c) = d.ll; D.pt(a,c) = d.pt; D.ul(a,c) = d.ul;
        ICCg.ll(a,c) = g.icc_ll; ICCg.pt(a,c) = g.icc_pt; ICCg.ul(a,c) = g.icc_ul;
        ICCd.ll(a,c) = d.icc_ll; ICCd.pt(a,c) = d.icc_pt; ICCd.ul(a,c) = d.icc_ul;

        if isfield(g,'admissibility')
            adm_all = [adm_all, g.admissibility, d.admissibility]; %#ok<AGROW>
        end

        %Residual CONTRAST SD at this z: sqrt(sigma_A^2 + sigma_B^2 - 2*cov_AB).
        %wp_cov_z IS the cross term and is exactly 0 on the non-concurrent path,
        %so cases 12 and 13 are byte-identical to the two-term form this replaced.
        %Under analysis 19 it is not 0, and omitting it overstated the reported
        %contrast SD - the events are measured on the same trials, so their
        %residuals covary. max(.,0) guards float noise only; psyrat_diffrel
        %already warns through psyrat:negerrvar on a genuinely negative contrast.
        sige_pt(a,c) = mean(sqrt(max(...
            vc_a.sigma_pi_e2 + vc_b.sigma_pi_e2 - 2.*wp_cov_z, 0)));
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
