function out = psyrat_rel_diffdynrel_trt_gamma_ls(varargin)
%Dynamic/conditional reliability of a two-event DIFFERENCE score in the TWO-FACET
%(trial x occasion) crossed design, for the LOCATION-SCALE scaled-chi-square /
%Gamma(log link) family, as a function of the standardized dimension
%predictor(s). This is the typical-person reference SURFACE of analysis 20
%(ic_diff_dynrel_sserrvar_trt_rescor, family = 'gamma', gammascale = 2), the
%two-facet arm of the MODULAR CUT-COPULA workflow.
%
% out = psyrat_rel_diffdynrel_trt_gamma_ls('b',b,'b_sigma',b_sigma,...
%   'b_dim',b_dim,'b_sigma_dim',b_sigma_dim,'sd_id',sd_id,'sd_trl',sd_trl,...
%   'sd_occ',sd_occ,'sd_tid',sd_tid,'sd_oid',sd_oid,'sd_to',sd_to,...
%   'cor_id',cor_id,'cor_t',cor_t,'cor_o',cor_o,'cor_pt',cor_pt,...
%   'cor_po',cor_po,'cor_ot',cor_ot,'ndim',ndim,'z1',z1axis,'z2',z2axis,...
%   'obs',ntrl,'nocc',nocc,'reltype',reltype,'CI',.95,'rho',rho)
%
%WHAT THIS IS. The TWO-FACET twin of PSYRAT_REL_DIFFDYNREL_GAMMA_LS. Read that
%file first: the event-crossed slope layout, the z conditioning, the
%typical-person estimand and the optional copula channel are all taken from it
%unchanged. What differs is the design, and the design changes the residual.
%
%   log mu_k(z)    = b_k       + b_dim_k       . X(z)
%   log sigma_k(z) = b_sigma_k + b_sigma_dim_k . X(z)
%   nu_k           = 2*(mu_k/sigma_k)^2   =>  Var(Y | mu, sigma) = sigma_k^2
%
%with the log-mean surface carrying six crossed effects per event (person,
%occasion, trial, and the three two-way interactions) rather than two. X(z) is
%[z1] for one dimension and [z1 z2 z1*z2] for two; ODD slope columns are event 1,
%EVEN columns event 2, the case-12/15 convention.
%
%THE RESIDUAL IS NOT THE ONE-FACET RESIDUAL, AND THE DIFFERENCE IS NOT COSMETIC.
%In the one-facet design person x trial and measurement error are confounded, so
%the toolbox reports them in a single slot and PSYRAT_DIFFREL divides that slot by
%n_i. Here person x trial is separately estimable, is reported as its own
%component, and PSYRAT_DIFFREL_TRT divides it by n_i while dividing the residual
%by n_i*n_o. The residual slot is therefore the HIGHEST-ORDER term only: the
%three-way mean-surface remainder plus the observation dispersion, which is
%exactly PSYRAT_GAMMA_VARCOMPS_TRT_LS's sigma_res2. Folding person x trial into it
%to "match" the one-facet would change the D-study arithmetic. The reference
%implementation splits it the same way (R/reliability.R:241-246 selects
%"pio_e"/"pio_mean" for the two-facet branch where the one-facet branch selects
%"pi_e"/"pi_mean"). See PSYRAT_GAMMA_CROSSCOV_TRT_RESCOR for the covariance half
%of the same statement.
%
%THE COPULA CHANNEL (rho), and why it is optional. Analysis 20's events are
%CONCURRENT - measured on the same trials - so their residuals covary. That
%covariance has no closed form under Gamma margins, which is why the concurrent
%gamma designs were barred; it is supplied here by converged Gauss-Hermite
%quadrature over the copula correlation (PSYRAT_GAMMA_COPULA_RESCOV), rescaled by
%the two events' expected means, and reaches PSYRAT_DIFFREL_TRT through er_cov.
%Omitting rho gives er_cov = 0, the non-concurrent reduction. That path exists
%for the tests, which use it to reduce this design to the one-facet one; the
%non-concurrent two-facet gamma dynamic difference design (a gamma analysis 15)
%is not enabled, so analysis 20 is the only production caller and it always
%supplies rho.
%
%THE QUADRATURE IS EVALUATED AT THE TYPICAL PERSON, and marginalizes over ALL
%FIVE non-person facets - not the trial facet alone as in the one-facet sibling.
%The typical person is u_p = 0 on both submodels, so
%E_k = exp(alpha_k(z) + 0.5*(s_t^2 + s_o^2 + s_pt^2 + s_po^2 + s_ot^2)) and
%nu_k = 2*(E_k/sigma_k(z))^2. Compare the per-participant calculator, which
%replaces alpha_k(z) with that participant's own log mean and keeps the same five
%facet variances (R/reliability.R:250-256).
%
%STRICT QUADRATURE HERE, NON-STRICT IN THE PER-PARTICIPANT CALCULATOR, and that
%asymmetry is deliberate. This surface evaluates at POPULATION nu(z), which stay
%in a well-behaved regime; the calculator evaluates at each participant's own
%nu, where a single tail draw can put the two events in wildly different regimes
%and defeat every node pair. A convergence failure here means something is wrong
%with the fit as a whole and should stop the surface.
%
%THE ESTIMAND OF THE RESIDUAL COVARIANCE. Evaluated at expected means, it is a
%MEAN-MATCHED DIAGNOSTIC quantity. It is NOT integrated over every realized
%location-facet effect and must not be described as a fully facet-marginalized
%residual covariance. That labeling requirement comes from the reference bundle
%and is carried through every output surface; see
%documentation/modular_cut_copula_workflow.md.
%
%DEVELOPMENT-ONLY. The two-facet configuration has no empirical anchor. The
%reference bundle's own two-facet runs were prepared against a SIMULATED second
%occasion and were never sampled, so unlike the one-facet design there is no
%completed external fit to replay against. The formulas here are pinned against
%the reference's R implementation via frozen fixtures, not against results.
%
%ESTIMAND OF THE OBSERVATION VARIANCE. Marginalized over the person population's
%residual-SD heterogeneity, per event, via the s_sd argument of
%PSYRAT_GAMMA_VARCOMPS_TRT_LS (E[sigma_p^2] = exp(2*log_sigma(z) + 2*s_sd^2)).
%NOTE the quadrature above deliberately uses the TYPICAL person's exp(lsig)
%rather than that marginalized value: the copula correlation describes a pair of
%observations from one person, so the point at which it is converted is that
%person's own scale, not a population average over scales.
%
%Inputs (name-value; all posterior-draw inputs are draws x N matrices, and every
%SD is on the LOG expected-score scale)
% b, b_sigma        - per-event log-mean / log residual-SD intercepts at z = 0
%                     (draws x 2).
% b_dim,
% b_sigma_dim       - event-crossed dimension slopes (draws x KDIM; KDIM = 2 for
%                     one dimension, 6 for two).
% sd_id             - person log-scale SDs (draws x 4):
%                     [mean_1, mean_2, logsigma_1, logsigma_2].
% sd_trl, sd_occ    - per-event trial / occasion main-effect SDs (draws x 2).
% sd_tid, sd_oid,
% sd_to             - per-event trial x person, occasion x person and
%                     trial x occasion SDs (draws x 2).
% cor_id            - person correlation draws (draws x 4 x 4). Supplies the
%                     cross-event person-mean correlation cor_id(:,1,2) only; the
%                     mean<->log-sigma blocks are not consumed (see the one-facet
%                     sibling's header for why adding them would be wrong).
% cor_t, cor_o,
% cor_pt, cor_po,
% cor_ot            - cross-event correlations of the five non-person facets
%                     (draws x 1 each), on the log scale.
% ndim              - number of dimensions: 1 or 2.
% z1                - grid axis for dimension 1.
% obs               - trial n' (positive scalar; applied to both events).
% nocc              - occasion n' (positive scalar).
% CI                - credible-interval width in (0,1].
%
%Optional Inputs
% reltype           - 1 (CE), 2 (CS) or 3 (CES); default 3. See PSYRAT_DIFFREL_TRT.
% z2                - grid axis for dimension 2 (required if ndim = 2).
% rho               - [draws x 1] copula correlation, one per retained draw.
%                     Empty (default) is the non-concurrent reduction. When
%                     supplied it must be on the SAME retained draw set as every
%                     other argument; the caller subsamples by
%                     REL.out.copula.sel before calling rather than recycling.
%
%Output (struct)
% out.ndim, out.obs, out.nocc, out.reltype, out.ci
% out.z1 (, out.z2)   - grid axes
% out.G, out.D        - difference-score surfaces (ll/pt/ul); ndim=1 -> 1 x nz1;
%                       ndim=2 -> nz1 x nz2 (rows index z1, columns index z2).
% out.ICCg, out.ICCd  - relative/absolute single-observation difference ICC.
% out.sige_pt         - posterior-mean residual CONTRAST SD at each grid point,
%                       sqrt(sigma_res,A^2 + sigma_res,B^2 - 2*cov_res_e,AB).
%
% See also PSYRAT_REL_DIFFDYNREL_GAMMA_LS, PSYRAT_GAMMA_VARCOMPS_TRT_LS,
% PSYRAT_GAMMA_CROSSCOV_TRT_RESCOR, PSYRAT_GAMMA_COPULA_RESCOV,
% PSYRAT_DIFFREL_TRT.

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
        'See help psyrat_rel_diffdynrel_trt_gamma_ls for more information about inputs'));
end

b           = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
sd_occ      = psyrat_req(varargin,'sd_occ','Please input sd_occ (per-event occasion log-SDs).');
sd_tid      = psyrat_req(varargin,'sd_tid','Please input sd_tid (per-event trial x person log-SDs).');
sd_oid      = psyrat_req(varargin,'sd_oid','Please input sd_oid (per-event occasion x person log-SDs).');
sd_to       = psyrat_req(varargin,'sd_to','Please input sd_to (per-event trial x occasion log-SDs).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_t       = psyrat_req(varargin,'cor_t','Please input cor_t (cross-event trial correlation).');
cor_o       = psyrat_req(varargin,'cor_o','Please input cor_o (cross-event occasion correlation).');
cor_pt      = psyrat_req(varargin,'cor_pt','Please input cor_pt (cross-event trial x person correlation).');
cor_po      = psyrat_req(varargin,'cor_po','Please input cor_po (cross-event occasion x person correlation).');
cor_ot      = psyrat_req(varargin,'cor_ot','Please input cor_ot (cross-event trial x occasion correlation).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1.');
obs         = psyrat_req(varargin,'obs','Please input obs (trial n'').');
nocc        = psyrat_req(varargin,'nocc','Please input nocc (occasion n'').');
reltype     = psyrat_opt(varargin,'reltype',3);
ciperc      = psyrat_req(varargin,'CI','Please input CI.');
z2          = psyrat_opt(varargin,'z2',[]);
%CONCURRENT (analysis 20) channel. Empty is the non-concurrent reduction; see the
%header for why no shipped analysis takes that path.
rho         = psyrat_opt(varargin,'rho',[]);
if ~isempty(rho)
    rho = rho(:);
    if numel(rho) ~= size(b,1)
        error('varargin:rho',...
            ['WARNING: rho has %d entries but the marginal draws have %d rows. '...
            'Subsample every draw array to REL.out.copula.sel before calling.\n'],...
            numel(rho),size(b,1));
    end
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
%~isfinite is tested explicitly for the same reason the one-facet sibling gives:
%NaN < 1 is FALSE in MATLAB, so a NaN n' would pass a bare magnitude check and
%then divide every error term by NaN, giving an all-NaN surface with no
%diagnostic.
if ~isscalar(obs) || ~isfinite(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a finite positive scalar.\n');
end
if ~isscalar(nocc) || ~isfinite(nocc) || nocc < 1
    error('varargin:nocc','WARNING: nocc must be a finite positive scalar.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1, 2, or 3.\n');
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
%Every crossed facet is per-event. A one-column matrix here would broadcast
%silently and give both events the same facet SD, which is a different model.
facetnames = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
facetvals  = {sd_trl,sd_occ,sd_tid,sd_oid,sd_to};
for f = 1:numel(facetnames)
    if size(facetvals{f},2) ~= 2
        error('varargin:facetsd',...
            ['WARNING: %s must have 2 columns, one per difference-score event; '...
            'got %d.\n'],facetnames{f},size(facetvals{f},2));
    end
end
if size(b,2) ~= 2 || size(b_sigma,2) ~= 2
    error('varargin:events',...
        ['WARNING: b and b_sigma must each have 2 columns, one per '...
        'difference-score event.\n']);
end

%event-crossed slopes: odd columns -> event 1, even columns -> event 2.
bdim1 = b_dim(:,1:2:end);      %draws x (KDIM/2)
bdim2 = b_dim(:,2:2:end);
bsd1  = b_sigma_dim(:,1:2:end);
bsd2  = b_sigma_dim(:,2:2:end);

%per-event log-scale SDs, all z-independent. Slot order follows
%psyrat_gamma_varcomps_trt_ls: p person, o occasion, t trial, po occasion x
%person, pt trial x person, ot trial x occasion. NOTE the toolbox's tid is
%"trial x person" and so supplies s_pt, while oid is "occasion x person" and
%supplies s_po - the names CROSS OVER.
%
%That cross-over has to be held consistent at FIVE places in this file, four of
%them positional and so unnamed at the call site: this unpacking; the two
%psyrat_gamma_varcomps_trt_ls calls (slots 5 and 6); the
%psyrat_gamma_crosscov_trt_rescor call (slots 5/6 for each event and its two
%interaction correlations); the varcov assembly below (sigma_pt2 -> tid,
%sigma_po2 -> oid); and the psyrat_diffrel_trt call ('bpi' <- tid, 'bpo' <- oid).
%A swap at any one of them is SILENT, because both indices are valid and the
%result is still a number. Change one, check all five.
sp_a  = sd_id(:,1);  sp_b  = sd_id(:,2);
ssd_a = sd_id(:,3);  ssd_b = sd_id(:,4);
st_a  = sd_trl(:,1); st_b  = sd_trl(:,2);
so_a  = sd_occ(:,1); so_b  = sd_occ(:,2);
spt_a = sd_tid(:,1); spt_b = sd_tid(:,2);
spo_a = sd_oid(:,1); spo_b = sd_oid(:,2);
sot_a = sd_to(:,1);  sot_b = sd_to(:,2);

%cross-event correlations. The person one comes out of the 4x4 block; the five
%facet ones arrive as their own draws.
cor_p  = cor_id(:,1,2);
cor_t  = cor_t(:);
cor_o  = cor_o(:);
cor_pt = cor_pt(:);
cor_po = cor_po(:);
cor_ot = cor_ot(:);

%Log-scale variance of everything EXCEPT the person effect, per event. This is
%what the typical person's expected observed score marginalizes over, and it is
%formed once because no facet SD carries a dimension slope.
v_other_a = st_a.^2 + so_a.^2 + spt_a.^2 + spo_a.^2 + sot_a.^2;
v_other_b = st_b.^2 + so_b.^2 + spt_b.^2 + spo_b.^2 + sot_b.^2;

z1 = z1(:)'; nz1 = numel(z1);
if ndim == 2; z2 = z2(:)'; nz2 = numel(z2); else; nz2 = 1; end

gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);  D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);  ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);

ndraws = size(b,1);
id_varcov  = zeros(ndraws,2,2);
trl_varcov = zeros(ndraws,2,2);
occ_varcov = zeros(ndraws,2,2);
tid_varcov = zeros(ndraws,2,2);
oid_varcov = zeros(ndraws,2,2);
to_varcov  = zeros(ndraws,2,2);

%One warning per surface, not one per grid point; see the one-facet sibling.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for a = 1:nz1
    for c = 1:nz2
        if ndim == 1
            zterms = z1(a);
        else
            zterms = [z1(a); z2(c); z1(a)*z2(c)];
        end
        zt = zterms(:);

        %z-conditioned per-event log-mean intercepts and log residual SDs. Each
        %event is evaluated at its OWN z-specific predictors; asymmetric slopes
        %on either submodel legitimately move the coefficient.
        alpha_a = b(:,1) + bdim1 * zt;
        alpha_b = b(:,2) + bdim2 * zt;
        lsig_a  = b_sigma(:,1) + bsd1 * zt;
        lsig_b  = b_sigma(:,2) + bsd2 * zt;

        %CONCURRENT residual covariance, at the typical person's expected means.
        %Strict quadrature - see the header for why this differs from the
        %per-participant calculator.
        if isempty(rho)
            cov_e = 0;
        else
            E_a = exp(alpha_a + 0.5.*v_other_a);
            E_b = exp(alpha_b + 0.5.*v_other_b);
            nu_a = 2 .* (E_a ./ exp(lsig_a)).^2;
            nu_b = 2 .* (E_b ./ exp(lsig_b)).^2;
            cov_e = zeros(ndraws,1);
            for q = 1:ndraws
                qq = psyrat_gamma_copula_rescov(nu_a(q), nu_b(q), rho(q));
                cov_e(q) = E_a(q) .* E_b(q) .* qq.value;
            end
        end

        %OBSERVED-score-scale components at this z, and the cross-event
        %observed-scale covariances. The converter is the LOCATION-SCALE one;
        %the cross-event helper wraps the shipped two-facet one and adds only
        %the concurrent residual term.
        vc_a = psyrat_gamma_varcomps_trt_ls(alpha_a, sp_a, so_a, st_a, ...
            spo_a, spt_a, sot_a, lsig_a, ssd_a);
        vc_b = psyrat_gamma_varcomps_trt_ls(alpha_b, sp_b, so_b, st_b, ...
            spo_b, spt_b, sot_b, lsig_b, ssd_b);
        cc = psyrat_gamma_crosscov_trt_rescor( ...
            alpha_a, sp_a, so_a, st_a, spo_a, spt_a, sot_a, ...
            alpha_b, sp_b, so_b, st_b, spo_b, spt_b, sot_b, ...
            cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot, cov_e);

        %assemble the six 2x2 (co)variance arrays psyrat_diffrel_trt expects
        %(id->bp, tid->bpi, oid->bpo, trl->bt, occ->bo, to->boi).
        id_varcov(:,1,1)  = vc_a.sigma_p2;   id_varcov(:,2,2)  = vc_b.sigma_p2;
        id_varcov(:,1,2)  = cc.cov_p;        id_varcov(:,2,1)  = cc.cov_p;
        trl_varcov(:,1,1) = vc_a.sigma_t2;   trl_varcov(:,2,2) = vc_b.sigma_t2;
        trl_varcov(:,1,2) = cc.cov_t;        trl_varcov(:,2,1) = cc.cov_t;
        occ_varcov(:,1,1) = vc_a.sigma_o2;   occ_varcov(:,2,2) = vc_b.sigma_o2;
        occ_varcov(:,1,2) = cc.cov_o;        occ_varcov(:,2,1) = cc.cov_o;
        tid_varcov(:,1,1) = vc_a.sigma_pt2;  tid_varcov(:,2,2) = vc_b.sigma_pt2;
        tid_varcov(:,1,2) = cc.cov_pt;       tid_varcov(:,2,1) = cc.cov_pt;
        oid_varcov(:,1,1) = vc_a.sigma_po2;  oid_varcov(:,2,2) = vc_b.sigma_po2;
        oid_varcov(:,1,2) = cc.cov_po;       oid_varcov(:,2,1) = cc.cov_po;
        to_varcov(:,1,1)  = vc_a.sigma_ot2;  to_varcov(:,2,2)  = vc_b.sigma_ot2;
        to_varcov(:,1,2)  = cc.cov_ot;       to_varcov(:,2,1)  = cc.cov_ot;

        %psyrat_diffrel_trt reads the residual as a log-scale SD and squares it
        %(exp(er_var).^2 = sigma_res2), so encode as 0.5*log(sigma_res2). The
        %covariance channel is read directly on the variance scale, hence no
        %transform there.
        er_var_z = [0.5*log(vc_a.sigma_res2(:)), 0.5*log(vc_b.sigma_res2(:))];
        er_cov_z = cc.cov_res_e_obs;

        g = psyrat_diffrel_trt('bp',id_varcov,'bpi',tid_varcov,'bpo',oid_varcov,...
            'bt',trl_varcov,'bo',occ_varcov,'boi',to_varcov,...
            'er_var',er_var_z,'obs',[obs obs],'nocc',[nocc nocc],...
            'reltype',reltype,'est','gen','er_cov',er_cov_z,'CI',ciperc);
        d = psyrat_diffrel_trt('bp',id_varcov,'bpi',tid_varcov,'bpo',oid_varcov,...
            'bt',trl_varcov,'bo',occ_varcov,'boi',to_varcov,...
            'er_var',er_var_z,'obs',[obs obs],'nocc',[nocc nocc],...
            'reltype',reltype,'est','dep','er_cov',er_cov_z,'CI',ciperc);

        G.ll(a,c) = g.ll; G.pt(a,c) = g.pt; G.ul(a,c) = g.ul;
        D.ll(a,c) = d.ll; D.pt(a,c) = d.pt; D.ul(a,c) = d.ul;
        ICCg.ll(a,c) = g.icc_ll; ICCg.pt(a,c) = g.icc_pt; ICCg.ul(a,c) = g.icc_ul;
        ICCd.ll(a,c) = d.icc_ll; ICCd.pt(a,c) = d.icc_pt; ICCd.ul(a,c) = d.icc_ul;

        if isfield(g,'admissibility')
            adm_all = [adm_all, g.admissibility, d.admissibility]; %#ok<AGROW>
        end

        %Residual CONTRAST SD at this z. er_cov_z is the pooled cross term and is
        %exactly the mean-surface three-way covariance when rho is absent, so the
        %non-concurrent reduction still subtracts the right thing. max(.,0)
        %guards float noise only.
        sige_pt(a,c) = mean(sqrt(max(...
            vc_a.sigma_res2 + vc_b.sigma_res2 - 2.*er_cov_z, 0)));
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

out = struct();
out.ndim = ndim; out.obs = obs; out.nocc = nocc; out.reltype = reltype;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2; out.z2 = z2; end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.G = gridops.squeeze(G,ndim);
out.D = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
