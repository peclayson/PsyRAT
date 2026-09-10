function [gen,dep,diagnostic] = psyrat_ssrel_diffdynrel_rescor_gamma_ls(varargin)
%Per-participant CONDITIONAL dynamic/conditional CONCURRENT DIFFERENCE-score
%reliability for the one-facet gamma location-scale MODULAR CUT-COPULA design
%(analysis 19, ic_diff_dynrel_sserrvar_rescor, family = 'gamma',
%gammascale = 2). Stage 3 of the modular workflow.
%
%[gen,dep,diagnostic] = psyrat_ssrel_diffdynrel_rescor_gamma_ls('b',b,...)
%
%WHAT THIS IS. The concurrent sibling of PSYRAT_SSREL_DIFFDYNREL_GAMMA_LS
%(analysis 13). Read that file first: the per-participant loop, the z
%conditioning, the one-row psyrat_ssrel_diff call and the id2 remap are all
%taken from it unchanged. THREE things differ, and each is a scientific choice
%rather than plumbing:
%
% (1) THE RESIDUALS ARE COUPLED. The two events are measured on the SAME trials,
%     so their residuals covary. That covariance has no closed form under Gamma
%     margins - which is why the concurrent designs were barred - and is
%     supplied here by converged Gauss-Hermite quadrature over the copula
%     correlation (PSYRAT_GAMMA_COPULA_RESCOV), rescaled by the two components'
%     expected means. It reaches psyrat_ssrel_diff through the wp_cov_ss channel
%     that analysis 13 fills with zeros.
%
% (2) THE RESIDUAL CHANNEL IS POOLED: the mean-surface person x trial term
%     PLUS each participant's conditional variance, which is what Rast &
%     Clayson's residual is (page 4: the person-by-task interaction together
%     with unsystematic measurement error; page 7 for the condition-specific
%     form - cite by PAGE, the display is unnumbered). Under a Gaussian model
%     there is nothing to pool, because with no log link the interaction IS the
%     residual; the log link is what splits it out. The pooled quantity is
%     exactly psyrat_gamma_varcomps_ls's sigma_pi_e2 when the participant's own
%     log residual SD is passed as log_sigma with s_sd = 0.
%     Analysis 13 POOLS THE SAME WAY since S17-FIX (owner-directed 2026-08-09;
%     before that it passed the conditional variance alone - finding S17, now
%     resolved). This file's construction is the one 13 adopted, so the two
%     designs' residual channels now agree by design.
%
% (3) THE PARTICIPANT'S OWN LOCATION EFFECT ENTERS - THROUGH THE ERROR, NOT THE
%     SIGNAL. The signal stays a POPULATION variance evaluated at the
%     participant's z, exactly as in analysis 13. But the degrees of freedom fed
%     to the quadrature are computed at that participant's OWN expected
%     component means, which carry their realized location effect (mu_ss). So
%     two participants at the same z with different amplitudes get different
%     residual covariances. This is the reference's estimand and is why mu_ss1/2
%     are stored at estimation time; without them the quadrature point cannot be
%     reconstructed.
%
%THE ESTIMAND OF THE RESIDUAL COVARIANCE, and what it is not. The quadrature is
%evaluated at each participant's expected means, making it a MEAN-MATCHED
%DIAGNOSTIC quantity. It is NOT integrated over every realized location-facet
%effect and must not be described as a fully facet-marginalized residual
%covariance. That labeling requirement comes from the reference bundle and is
%carried through every output surface; see
%documentation/modular_cut_copula_workflow.md.
%
%DRAW ALIGNMENT IS THE CALLER'S JOB. The copula correlation exists only for the
%retained (thinned) stage-1 draws, so every array passed here must already be
%subsampled to that same set and share one draw dimension. Passing full-length
%marginal draws with a short rho is rejected rather than recycled: silently
%pairing draw 1's rho with draw 5000's margins would produce a plausible number
%from an incoherent posterior.
%
%INPUTS (name-value; every SD is on the LOG expected-score scale)
%   b, b_sigma, b_dim, b_sigma_dim, sd_id, sd_trl, cor_id, cor_i,
%   er_ss1, er_ss2, idtable, ndim, CI
%                - as PSYRAT_SSREL_DIFFDYNREL_GAMMA_LS; see that file
%   mu_ss1/mu_ss2- [ndraws x NSUB] per-subject ABSOLUTE log-mean at z = 0
%                  (b[e] + r_id[:,e]); the dimension contribution is added here
%   rho          - [ndraws x 1] copula correlation, one per retained draw
%
%OUTPUTS
%   gen, dep   - stacked one-row-per-participant tables, as psyrat_ssrel_diff
%                returns for est = 'gen' and est = 'dep'
%   diagnostic - table of per-participant quadrature provenance: the worst node
%                count and error over draws, and the posterior-mean residual
%                covariance and correlation. Surfaced so a reader can see where
%                the quadrature had to work hardest.
%
%See also PSYRAT_SSREL_DIFFDYNREL_GAMMA_LS, PSYRAT_GAMMA_COPULA_RESCOV,
%PSYRAT_GAMMA_CROSSCOV_RESCOR, PSYRAT_GAMMA_VARCOMPS_LS, PSYRAT_SSREL_DIFF.

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

b           = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_i       = psyrat_req(varargin,'cor_i','Please input cor_i (cross-event trial correlation).');
er_ss1      = psyrat_req(varargin,'er_ss1','Please input er_ss1 (per-subject log residual SD at z=0, event 1).');
er_ss2      = psyrat_req(varargin,'er_ss2','Please input er_ss2 (per-subject log residual SD at z=0, event 2).');
mu_ss1      = psyrat_req(varargin,'mu_ss1','Please input mu_ss1 (per-subject log mean at z=0, event 1).');
mu_ss2      = psyrat_req(varargin,'mu_ss2','Please input mu_ss2 (per-subject log mean at z=0, event 2).');
rho         = psyrat_req(varargin,'rho','Please input rho (copula correlation per retained draw).');
idtable     = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if size(sd_id,2) ~= 4
    error('varargin:sd_id',...
        ['WARNING: sd_id must have 4 columns, ordered '...
        '[mean_1, mean_2, logsigma_1, logsigma_2]; got %d.\n'],size(sd_id,2));
end
if size(b,2) ~= 2 || size(b_sigma,2) ~= 2 || size(sd_trl,2) ~= 2
    error('varargin:events',...
        ['WARNING: b, b_sigma and sd_trl must each have 2 columns, one per '...
        'difference-score event.\n']);
end
kdim = size(b_dim,2);
if size(b_sigma_dim,2) ~= kdim || mod(kdim,2) ~= 0 || kdim ~= 2*(2*ndim-1)
    error('varargin:b_dim',...
        ['WARNING: b_dim and b_sigma_dim must both have %d columns for '...
        'ndim = %d (2 per dimension term, event-crossed); got %d and %d.\n'],...
        2*(2*ndim-1),ndim,kdim,size(b_sigma_dim,2));
end

ndraws = size(b,1);
rho = rho(:);
%DRAW ALIGNMENT. Every array must already be on the retained-draw set; see the
%header. Recycling a short rho would pair unrelated draws.
if numel(rho) ~= ndraws
    error('varargin:rho',...
        ['WARNING: rho has %d entries but the marginal draws have %d rows. '...
        'The copula correlation exists only for the RETAINED stage-1 draws, so '...
        'the caller must subsample every marginal array to that same set '...
        'before calling. How to fix: index the draw arrays by REL.out.copula '...
        '.sel.\n'],numel(rho),ndraws);
end
if any(~isfinite(rho)) || any(abs(rho) >= 1)
    error('varargin:rho',...
        'WARNING: every copula correlation must be finite with magnitude below one.\n');
end

nsub = height(idtable);
if size(er_ss1,2) < max(idtable.id2) || size(er_ss2,2) < max(idtable.id2)
    %FAIL LOUDLY, for the reason analysis 13 documents: psyrat_ssrel_diff's own
    %range check silently falls back to the POPULATION residual for an id2 it
    %cannot index, turning a subject-level table into a population one.
    error('varargin:er_ss',...
        ['WARNING: er_ss1/er_ss2 have %d and %d subject columns but idtable '...
        'references id2 up to %d.\n'],...
        size(er_ss1,2),size(er_ss2,2),max(idtable.id2));
end
if size(mu_ss1,2) < max(idtable.id2) || size(mu_ss2,2) < max(idtable.id2)
    error('varargin:mu_ss',...
        ['WARNING: mu_ss1/mu_ss2 have %d and %d subject columns but idtable '...
        'references id2 up to %d. These carry each participant''s own location '...
        'effect and set the point at which the copula covariance is '...
        'evaluated, so a short matrix would silently move the estimand.\n'],...
        size(mu_ss1,2),size(mu_ss2,2),max(idtable.id2));
end

%event-crossed slopes: odd columns -> event 1, even -> event 2.
bdim1 = b_dim(:,1:2:end);
bdim2 = b_dim(:,2:2:end);
bsd1  = b_sigma_dim(:,1:2:end);
bsd2  = b_sigma_dim(:,2:2:end);

sp_a = sd_id(:,1); sp_b = sd_id(:,2);
si_a = sd_trl(:,1); si_b = sd_trl(:,2);
cor_p = cor_id(:,1,2);
cor_i = cor_i(:);

%The only non-person facet in this design is the trial main effect, so the
%marginalization that turns a participant's log-mean into an expected observed
%score is over that alone.
v_other_a = si_a.^2;
v_other_b = si_b.^2;

id_varcov  = zeros(ndraws,2,2);
trl_varcov = zeros(ndraws,2,2);

genRows = cell(nsub,1);
depRows = cell(nsub,1);
diagRows = cell(nsub,1);

%One warning per table rather than one per participant; see analysis 13.
adm_ws = warning('off','psyrat:negerrvar');
cleanupWs = onCleanup(@() warning(adm_ws));

for c = 1:nsub

    if ndim == 1
        zt = idtable.z1(c);
    else
        zt = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end
    zt = zt(:);
    col = idtable.id2(c);

    %POPULATION log-mean and log residual SD at this participant's z. The signal
    %components are built from these, not from the participant's own effects:
    %a participant's amplitude says where they sit, not how spread the
    %population is.
    alpha_a = b(:,1) + bdim1 * zt;
    alpha_b = b(:,2) + bdim2 * zt;
    lsig_a  = b_sigma(:,1) + bsd1 * zt;
    lsig_b  = b_sigma(:,2) + bsd2 * zt;

    %THIS PARTICIPANT'S OWN log residual SD and log mean at their z.
    ss_a = er_ss1(:,col) + bsd1 * zt;
    ss_b = er_ss2(:,col) + bsd2 * zt;
    mu_a = mu_ss1(:,col) + bdim1 * zt;
    mu_b = mu_ss2(:,col) + bdim2 * zt;

    %Expected observed score for this participant, marginalizing over the trial
    %facet but NOT over their own location effect - that effect is realized, not
    %a population to average across. These set the point at which the copula
    %correlation is converted to an observed-scale covariance.
    E_a = exp(mu_a + 0.5.*v_other_a);
    E_b = exp(mu_b + 0.5.*v_other_b);
    sig_a = exp(ss_a);
    sig_b = exp(ss_b);
    nu_a = 2 .* (E_a ./ sig_a).^2;
    nu_b = 2 .* (E_b ./ sig_b).^2;

    %QUADRATURE, per draw. The converter is scalar because its node escalation
    %is per-argument: two draws can need different node counts, and accepting
    %the looser one for both would silently reduce accuracy.
    %
    %NON-STRICT, and the reason is a failure mode measured rather than imagined.
    %A single draw can put the two components in wildly different regimes - one
    %residual SD far below its mean, the other far above - where the integrand
    %is genuinely hard and no node pair meets tolerance. Aborting there would
    %destroy an entire multi-hour analysis over a handful of tail draws. The
    %best estimate is kept, each such draw is COUNTED, and the count is surfaced
    %per participant so a reader can see it. A high count is a statement about
    %the fit's margins, not about the quadrature.
    cov_e = zeros(ndraws,1);
    nodes = zeros(ndraws,1);
    qerr = zeros(ndraws,1);
    unconverged = 0;
    for d = 1:ndraws
        q = psyrat_gamma_copula_rescov(nu_a(d), nu_b(d), rho(d), [], [], false);
        cov_e(d) = E_a(d) .* E_b(d) .* q.value;
        nodes(d) = q.nodes;
        qerr(d) = q.error;
        if ~q.converged
            unconverged = unconverged + 1;
        end
    end

    %OBSERVED-scale components at this participant's z. Passing the PARTICIPANT'S
    %log residual SD as log_sigma (with s_sd = 0) makes sigma_pi_e2 the pooled
    %channel this design reports: the mean-surface person x trial term, built
    %from the population alpha at z, plus that participant's own conditional
    %variance. See the header, point (2).
    vc_a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, si_a, ss_a, 0);
    vc_b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, si_b, ss_b, 0);
    cc   = psyrat_gamma_crosscov_rescor(alpha_a, sp_a, si_a, ...
                                        alpha_b, sp_b, si_b, cor_p, cor_i, cov_e);

    id_varcov(:,1,1) = vc_a.sigma_p2;  id_varcov(:,2,2) = vc_b.sigma_p2;
    id_varcov(:,1,2) = cc.cov_p_obs;   id_varcov(:,2,1) = cc.cov_p_obs;
    trl_varcov(:,1,1) = vc_a.sigma_i2; trl_varcov(:,2,2) = vc_b.sigma_i2;
    trl_varcov(:,1,2) = cc.cov_i_obs;  trl_varcov(:,2,1) = cc.cov_i_obs;

    %psyrat_ssrel_diff reads both residual channels as LOG SDs and squares them,
    %so the pooled variance is handed over as 0.5*log(.). The covariance channel
    %is read directly on the variance scale, hence no transform there.
    er_ss_a = 0.5*log(vc_a.sigma_pi_e2);
    er_ss_b = 0.5*log(vc_b.sigma_pi_e2);
    vcp_a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, si_a, lsig_a, 0);
    vcp_b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, si_b, lsig_b, 0);
    er_var_z = [0.5*log(vcp_a.sigma_pi_e2), 0.5*log(vcp_b.sigma_pi_e2)];

    row = idtable(c,:);
    row.id2 = 1;

    g = psyrat_ssrel_diff('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
        'idtable',row,'CI',ciperc,'est','gen',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,...
        'wp_cov_ss',cc.cov_pi_e_obs,'wp_cov',0);
    d_ = psyrat_ssrel_diff('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
        'idtable',row,'CI',ciperc,'est','dep',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,...
        'wp_cov_ss',cc.cov_pi_e_obs,'wp_cov',0);

    g.id2(:) = idtable.id2(c);
    d_.id2(:) = idtable.id2(c);
    genRows{c} = g;
    depRows{c} = d_;

    %Provenance for the quadrature, per participant. res_cor is the residual
    %CORRELATION on the observed scale, which is NOT the copula correlation and
    %is reported so the two are not confused.
    res_cor = cov_e ./ (sig_a .* sig_b);
    diagRows{c} = table(idtable.id(c), idtable.id2(c), ...
        max(nodes), max(qerr), unconverged, unconverged/ndraws, ...
        mean(cov_e), mean(res_cor), mean(rho), ...
        'VariableNames',{'id','id2','quad_nodes_max','quad_err_max',...
        'quad_unconverged_n','quad_unconverged_frac',...
        'res_cov_pt','res_cor_pt','rho_e_pt'});
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});
diagnostic = vertcat(diagRows{:});

end
