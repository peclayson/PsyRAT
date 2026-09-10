function [gen,dep] = psyrat_ssrel_dynrel_trt_gamma_ls(varargin)
%Per-participant CONDITIONAL dynamic/conditional reliability for the TWO-FACET
%(trial + occasion) LOCATION-SCALE Gamma / scaled-chi-square subject-level design
%(analysis 27, ic_dynrel_sserrvar_trt, family = 'gamma', gammascale = 2).
%
%[gen,dep] = psyrat_ssrel_dynrel_trt_gamma_ls('alpha0',a0,'b',b,...
%  'logsig0',ls0,'b_sigma',bs,'sig_p',sp,'sig_occ',so,'sig_trl',st,...
%  'sig_occxid',spo,'sig_trlxid',spt,'sig_trlxocc',sot,'ind_sd',w,...
%  'idtable',ssinfo,'ndim',1,'reltype',3,'nocc',1,'CI',.95)
%
%WHAT THIS IS. The subject-level sibling of psyrat_rel_dynrel_trt_gamma_ls, and
%the location-scale twin of psyrat_ssrel_dynrel_trt_gamma. The surface function
%evaluates a POPULATION coefficient over a z grid, marginalizing over the person
%residual-SD population; this function evaluates one coefficient PER PARTICIPANT,
%at that participant's own z, CONDITIONING on that participant's own log residual
%SD.
%
%WHY IT IS SIMPLER THAN ITS LOG-NU TWIN, AND WHY THAT IS THE POINT.
%Under gammascale = 1 the conditional variance is 2*mu^2/nu, so the twin must
%DERIVE each participant's residual from two person effects (ind_bs = u_p and
%ind_nu = v_p) and carry the within-person log-mean variance
%V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2 as an exp(2*V_w) factor. Here the
%model estimates the quantity directly:
%
%   log sigma_s(z) = logsig0 + b_sigma . x(z) + w_s        (w_s = ind_sd)
%   Var(Y | s)     = sigma_s(z)^2 = exp(2*(logsig0 + b_sigma . x(z) + w_s))
%
%with no mu anywhere. **V_w HAS NO COUNTERPART HERE AND MUST NOT BE ADDED.** It
%exists in the twin only because the log-nu conditional variance is an
%expectation of 2*mu^2/nu over the within-person facets; sigma_s is a single
%per-participant quantity carrying no facet terms, so there is no expectation to
%take. This is the two-facet form of the same rule
%psyrat_gamma_varcomps_trt_ls's header states for the group path, and the
%two-facet form of the one-facet twin's dropped exp(2*s_i^2).
%
%CONSEQUENCE: ind_bs IS NOT AN INPUT, AND MUST NOT BECOME ONE.
%The participant's own log-mean u_s does not appear in this calculation at all:
%
% (a) The residual is exp(2*(log sigma_s)), which never touches mu. Under the
%     log-nu twin it is 2*mu(z,u_s)^2*exp(2*V_w)/nu, so u_s DOES enter there and
%     a higher-mean participant receives a strictly lower conditional
%     coefficient.
% (b) sigma_p^2 is a BETWEEN-person variance, so it is evaluated at alpha(z_s)
%     and NOT at alpha(z_s) + u_s - conditioning on a participant tells you where
%     they sit, not how spread the population is. This is the same rule the
%     log-nu twin states, and it is unchanged here (owner ruling 11: the signal
%     components remain population quantities evaluated at that participant's z).
%
%So under gammascale = 2 the per-participant coefficient is a FUNCTION of
%(z_s, w_s) only, with the participant's own log-mean appearing nowhere in it.
%That decoupling is the reason this parameterization exists. It is a property,
%not an omission: a future reader who adds ind_bs "for parity" would reintroduce
%exactly the mean-dispersion blending the owner ruled against.
%
%THAT IS FUNCTIONAL INDEPENDENCE, NOT STATISTICAL INDEPENDENCE - do not restate
%it as "the coefficient is independent of the participant's amplitude". u_p and
%w_p are columns 1 and 2 of ONE correlated 2-D LKJ Cholesky block in the fitted
%model (grep chol_corrmat in the gammascale == 2 branch of
%psyrat_build_stan_dynrel_trt_chisq), and rho_p,sigma is estimated and reported
%in the component table. Whenever it is nonzero, a participant's own log-mean
%shifts the POSTERIOR of their w_s, so per-participant coefficients will still
%track amplitude across the sample. What this parameterization removes is the
%amplitude term from the FORMULA; it does not and cannot remove an estimated
%correlation from the data. This is the same formula-versus-fitted-numbers
%distinction PSYRAT_AUDIT_FINDINGS.md already draws for S16 on the static
%subject-level designs, and it is the honest version of the claim.
%
%Unlike the one-facet location-scale twin, the estimator does not lift ind_bs
%into REL.out for this design - nothing reads it. Note that is a statement about
%REL.out only: the Stan model still declares, computes and samples ind_bs (it is
%in the emitted golden and in the mean submodel), it simply is not extracted.
%
%WHAT DOES STILL MOVE WITH z. Each participant sits at their own z_s, and both
%submodels are evaluated there: the six signal components carry exp(2*alpha(z_s))
%while the residual carries exp(2*(logsig0 + b_sigma.x(z_s))). Unlike the log-nu
%twin, those do NOT share a common factor, so the log-mean slope b does not cancel
%out of the coefficient here either (see psyrat_rel_dynrel_trt_gamma_ls for the
%algebra). Both slopes matter.
%
%HOW THE PER-PARTICIPANT CALL WORKS. As in both twins, psyrat_ssrel_trt is called
%ONCE PER PARTICIPANT on a one-row idtable and the rows are stacked. Under gamma
%every observed-scale component carries exp(2*alpha(z)) and each participant sits
%at a different z, so bp/bo/bt/txp/oxp/txo differ per participant and cannot be
%passed as single shared columns in one call. Nothing in psyrat_ssrel_trt changes,
%and the per-subject coefficient, ICC and SEM formulas stay in exactly one place.
%
%(s_sd) is passed to psyrat_gamma_varcomps_trt_ls as ZERO. Only the six signal
%components are consumed here and none depends on it, so the choice is
%numerically inconsequential; zero is passed because it states the conditional
%intent. THREE returned fields do depend on s_sd and none of them is read here:
%vc.e_cond_var (where s_sd acts directly), vc.total_var (= var_mu + e_cond_var),
%and vc.sigma_res2 (derived from total_var by subtraction, NOT from e_cond_var
%directly). The residual comes from the per-person plug-in above, and THAT is
%what makes this calculation conditional rather than marginal; substituting the
%group vc.sigma_res2 would silently undo the conditioning. If you ever add a read
%here, check it against that three-field list rather than against e_cond_var
%alone.
%
%THIS IS AN EXTRAPOLATION, NOT AN IMPLEMENTATION OF RAST & CLAYSON. That paper's
%location-scale model is Gaussian throughout. Cite it for the idea; do not cite
%it for this model's behavior.
%
%Inputs (all draw vectors are columns of posterior draws)
% alpha0     - log-mean grand intercept at z = 0 (Intercept)       [ndraws x 1]
% b          - log-mean dimension slopes                        [ndraws x KDIM]
% logsig0    - log residual-SD grand intercept at z = 0 (Intercept_sigma)
%                                                                 [ndraws x 1]
% b_sigma    - log residual-SD dimension slopes                 [ndraws x KDIM]
% sig_p      - person log-mean SD (s_p)                            [ndraws x 1]
% sig_occ    - occasion main-effect log-mean SD                    [ndraws x 1]
% sig_trl    - trial main-effect log-mean SD                       [ndraws x 1]
% sig_occxid - occasion x person log-mean SD (s_po)                [ndraws x 1]
% sig_trlxid - trial x person log-mean SD (s_pt)                   [ndraws x 1]
% sig_trlxocc- trial x occasion log-mean SD (s_ot)                 [ndraws x 1]
% ind_sd     - per-participant log residual-SD effects w_p      [ndraws x NSUB]
% idtable    - ssinfo: one row per participant with id, id2, z1[, z2], trls,
%              nocc. id2 indexes the ind_sd COLUMNS.
% ndim       - 1 or 2 dimensions
% reltype    - 1, 2, or 3 (see psyrat_rel_trt)
% nocc       - n' (occasion count)
% CI         - credible interval width in decimal format (.95 = 95%)
%
%Outputs
% gen - psyrat_ssrel_trt output table for generalizability (gcoeff = 2), one row
%       per participant in idtable row order
% dep - the same for dependability (gcoeff = 1)
%
%See also PSYRAT_SSREL_DYNREL_TRT_GAMMA, PSYRAT_REL_DYNREL_TRT_GAMMA_LS,
%PSYRAT_SSREL_DYNREL_GAMMA_LS, PSYRAT_GAMMA_VARCOMPS_TRT_LS, PSYRAT_SSREL_TRT.

% Copyright (C) 2016-2026 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

alpha0  = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b       = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
logsig0 = psyrat_req(varargin,'logsig0','Please input logsig0 (Intercept_sigma draws).');
b_sigma = psyrat_req(varargin,'b_sigma','Please input b_sigma (log residual-SD dimension slopes).');
sig_p   = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
s_o     = psyrat_req(varargin,'sig_occ','Please input sig_occ (occasion log-mean SD draws).');
s_t     = psyrat_req(varargin,'sig_trl','Please input sig_trl (trial log-mean SD draws).');
s_po    = psyrat_req(varargin,'sig_occxid','Please input sig_occxid (person x occasion SD draws).');
s_pt    = psyrat_req(varargin,'sig_trlxid','Please input sig_trlxid (person x trial SD draws).');
s_ot    = psyrat_req(varargin,'sig_trlxocc','Please input sig_trlxocc (trial x occasion SD draws).');
ind_sd  = psyrat_req(varargin,'ind_sd','Please input ind_sd (per-subject log residual-SD effects).');
idtable = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim    = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
reltype = psyrat_req(varargin,'reltype','Please input reltype (1, 2, or 3).');
nocc    = psyrat_req(varargin,'nocc','Please input nocc (occasion count).');
ciperc  = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end

%column vectors of draws
alpha0  = alpha0(:);
logsig0 = logsig0(:);
sig_p   = sig_p(:);
s_o     = s_o(:);
s_t     = s_t(:);
s_po    = s_po(:);
s_pt    = s_pt(:);
s_ot    = s_ot(:);

%b/b_sigma must be draws x KDIM (KDIM = 1 for one dimension; 3 for two), matching
%the Xdim layout the estimator builds. Same shape check
%psyrat_rel_dynrel_trt_gamma_ls applies, so a transposed input fails loudly rather
%than broadcasting into a wrong answer.
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if size(b,2) ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        size(b,2),ndim,expectedk);
end
if size(b_sigma,2) ~= expectedk
    error('varargin:bsigma',...
        'WARNING: b_sigma has %d columns but ndim=%d expects %d.\n',...
        size(b_sigma,2),ndim,expectedk);
end

if any(idtable.id2 > size(ind_sd,2)) || any(idtable.id2 < 1)
    error('varargin:id2',...
        ['WARNING: idtable.id2 indexes outside the per-subject draw columns '...
        '(max id2 = %d, columns = %d).\n'],max(idtable.id2),size(ind_sd,2));
end

nsub = height(idtable);
genRows = cell(nsub,1);
depRows = cell(nsub,1);

for c = 1:nsub
    col = idtable.id2(c);

    %dimension design row for THIS participant's z (KDIM x 1), matching Xdim:
    %[z1] for one dimension, [z1; z2; z1*z2] for two.
    if ndim == 1
        xrow = idtable.z1(c);
    else
        xrow = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end

    %z-conditioned log-mean intercept and population log residual SD at this
    %participant's z
    alpha_z  = alpha0  + b       * xrow(:);   % [ndraws x 1]
    logsig_z = logsig0 + b_sigma * xrow(:);   % [ndraws x 1]

    %observed-scale crossed SIGNAL components at this participant's z. s_sd = 0:
    %the conditional estimand; the six components do not depend on it anyway, and
    %vc.e_cond_var / vc.sigma_res2 (the population residual) are deliberately not
    %read - the residual below is this participant's own.
    %
    %*** logsig_z GOES IN UNEXPONENTIATED HERE AND EXPONENTIATED FIVE LINES DOWN.
    %BOTH ARE CORRECT. *** Slot 8 of psyrat_gamma_varcomps_trt_ls is log_sigma, a
    %LOG-scale quantity that the converter exponentiates internally (grep
    %e_cond_var = exp(2.*log_sigma ...)), whereas err_ss below must be an SD on the
    %OBSERVED scale because psyrat_ssrel_trt squares its 'err' input. Note this is
    %also where the LS call stops being a rename of the log-nu one: that sibling
    %passes nu_z = EXP(lognu0 + b_nu*xrow), because ITS slot 8 is a natural-scale
    %precision. Do not "harmonize" these three sites - each takes the scale its
    %callee documents, and getting one wrong yields finite, plausible, silently
    %wrong numbers rather than an error.
    vc = psyrat_gamma_varcomps_trt_ls(alpha_z, sig_p, s_o, s_t, ...
        s_po, s_pt, s_ot, logsig_z, 0);

    %per-participant CONDITIONAL residual SD at their own z and their own w_s.
    %Exponentiated (unlike the converter argument above) because psyrat_ssrel_trt
    %takes 'err' as an observed-scale SD. No mu term and no exp(2*V_w)
    %within-person factor - see the header.
    err_ss = exp(logsig_z + ind_sd(:,col));   % [ndraws x 1]

    row = idtable(c,:);

    %Slot mapping follows psyrat_gamma_extract_trt: oxp <- sigma_po2 (person x
    %occasion), txp <- sigma_pt2 (person x trial). The varcomps names and the Stan
    %names are transposed relative to each other, so this pairing is deliberate.
    %
    %NOTE the *_ss arguments are deliberately NOT passed, for the same reason the
    %log-nu twin gives: psyrat_ssrel_trt's per-subject path indexes them by
    %idtable.id2, which on a one-row call is that participant's ORIGINAL index
    %rather than 1, so the bounds check would fail and it would silently fall back
    %to the non-_ss values anyway. Every argument below already carries THIS
    %participant's own components, so the fallback is the correct path.
    base = {'reltype',reltype,...
        'bp',sqrt(vc.sigma_p2),'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2),...
        'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2),...
        'txo',sqrt(vc.sigma_ot2),'err',err_ss,...
        'idtable',row,'CI',ciperc,'nocc',nocc};

    %gcoeff = 2 -> generalizability G_s; gcoeff = 1 -> dependability phi_s
    genRows{c} = psyrat_ssrel_trt('gcoeff',2,base{:});
    depRows{c} = psyrat_ssrel_trt('gcoeff',1,base{:});
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});

end
