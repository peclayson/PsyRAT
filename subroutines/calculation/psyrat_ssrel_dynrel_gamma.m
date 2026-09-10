function [gen,dep] = psyrat_ssrel_dynrel_gamma(varargin)
%Per-participant CONDITIONAL dynamic/conditional reliability for the ONE-FACET
%Gamma / scaled-chi-square subject-level design (analysis 26,
%ic_dynrel_sserrvar, family = 'gamma').
%
%[gen,dep] = psyrat_ssrel_dynrel_gamma('alpha0',a0,'b',b,'lognu0',l0,'b_nu',bnu,
%  'sig_p',sp,'sig_i',si,'ind_bs',u,'ind_nu',v,'idtable',ssinfo,'ndim',1,'CI',.95)
%
%WHAT THIS IS. The subject-level sibling of psyrat_rel_dynrel_gamma. That
%function draws the POPULATION surface over a z GRID and MARGINALIZES over the
%person dispersion population (gamma estimand #2, the maintainer-directed
%group-level choice). This function instead evaluates one coefficient PER
%PARTICIPANT, at that participant's own z, CONDITIONING on their own log-mean
%(u_p = ind_bs) and log-nu (v_p = ind_nu). Structurally it is the same
%calculation with the grid loop replaced by a subject loop and the marginal
%observation term replaced by the per-person plug-in.
%
%WHY A SEPARATE FUNCTION, AND WHY NOT REUSE psyrat_ssrel WHOLESALE. Under the
%Gaussian location-scale dynrel only the RESIDUAL depends on z; sigma_p and
%sigma_i are z-invariant, so local_subj_ssrel_table_nondiff can hand psyrat_ssrel
%a single shared bp and i and vary only wp_ss per subject. Under gamma every
%observed-scale component carries exp(2*alpha(z)) (see psyrat_gamma_varcomps),
%so sigma_p2(z) and sigma_i2(z) move with z as well - and each participant sits
%at a DIFFERENT z. psyrat_ssrel takes bp and i as single columns of draws
%(psyrat_ssrel.m: 'bp = cell2mat(bp).^2', 'i2 = i(:).^2' plus a hard size guard),
%so it cannot carry a per-subject bp/i in one call.
%
%The resolution keeps the vetted arithmetic: psyrat_ssrel is called ONCE PER
%PARTICIPANT with that participant's own bp, i and residual on a one-row
%idtable, and the one-row results are stacked. Nothing in psyrat_ssrel changes,
%and the per-subject coefficient, ICC and SEM formulas stay in exactly one place.
%
%THE ESTIMAND (owner decision, 2026-07-31: "own-z, all components"). For
%participant s at their own z_s:
%
%  alpha(z_s) = alpha0 + b   . x(z_s)          nu(z_s) = exp(lognu0 + b_nu . x(z_s))
%  sigma_p2(z_s), sigma_i2(z_s)  <- psyrat_gamma_varcomps at that z
%  sigma_pi,e^2(s,z_s) = pi_mean(z_s) + 2*exp( 2*(alpha(z_s) + u_s) + 2*sig_i^2
%                               - (lognu0 + b_nu . x(z_s) + v_s) )
%  with pi_mean(z) = exp(2*alpha(z))*K, K = (exp(2 sig_p^2) - exp(sig_p^2))
%  (exp(2 sig_i^2) - exp(sig_i^2)): the population person x trial term of the
%  log-linked mean surface, pooled as under location-scale (S17-FIX; S20)
%
%and then G_s / phi_s / ICC / SEM through psyrat_ssrel exactly as the Gaussian
%path does. Two points that are easy to get wrong:
%
% (a) sigma_p2 is a BETWEEN-person variance, so it is evaluated at alpha(z_s) and
%     NOT at alpha(z_s) + u_s. Conditioning on a participant tells you where they
%     sit, not how spread the population is.
% (b) The residual drops the person facet variance sig_p^2 that the marginal
%     converter carries, because u_s is conditioned on rather than integrated
%     out. Only the WITHIN-person mean facet (2*sig_i^2, the trial main effect)
%     remains. This is the case-6 formula (psyrat_gamma_extract_sserr) with the
%     dimension slopes added; it is EXACT under the ID-only nu structure, since
%     nu is subject-constant and E_t[2*mu^2/nu_s | s] carries no lognormal
%     approximation error.
%
%(s_v, cov_pv) are passed to psyrat_gamma_varcomps as ZERO. Only sigma_p2 and
%sigma_i2 are consumed here and NEITHER depends on those two arguments, so the
%choice is numerically inconsequential; zero is passed because it states the
%conditional intent. Note the marginal-side extractors (psyrat_gamma_extract_sserr,
%psyrat_gamma_extract_trt_sserr) instead pass the estimated values through "for
%parity" with their group-level twins - also inconsequential, for the same reason.
%The two conventions differ in wording only. Estimand #2 itself lives in
%vc.e_cond_var / vc.sigma_pi_e2, which this function never reads: it forms the
%residual from the per-person plug-in below instead, and THAT is what makes the
%calculation conditional rather than marginal.
%
%A BEHAVIORAL CONSEQUENCE WORTH KNOWING, AND A DISTINCTION THAT IS EASY TO BOTCH.
%The log-mean SLOPE b still cancels here exactly as it does on the group surface:
%every observed-scale component AND the conditional residual carry exp(2*alpha(z)),
%so it divides out of the coefficient and only b_nu moves it. Verified by execution
%on this function: at b = 0 vs b = 0.5 with b_nu = 0 the per-participant G is
%identical to 0.0e+00.
%
%What does NOT cancel is the participant's OWN log-mean effect u_s. It enters the
%residual but not the population between-person numerator (which is evaluated at
%alpha(z_s), see (a) above), so a participant with a higher mean receives a
%strictly lower conditional coefficient - verified: G 0.974 -> 0.932 for
%u_s = 0 -> 0.5. That is a property of the mean-coupled variance Var = 2*mu^2/nu,
%inherited from case 6, not a defect. It means the per-participant SPREAD under
%gamma behaves qualitatively unlike the Gaussian analysis 26, whose per-subject
%residual is a pure scale effect.
%
%These are two different statements about two different quantities (a z-slope
%versus a person effect). An earlier version of this note conflated them and
%claimed the slope cancellation "does not hold" at the subject level, which is
%false; do not reintroduce that claim.
%
%Inputs (all draw vectors are columns of posterior draws)
% alpha0  - log-mean grand intercept at z = 0 (Intercept)          [ndraws x 1]
% b       - log-mean dimension slopes                           [ndraws x KDIM]
% lognu0  - log-nu grand intercept at z = 0 (Intercept_nu)         [ndraws x 1]
% b_nu    - log-nu dimension slopes                             [ndraws x KDIM]
% sig_p   - person log-mean SD (s_p)                               [ndraws x 1]
% sig_i   - trial main-effect log-mean SD (s_i)                    [ndraws x 1]
% ind_bs  - per-participant log-mean effects u_p                [ndraws x NSUB]
% ind_nu  - per-participant log-nu effects v_p                  [ndraws x NSUB]
% idtable - ssinfo: one row per participant with id, id2, z1[, z2], trls.
%           id2 indexes the ind_bs/ind_nu COLUMNS.
% ndim    - 1 or 2 dimensions
% CI      - credible interval width in decimal format (.95 = 95%)
%
%Outputs
% gen - psyrat_ssrel output table for generalizability (gcoeff = 2), one row per
%       participant in idtable row order
% dep - the same for dependability (gcoeff = 1)
%
%Both carry psyrat_ssrel's dep_*/icc_*/sem_* column names regardless of gcoeff;
%psyrat_dynrel_summary's local_pack_ssrel_nondiff renames them into the shared
%ssrel_table schema.

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
lognu0  = psyrat_req(varargin,'lognu0','Please input lognu0 (Intercept_nu draws).');
b_nu    = psyrat_req(varargin,'b_nu','Please input b_nu (log-nu dimension slopes).');
sig_p   = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_i   = psyrat_req(varargin,'sig_i','Please input sig_i (trial log-mean SD draws).');
ind_bs  = psyrat_req(varargin,'ind_bs','Please input ind_bs (per-subject log-mean effects).');
ind_nu  = psyrat_req(varargin,'ind_nu','Please input ind_nu (per-subject log-nu effects).');
idtable = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim    = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
ciperc  = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end

%column vectors of draws
alpha0 = alpha0(:);
lognu0 = lognu0(:);
sig_p  = sig_p(:);
sig_i  = sig_i(:);

%b/b_nu must be draws x KDIM (KDIM = 1 for one dimension; 3 for two), matching
%the Xdim layout the estimator builds. Same shape check psyrat_rel_dynrel_gamma
%applies, so a transposed input fails loudly rather than broadcasting.
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if size(b,2) ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        size(b,2),ndim,expectedk);
end
if size(b_nu,2) ~= expectedk
    error('varargin:bnu',...
        'WARNING: b_nu has %d columns but ndim=%d expects %d.\n',...
        size(b_nu,2),ndim,expectedk);
end

nsub = height(idtable);
if size(ind_bs,2) ~= size(ind_nu,2)
    error('varargin:indblocks',...
        'WARNING: ind_bs has %d columns but ind_nu has %d.\n',...
        size(ind_bs,2),size(ind_nu,2));
end
if any(idtable.id2 > size(ind_bs,2)) || any(idtable.id2 < 1)
    error('varargin:id2',...
        ['WARNING: idtable.id2 indexes outside the per-subject draw columns '...
        '(max id2 = %d, columns = %d).\n'],max(idtable.id2),size(ind_bs,2));
end

genRows = cell(nsub,1);
depRows = cell(nsub,1);

for c = 1:nsub
    col = idtable.id2(c);

    %dimension design row for THIS participant's z (KDIM x 1), matching the
    %Xdim layout: [z1] for one dimension, [z1; z2; z1*z2] for two.
    if ndim == 1
        xrow = idtable.z1(c);
    else
        xrow = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end

    %z-conditioned log-mean intercept and log-nu at this participant's z
    alpha_z = alpha0 + b * xrow(:);          % [ndraws x 1]
    lognu_z = lognu0 + b_nu * xrow(:);       % [ndraws x 1]
    nu_z    = exp(lognu_z);

    %observed-scale GROUP components at this participant's z. (s_v, cov_pv) = 0:
    %the conditional estimand. sigma_p2/sigma_i2 do not depend on them anyway.
    vc = psyrat_gamma_varcomps(alpha_z, sig_p, sig_i, nu_z, 0, 0);

    %S20 (pre-beta batch, 2026-09-05): pool the POPULATION person x trial
    %term of the mean surface at this z, pi_mean = var_mu - sigma_p2 -
    %sigma_i2 = exp(2*alpha(z))*K, into the participant's residual, as the
    %location-scale twin (psyrat_ssrel_dynrel_gamma_ls) does under S17-FIX.
    %The exponential alone is the conditional observation-level part; the
    %group residual carries pi_mean by subtraction, so without it every
    %per-participant coefficient was optimistically biased. Population
    %quantity: it does not touch u_s.
    pi_mean = max(vc.var_mu - vc.sigma_p2 - vc.sigma_i2, 0);    % [ndraws x 1]

    %per-participant CONDITIONAL residual at their own z and their own (u,v).
    %sig_p^2 is absent by construction: u_p is conditioned on, not integrated out.
    sig_pi_e2 = pi_mean + 2 .* exp(2.*(alpha_z + ind_bs(:,col)) + 2.*sig_i.^2 ...
        - (lognu_z + ind_nu(:,col)));        % [ndraws x 1]

    %Hand psyrat_ssrel the residual through its log-SD-equivalent contract
    %(wp = exp(wp_pop + wp_ss)^2). The whole residual is carried in wp_pop and
    %wp_ss is zero, because this call covers exactly one participant: there is no
    %population-vs-person split left to represent once z and (u,v) are both fixed.
    wp_pop = 0.5 .* log(sig_pi_e2);
    wp_ss  = zeros(size(wp_pop));

    row = idtable(c,:);

    %gcoeff = 2 -> generalizability G_s (relative error, excludes sigma_i^2)
    %gcoeff = 1 -> dependability  phi_s (absolute error, adds sigma_i^2)
    genRows{c} = psyrat_ssrel('bp',{sqrt(vc.sigma_p2)},'wp_pop',wp_pop,...
        'wp_ss',{wp_ss},'idtable',row,'CI',ciperc,...
        'i',sqrt(vc.sigma_i2),'gcoeff',2);
    depRows{c} = psyrat_ssrel('bp',{sqrt(vc.sigma_p2)},'wp_pop',wp_pop,...
        'wp_ss',{wp_ss},'idtable',row,'CI',ciperc,...
        'i',sqrt(vc.sigma_i2),'gcoeff',1);
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});

end
