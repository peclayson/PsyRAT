function [ssrel,drawlevel] = psyrat_ssrel_dodiffdynrel(varargin)
%Per-participant person-specific dynamic NONCONCURRENT difference-of-
%differences reliability for the one-facet GAUSSIAN design (analysis 28,
%ic_dodiff_dynrel_sserrvar, family = 'gaussian'). The Gaussian sibling of
%PSYRAT_SSREL_DODIFFDYNREL_GAMMA_LS (no family marker in the name, the
%toolbox convention). Derivation: SCIENTIFIC_FORMULA_AUDIT.md section 26.
%
%[ssrel,drawlevel] = psyrat_ssrel_dodiffdynrel('b',b,...)
%
%WHAT THIS IS. The same estimand as the gamma sibling - the reliability of
%C = (E1 - E2) - (E3 - E4), w = (+1,-1,-1,+1), with the population signal /
%person residual split, the nonconcurrent boundary (all six observation-level
%residual covariances fixed to exactly zero), and two D-study designs per
%participant - read that file first. What the IDENTITY LINK changes
%(section 26):
%
% (1) COMPONENTS ARE THE MODEL BLOCKS DIRECTLY. No lognormal moment
%     conversion exists or may be applied; the signal matrices are constant
%     in z and the residual channel is exactly diagonal, so the per-
%     participant relative error is EXACTLY sum_q c_q^2*sigma_qp^2/n_q.
% (2) EXPECTED SCORES ARE THE ADDITIVE PREDICTOR. E_own = the participant's
%     own location effect plus the slopes at their z - no +0.5*variance
%     lognormal mean correction (a log-link artifact).
% (3) NO nu READ-OUT. nu = 2*(mu/sigma)^2 is a scaled-chi-square quantity;
%     the table drops the nu1..4_pt columns and keeps ressd1..4_pt. No
%     downstream consumer reads the nu columns (verified 2026-08-20).
%
%INPUTS (name-value; cells q = 1..4 in dodmap CONSTITUENT order)
% b            - [ndraws x 4] per-cell NATURAL-SCALE mean intercepts at z = 0
%                (identity link; no offsets exist for this family)
% b_sigma      - [ndraws x 4] ABSOLUTE per-cell log residual-SD intercepts
% b_dim        - [ndraws x 4 x KDIM] per-cell mean dimension slopes
% b_sigma_dim  - [ndraws x 4 x KDIM] per-cell log residual-SD dimension slopes
% sd_id        - [ndraws x 8] joint person SDs [loc_1..4, logsd_1..4]
%                (locations on the natural scale)
% cor_id       - [ndraws x 8 x 8] joint person correlation draws
% sd_trl       - [ndraws x 4] per-cell trial SDs
% cor_trl      - [ndraws x 4 x 4] cross-cell trial correlation draws
% mu_ss1..mu_ss4 - [ndraws x NSUB] per-subject NATURAL-SCALE mean at z = 0
%                (b(q) + person location effect)
% er_ss1..er_ss4 - [ndraws x NSUB] per-subject ABSOLUTE log residual SD at
%                z = 0 (b_sigma(q) + person scale effect)
% idtable      - ssinfo table: id, id2, z1 (,z2), trls1..trls4
% nprime       - [1 x 4] standardized per-cell trial counts n'
% ndim         - 1 or 2
% CI           - credible interval width in (0,1)
% cvec         - optional signed contrast (default [1 -1 -1 1])
%
%OUTPUTS - as the gamma sibling minus the nu columns:
% ssrel     - one row per participant: id, id2, z1(,z2), trls1..trls4, the
%             shared coefficient columns (gen_*/dep_*/iccg_*/iccd_*/
%             sem_gen_*/sem_dep_* with *_std twins), uni_pt, relerr_pt,
%             abserr_pt (+_std), expected_dod_pop_pt, expected_dod_person_pt,
%             residual_var_dod_person_pt, ressd1_pt..ressd4_pt,
%             invalid_frac, invalid_frac_std, n_floored.
% drawlevel - OPTIONAL per-draw arrays: u, rel, abs_, G, D, iccg, iccd,
%             rel_std, abs_std, G_std, D_std, invalid, invalid_std, plus
%             per-cell [ndraws x 4 x nsub] expected_person, expected_pop,
%             sigma_person.
%
%See also PSYRAT_SSREL_DODIFFDYNREL_GAMMA_LS, PSYRAT_DOD_COMPOSITE,
%PSYRAT_DOD_COMPONENTS, PSYRAT_REL_DODIFFDYNREL.

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

%The nonconcurrent boundary, enforced at the signature exactly as in
%PSYRAT_DOD_COMPONENTS: coupling vocabulary is rejected, never silently
%swallowed.
for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
    if any(strcmpi(coupled{1},varargin(1:2:end)))
        error('psyrat_ssrel_dodiffdynrel:noresidualcoupling',... %Error code and associated error
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
cor_trl     = psyrat_req(varargin,'cor_trl','Please input cor_trl (cross-cell trial correlations, draws x 4 x 4).');
mu_ss       = local_ss_block(varargin,'mu_ss','per-subject natural-scale means');
er_ss       = local_ss_block(varargin,'er_ss','per-subject absolute log residual SDs');
idtable     = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
nprime      = psyrat_req(varargin,'nprime','Please input nprime (standardized per-cell trial counts).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
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
ndraws = size(b,1);
if size(b,2) ~= 4 || size(b_sigma,2) ~= 4 || size(sd_trl,2) ~= 4
    error('varargin:cells',... %Error code and associated error
        'WARNING: b, b_sigma and sd_trl must each have 4 columns, one per cell.\n');
end
if size(sd_id,2) ~= 8
    error('varargin:sd_id',... %Error code and associated error
        ['WARNING: sd_id must have 8 columns, ordered [loc_1..4, logsd_1..4]; '...
        'got %d.\n'],size(sd_id,2));
end
kdim = 2*ndim - 1;
if size(b_dim,2) ~= 4 || size(b_dim,3) ~= kdim || ...
        size(b_sigma_dim,2) ~= 4 || size(b_sigma_dim,3) ~= kdim
    error('varargin:b_dim',... %Error code and associated error
        ['WARNING: b_dim and b_sigma_dim must be [ndraws x 4 x %d] for '...
        'ndim = %d.\n'],kdim,ndim);
end
if numel(nprime) ~= 4 || any(nprime <= 0) || any(isnan(nprime))
    error('varargin:nprime',... %Error code and associated error
        'WARNING: nprime must be four positive per-cell trial counts.\n');
end
nprime = nprime(:)';

nsub = height(idtable);
for q = 1:4
    if size(mu_ss{q},2) < max(idtable.id2) || size(er_ss{q},2) < max(idtable.id2)
        %FAIL LOUDLY: a short per-subject matrix would silently move the
        %estimand from this participant's own residual to another's.
        error('varargin:ss',... %Error code and associated error
            ['WARNING: mu_ss%d/er_ss%d have %d/%d subject columns but idtable '...
            'references id2 up to %d.\n'],q,q,...
            size(mu_ss{q},2),size(er_ss{q},2),max(idtable.id2));
    end
end

%person LOCATION block: SDs 1..4 and the 4x4 location sub-block of the joint
%8x8 correlation
sd_p = sd_id(:,1:4);
cor_p = cor_id(:,1:4,1:4);

wantdraw = nargout > 1;
if wantdraw
    z1nan = nan(ndraws,nsub);
    drawlevel = struct('u',z1nan,'rel',z1nan,'abs_',z1nan,'G',z1nan,...
        'D',z1nan,'iccg',z1nan,'iccd',z1nan,'rel_std',z1nan,...
        'abs_std',z1nan,'G_std',z1nan,'D_std',z1nan,...
        'invalid',z1nan,'invalid_std',z1nan,...
        'expected_person',nan(ndraws,4,nsub),'expected_pop',nan(ndraws,4,nsub),...
        'sigma_person',nan(ndraws,4,nsub));
end

rows = cell(nsub,1);
ciedge = (1-ciperc)/2;

for c = 1:nsub

    if ndim == 1
        zt = idtable.z1(c);
    else
        zt = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end
    zt = zt(:);
    col = idtable.id2(c);

    %POPULATION means at this participant's z, and the PARTICIPANT'S OWN
    %values (their realized effects + the same slopes). The population log
    %residual SD is deliberately NOT formed here: nothing in this estimand
    %consumes it (the signal is population-at-z, the residual is the
    %participant's own).
    alpha = zeros(ndraws,4);
    lsig_own = zeros(ndraws,4);
    mu_own = zeros(ndraws,4);
    for q = 1:4
        dimterm = reshape(b_dim(:,q,:),ndraws,kdim) * zt;
        sigterm = reshape(b_sigma_dim(:,q,:),ndraws,kdim) * zt;
        alpha(:,q) = b(:,q) + dimterm;
        mu_own(:,q) = mu_ss{q}(:,col) + dimterm;
        lsig_own(:,q) = er_ss{q}(:,col) + sigterm;
    end

    %Components: population signal at z, the participant's OWN residual.
    mats = psyrat_dod_components('facet','trial',...
        'alpha',alpha,'log_sigma',lsig_own,...
        'sd_p',sd_p,'cor_p',cor_p,'sd_trl',sd_trl,'cor_trl',cor_trl);

    %floor exactly the empty-cell ZERO (B16): a participant with no trials
    %in a cell would trip psyrat_dod_composite's n > 0 guard and kill the
    %whole per-participant summary. ONLY zero is floored -- deliberately
    %NOT max(1,...): MATLAB's two-input max omits NaN, so a max floor
    %would silently convert corrupt counts (NaN, negative) to n = 1 and
    %hide them from the composite guard, which must keep erroring loudly
    %on them (the reject-before-floor rationale psyrat_dynrel_summary
    %documents for its own n' floor). The floor applies to the composite
    %input only; ssrel keeps the true counts, and preflight
    %(dodPersonMissingCell) remains the gate on the GUI path.
    n_act = [idtable.trls1(c) idtable.trls2(c) idtable.trls3(c) idtable.trls4(c)];
    %B16 residual: record HOW MANY cells the floor substituted, so the export
    %discloses it. Counted BEFORE the floor is applied, and the counts
    %themselves stay untouched -- the true zero is what the table must keep.
    n_floored = nnz(n_act == 0);
    n_act(n_act == 0) = 1;

    u    = psyrat_dod_composite('varcov',mats.p,'cvec',cvec,'rule','unscaled');
    relu = psyrat_dod_composite('varcov',mats.pi_e,'cvec',cvec,'rule','unscaled');
    iu   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','unscaled');
    rel_a = psyrat_dod_composite('varcov',mats.pi_e,'cvec',cvec,'rule','harmonic','n',n_act);
    i_a   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','harmonic','n',n_act);
    rel_s = psyrat_dod_composite('varcov',mats.pi_e,'cvec',cvec,'rule','harmonic','n',nprime);
    i_s   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','harmonic','n',nprime);

    [coefA,invA] = local_coeffs(u,rel_a,i_a);
    [coefS,invS] = local_coeffs(u,rel_s,i_s);
    iccg = local_ratio(u.value, u.value + relu.value);
    iccd = local_ratio(u.value, u.value + relu.value + iu.value);

    %Per-participant expected scores: the additive predictor itself under
    %the identity link (no lognormal mean correction; trial and residual
    %effects have mean zero additively). No nu: a scaled-chi-square
    %quantity with no Gaussian counterpart.
    sig_own = exp(lsig_own);
    E_own = mu_own;
    expected_dod_pop = mats.expected * cvec(:);
    expected_dod_own = E_own * cvec(:);
    %the four-term identity (section 26): exactly the sum, no covariance
    %corrections, because the six residual covariances are zero by estimand
    res_var_dod = sum(sig_own.^2, 2);

    if wantdraw
        drawlevel.u(:,c) = u.value;
        drawlevel.rel(:,c) = coefA.rel;   drawlevel.abs_(:,c) = coefA.abs_;
        drawlevel.G(:,c) = coefA.G;       drawlevel.D(:,c) = coefA.D;
        drawlevel.iccg(:,c) = iccg;       drawlevel.iccd(:,c) = iccd;
        drawlevel.rel_std(:,c) = coefS.rel; drawlevel.abs_std(:,c) = coefS.abs_;
        drawlevel.G_std(:,c) = coefS.G;   drawlevel.D_std(:,c) = coefS.D;
        drawlevel.invalid(:,c) = invA;    drawlevel.invalid_std(:,c) = invS;
        drawlevel.expected_person(:,:,c) = E_own;
        drawlevel.expected_pop(:,:,c) = mats.expected;
        drawlevel.sigma_person(:,:,c) = sig_own;
    end

    row = idtable(c,:);
    row = local_summarize(row,'gen',coefA.G,ciedge);
    row = local_summarize(row,'dep',coefA.D,ciedge);
    row = local_summarize(row,'iccg',iccg,ciedge);
    row = local_summarize(row,'iccd',iccd,ciedge);
    row = local_summarize(row,'sem_gen',local_sem(coefA.rel),ciedge);
    row = local_summarize(row,'sem_dep',local_sem(coefA.abs_),ciedge);
    row = local_summarize(row,'gen',coefS.G,ciedge,'_std');
    row = local_summarize(row,'dep',coefS.D,ciedge,'_std');
    row = local_summarize(row,'sem_gen',local_sem(coefS.rel),ciedge,'_std');
    row = local_summarize(row,'sem_dep',local_sem(coefS.abs_),ciedge,'_std');
    row.uni_pt = mean(u.value,'omitnan');
    row.relerr_pt = mean(coefA.rel,'omitnan');
    row.abserr_pt = mean(coefA.abs_,'omitnan');
    row.relerr_pt_std = mean(coefS.rel,'omitnan');
    row.abserr_pt_std = mean(coefS.abs_,'omitnan');
    row.expected_dod_pop_pt = mean(expected_dod_pop,'omitnan');
    row.expected_dod_person_pt = mean(expected_dod_own,'omitnan');
    row.residual_var_dod_person_pt = mean(res_var_dod,'omitnan');
    for q = 1:4
        row.(sprintf('ressd%d_pt',q)) = mean(sig_own(:,q),'omitnan');
    end
    row.invalid_frac = mean(invA);
    row.invalid_frac_std = mean(invS);
    row.n_floored = n_floored;
    rows{c} = row;
end

ssrel = vertcat(rows{:});

end

function ss = local_ss_block(args,stem,what)
%Collect mu_ss1..mu_ss4 / er_ss1..er_ss4 into a 4-cell block.
ss = cell(1,4);
for q = 1:4
    ss{q} = psyrat_req(args,sprintf('%s%d',stem,q),...
        sprintf('Please input %s%d (%s, cell %d).',stem,q,what,q));
end
end

function [coef,invalid] = local_coeffs(u,relc,ic)
%Coefficient assembly shared with the gamma sibling: error terms from the
%scaled composites, coefficients as ratios, and a per-draw invalid flag
%pooling every component's material-negative flag with the coefficient
%range check (tolerance 1e-10).
coef.rel = relc.value;
coef.abs_ = relc.value + ic.value;
coef.G = local_ratio(u.value, u.value + coef.rel);
coef.D = local_ratio(u.value, u.value + coef.abs_);
tol = 1e-10;
badcoef = @(x) ~isfinite(x) | x < -tol | x > 1 + tol;
invalid = double(u.invalid_negative | relc.invalid_negative | ...
    ic.invalid_negative | badcoef(coef.G) | badcoef(coef.D));
end

function r = local_ratio(num,den)
%NaN where undefined or structurally negative, never a repaired value.
r = num ./ den;
r(~isfinite(num) | ~isfinite(den) | num < 0 | den <= 0) = NaN;
end

function s = local_sem(v)
%The square root where the error variance is finite and nonnegative, NaN
%otherwise - a materially negative error term yields NaN + a flag, never a
%silent zero.
s = nan(size(v));
ok = isfinite(v) & v >= 0;
s(ok) = sqrt(v(ok));
end

function row = local_summarize(row,stem,draws,ciedge,suffix)
%NaN draws (structurally undefined coefficients) are omitted; the
%invalid_frac columns carry the story. The optional suffix places a design
%tag AFTER the statistic (gen_pt_std), matching the relerr_pt_std convention.
if nargin < 5, suffix = ''; end
q = quantile(draws,[ciedge 1-ciedge]);   %MATLAB quantile ignores NaNs
row.(sprintf('%s_pt%s',stem,suffix)) = mean(draws,'omitnan');
row.(sprintf('%s_ll%s',stem,suffix)) = q(1);
row.(sprintf('%s_ul%s',stem,suffix)) = q(2);
end
