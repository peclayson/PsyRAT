function mats = psyrat_dod_components(varargin)
%PSYRAT_DOD_COMPONENTS GAUSSIAN observed-scale 4x4 component covariance
% matrices for the person-specific dynamic NONCONCURRENT difference-of-
% differences designs (analyses 28/29), family = 'gaussian'. This is the
% GAUSSIAN sibling of PSYRAT_GAMMA_DOD_COMPONENTS_LS (no family marker in the
% name, the toolbox convention); the family-generic composite kernel both
% families feed is PSYRAT_DOD_COMPOSITE.
%
%   mats = psyrat_dod_components('facet','trial', ...
%       'alpha',A,'log_sigma',LS,'sd_p',SP,'cor_p',CP,'sd_trl',ST,'cor_trl',CT)
%   mats = psyrat_dod_components('facet','trial_occasion', ..., ...
%       'sd_occ',SO,'cor_occ',CO,'sd_tid',SPI,'cor_tid',CPI, ...
%       'sd_oid',SPO,'cor_oid',CPO,'sd_to',SIO,'cor_to',CIO)
%
% Required Inputs (nd = number of posterior draws; cells q = 1..4 are in
% dodmap CONSTITUENT order, the same order as the contrast weights):
%  facet     - 'trial' (analysis 28) or 'trial_occasion' (analysis 29)
%  alpha     - nd x 4 NATURAL-SCALE population expected score at the
%              evaluation point z: b(q) + Xdim(z)*b_dim(q,:)'. Identity link:
%              alpha reaches mats.expected UNTOUCHED and enters no second
%              moment (see THE IDENTITY-LINK THEOREMS below).
%  log_sigma - nd x 4 ABSOLUTE log residual SD at z (the scale submodel is
%              log-linked in both families). For the per-participant estimand
%              pass the person's own value (population + u_logsigma); for the
%              typical-person surface pass the population value. This function
%              is agnostic - person heterogeneity enters HERE, through the
%              caller's choice (same contract as the gamma sibling).
%  sd_p      - nd x 4 person LOCATION SDs on the NATURAL scale (sd_id(:,1:4);
%              log scale in the gamma sibling - the link, not the slot, sets
%              the scale)
%  cor_p     - nd x 4 x 4 person location correlation block (the [1:4,1:4]
%              sub-block of the joint 8x8 person correlation)
%  sd_trl    - nd x 4 trial SDs; cor_trl - nd x 4 x 4 trial correlations
%  (two-facet adds:) sd_occ/cor_occ  occasion block
%                    sd_tid/cor_tid  person x trial block
%                    sd_oid/cor_oid  person x occasion block
%                    sd_to /cor_to   trial x occasion block
%
% Output struct mats. Every component is an nd x 4 x 4 SYMMETRIC matrix
% (diagonal = per-cell variances, off-diagonals = signed cross-cell
% covariances), ready for PSYRAT_DOD_COMPOSITE. Same keys as the gamma
% sibling:
%
%  one-facet:  mats.p, mats.i, mats.pi_e
%  two-facet:  mats.p, mats.i, mats.o, mats.pi, mats.po, mats.io, mats.pio_e
%  both:       mats.expected (nd x 4 population expected scores at z = alpha
%              itself under the identity link),
%              mats.v_total  (nd x 4 summed random-effect variance per cell,
%              residual excluded - the same structural meaning as the gamma
%              sibling's log-scale v_total)
%
% THE IDENTITY-LINK THEOREMS (SCIENTIFIC_FORMULA_AUDIT.md section 26), which
% are what make this assembler nearly arithmetic-free:
%  (1) Signal components are the model covariance blocks DIRECTLY:
%      variances sd_*.^2, covariances sd_a.*sd_b.*cor_ab. Location parameters
%      (alpha, and z through it) enter no second moment, so every signal
%      matrix is constant in z; the reliability surface varies with z only
%      through the residual diagonal. The section-25 lognormal moment
%      identities are log-link arithmetic and MUST NOT be applied here.
%  (2) The residual channel is EXACTLY diagonal: diag(exp(2*log_sigma)). The
%      identity link induces no mean-surface interaction (the gamma pi_mean
%      remainder is identically zero here), and the six observation-level
%      residual covariances are zero by estimand, so the pooled cross channel
%      vanishes identically. The off-diagonals below are LITERAL zeros, never
%      a computed difference - the bit-identity the conversion tests pin.
%
% THE DEFINING ASSUMPTION carries over from the gamma sibling: all six
% pairwise OBSERVATION-level residual covariances are fixed to exactly zero
% (the nonconcurrent, counterfactual estimand). Person/facet covariances are a
% different level of variation and are retained with their signs.
%
% WHY THERE IS NO s_sd INPUT. The person-specific estimand conditions on each
% participant's own residual SD; the caller passes it through log_sigma.
% Marginalizing over person scale heterogeneity is a DIFFERENT estimand and is
% deliberately unreachable from this function (same rule as the gamma
% sibling).
%
% See also PSYRAT_DOD_COMPOSITE, PSYRAT_GAMMA_DOD_COMPONENTS_LS.

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
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_dod_components for more information'));
end

%The nonconcurrent boundary, enforced at the signature. This codebase's
%name-value pattern ignores unknown names, which here would let a residual
%coupling argument be SILENTLY swallowed - the one confusion this design must
%never permit. Reject the coupling vocabulary by name instead.
for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
    if any(strcmpi(coupled{1},varargin(1:2:end)))
        error('psyrat_dod_components:noresidualcoupling',... %Error code and associated error
            ['''%s'' is not an input: the nonconcurrent DoD estimand fixes '...
            'all six observation-level residual covariances to exactly zero '...
            'by assumption. There is no residual coupling parameter to pass. '...
            'For concurrent designs use the analysis-19/20 machinery.'],...
            coupled{1});
    end
end

facet = local_req(varargin,'facet','The facet structure (facet) was not specified.');
alpha = local_req(varargin,'alpha','Expected scores (alpha) were not specified.');
log_sigma = local_req(varargin,'log_sigma','Absolute log residual SDs (log_sigma) were not specified.');
sd_p = local_req(varargin,'sd_p','Person location SDs (sd_p) were not specified.');
cor_p = local_req(varargin,'cor_p','Person location correlations (cor_p) were not specified.');
sd_trl = local_req(varargin,'sd_trl','Trial SDs (sd_trl) were not specified.');
cor_trl = local_req(varargin,'cor_trl','Trial correlations (cor_trl) were not specified.');

if ~any(strcmpi(facet,{'trial','trial_occasion'}))
    error('varargin:facet',... %Error code and associated error
        'facet must be ''trial'' or ''trial_occasion''.');
end
istwofacet = strcmpi(facet,'trial_occasion');

if istwofacet
    sd_occ = local_req(varargin,'sd_occ','Occasion SDs (sd_occ) were not specified.');
    cor_occ = local_req(varargin,'cor_occ','Occasion correlations (cor_occ) were not specified.');
    sd_tid = local_req(varargin,'sd_tid','Person x trial SDs (sd_tid) were not specified.');
    cor_tid = local_req(varargin,'cor_tid','Person x trial correlations (cor_tid) were not specified.');
    sd_oid = local_req(varargin,'sd_oid','Person x occasion SDs (sd_oid) were not specified.');
    cor_oid = local_req(varargin,'cor_oid','Person x occasion correlations (cor_oid) were not specified.');
    sd_to = local_req(varargin,'sd_to','Trial x occasion SDs (sd_to) were not specified.');
    cor_to = local_req(varargin,'cor_to','Trial x occasion correlations (cor_to) were not specified.');
end

nd = size(alpha,1);
local_check_draws(alpha,nd,'alpha');
local_check_draws(log_sigma,nd,'log_sigma');
local_check_draws(sd_p,nd,'sd_p');
local_check_corr(cor_p,nd,'cor_p');
local_check_draws(sd_trl,nd,'sd_trl');
local_check_corr(cor_trl,nd,'cor_trl');
if istwofacet
    local_check_draws(sd_occ,nd,'sd_occ');
    local_check_corr(cor_occ,nd,'cor_occ');
    local_check_draws(sd_tid,nd,'sd_tid');
    local_check_corr(cor_tid,nd,'cor_tid');
    local_check_draws(sd_oid,nd,'sd_oid');
    local_check_corr(cor_oid,nd,'cor_oid');
    local_check_draws(sd_to,nd,'sd_to');
    local_check_corr(cor_to,nd,'cor_to');
end

%Identity link: expected is alpha itself, and the signal blocks assemble
%directly from their SDs and correlations (theorem 1 in the header).
mats = struct;
mats.expected = alpha;
if istwofacet
    blocks = {'p',sd_p,cor_p; 'i',sd_trl,cor_trl; 'o',sd_occ,cor_occ; ...
        'pi',sd_tid,cor_tid; 'po',sd_oid,cor_oid; 'io',sd_to,cor_to};
    reskey = 'pio_e';
    mats.v_total = sd_p.^2 + sd_trl.^2 + sd_occ.^2 + sd_tid.^2 + ...
        sd_oid.^2 + sd_to.^2;
else
    blocks = {'p',sd_p,cor_p; 'i',sd_trl,cor_trl};
    reskey = 'pi_e';
    mats.v_total = sd_p.^2 + sd_trl.^2;
end

for k = 1:size(blocks,1)
    key = blocks{k,1}; sd = blocks{k,2}; cor = blocks{k,3};
    M = zeros(nd,4,4);
    for q = 1:4
        M(:,q,q) = sd(:,q).^2;
    end
    for a = 1:3
        for b = (a+1):4
            cv = sd(:,a) .* sd(:,b) .* cor(:,a,b);
            M(:,a,b) = cv;
            M(:,b,a) = cv;
        end
    end
    mats.(key) = M;
end

%The residual channel (theorem 2 in the header): diagonal = the cell's own
%sigma^2; off-diagonals stay the LITERAL zeros this matrix was born with -
%never a computed difference. The conversion tests pin that as a bit-identity.
R = zeros(nd,4,4);
for q = 1:4
    R(:,q,q) = exp(2 .* log_sigma(:,q));
end
mats.(reskey) = R;

end

function val = local_req(args,name,msg)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    error(['varargin:' lower(name)],... %Error code and associated error
        msg);
end
val = args{ind+1};
end

function local_check_draws(x,nd,name)
if ~isnumeric(x) || ~ismatrix(x) || size(x,1) ~= nd || size(x,2) ~= 4
    error('psyrat_dod_components:size',... %Error code and associated error
        '%s must be an ndraws x 4 numeric matrix (%d draws expected).',name,nd);
end
end

function local_check_corr(x,nd,name)
if ~isnumeric(x) || ndims(x) ~= 3 || size(x,1) ~= nd || ...
        size(x,2) ~= 4 || size(x,3) ~= 4
    error('psyrat_dod_components:size',... %Error code and associated error
        '%s must be an ndraws x 4 x 4 numeric array (%d draws expected).',name,nd);
end
%A covariance matrix passed where a CORRELATION matrix belongs is the likeliest
%scale error on this signature; unit diagonals catch it cheaply.
for q = 1:4
    if any(abs(x(:,q,q) - 1) > 1e-8)
        error('psyrat_dod_components:notcorrelation',... %Error code and associated error
            ['%s must be a CORRELATION matrix (unit diagonal); its [%d,%d] '...
            'diagonal is not 1. How to fix: pass the correlation block, not '...
            'the covariance block; SDs enter through the sd_* inputs.'],...
            name,q,q);
    end
end
end
