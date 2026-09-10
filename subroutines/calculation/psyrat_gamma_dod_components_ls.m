function mats = psyrat_gamma_dod_components_ls(varargin)
%PSYRAT_GAMMA_DOD_COMPONENTS_LS Observed-scale 4x4 component covariance
% matrices for the person-specific dynamic NONCONCURRENT difference-of-
% differences designs (analyses 28/29), location-scale gamma family only.
%
%   mats = psyrat_gamma_dod_components_ls('facet','trial', ...
%       'alpha',A,'log_sigma',LS,'sd_p',SP,'cor_p',CP,'sd_trl',ST,'cor_trl',CT)
%   mats = psyrat_gamma_dod_components_ls('facet','trial_occasion', ..., ...
%       'sd_occ',SO,'cor_occ',CO,'sd_tid',SPI,'cor_tid',CPI, ...
%       'sd_oid',SPO,'cor_oid',CPO,'sd_to',SIO,'cor_to',CIO)
%
% Required Inputs (nd = number of posterior draws; cells q = 1..4 are in
% dodmap CONSTITUENT order, the same order as the contrast weights):
%  facet     - 'trial' (analysis 28) or 'trial_occasion' (analysis 29)
%  alpha     - nd x 4 ABSOLUTE log expected score at the evaluation point z:
%              mu_offset(q) + b(q) + Xdim(z)*b_dim(q,:)'
%  log_sigma - nd x 4 ABSOLUTE log residual SD at z. For the per-participant
%              estimand pass the person's own value (population + u_logsigma);
%              for the typical-person surface pass the population value. This
%              function is agnostic - person heterogeneity enters HERE, through
%              the caller's choice, never through a marginalizing s_sd (see
%              WHY THERE IS NO s_sd below).
%  sd_p      - nd x 4 person LOCATION SDs on the log scale (sd_id(:,1:4))
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
% covariances), ready for PSYRAT_DOD_COMPOSITE. Component keys follow the
% REFERENCE bundle's vocabulary (reliability.R; frozen in
% tests/baselines/nonconcurrent_dod), with the toolbox facet-name mapping
% applied HERE and only here:
%
%   toolbox input   reference component      meaning
%   sd_trl          i                        trial main effect
%   sd_occ          o                        occasion main effect
%   sd_tid          pi                       person x trial
%   sd_oid          po                       person x occasion
%   sd_to           io                       trial x occasion
%
%  one-facet:  mats.p, mats.i, mats.pi_e
%  two-facet:  mats.p, mats.i, mats.o, mats.pi, mats.po, mats.io, mats.pio_e
%  both:       mats.expected (nd x 4 population expected scores at z),
%              mats.v_total  (nd x 4 total log-scale variance per cell)
%
% THE DEFINING ASSUMPTION, IMPLEMENTED AS A BIT-IDENTITY. All six pairwise
% OBSERVATION-level residual covariances are fixed to exactly zero (the
% nonconcurrent estimand). This function therefore calls the concurrent cross
% helpers with cov_e_obs = 0, which makes the pooled residual cross channel
% (cov_pi_e_obs / cov_res_e_obs) bit-identical to the mean-surface-only cross
% term - the residual off-diagonals carry ONLY the person/facet covariance the
% model retains, never observation-level coupling. The diagonal residual is
% the pooled channel of PSYRAT_GAMMA_VARCOMPS(_TRT)_LS: the mean-surface
% interaction remainder plus the cell's own sigma^2 (the S17 pooling
% convention).
%
% WHY THERE IS NO s_sd INPUT. The person-specific estimand conditions on each
% participant's own residual SD; the caller passes it through log_sigma. The
% marginal-over-persons residual (s_sd > 0 in the underlying converters) is a
% DIFFERENT estimand and is deliberately unreachable from this function -
% adding s_sd "for symmetry" would silently blend the two.
%
% NUMERICS ARE THE CONVERTERS'. Every diagonal comes from
% PSYRAT_GAMMA_VARCOMPS_LS / _TRT_LS and every off-diagonal from
% PSYRAT_GAMMA_CROSSCOV_RESCOR / _TRT_RESCOR, so all shared quantities are
% bit-identical to the shipped analysis-19/20 path, including the converters'
% tiny-negative handling on the subtraction components. This file adds no
% arithmetic of its own beyond assembling matrices.
%
% See also PSYRAT_DOD_COMPOSITE, PSYRAT_GAMMA_VARCOMPS_LS,
% PSYRAT_GAMMA_VARCOMPS_TRT_LS, PSYRAT_GAMMA_CROSSCOV_RESCOR,
% PSYRAT_GAMMA_CROSSCOV_TRT_RESCOR.

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
        'See help psyrat_gamma_dod_components_ls for more information'));
end

%The nonconcurrent boundary, enforced at the signature. This codebase's
%name-value pattern ignores unknown names, which here would let a residual
%coupling argument be SILENTLY swallowed - the one confusion this design must
%never permit. Reject the coupling vocabulary by name instead.
for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
    if any(strcmpi(coupled{1},varargin(1:2:end)))
        error('psyrat_gamma_dod_components_ls:noresidualcoupling',... %Error code and associated error
            ['''%s'' is not an input: the nonconcurrent DoD estimand fixes '...
            'all six observation-level residual covariances to exactly zero '...
            'by assumption. There is no residual coupling parameter to pass. '...
            'For concurrent designs use the analysis-19/20 machinery.'],...
            coupled{1});
    end
end

facet = local_req(varargin,'facet','The facet structure (facet) was not specified.');
alpha = local_req(varargin,'alpha','Absolute log means (alpha) were not specified.');
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

z44 = zeros(nd,4,4);
mats = struct;
mats.expected = zeros(nd,4);
mats.v_total = zeros(nd,4);
if istwofacet
    keys = {'p','i','o','pi','po','io','pio_e'};
else
    keys = {'p','i','pi_e'};
end
for k = 1:numel(keys)
    mats.(keys{k}) = z44;
end

%Diagonals: one converter call per cell, all shared quantities bit-identical to
%the shipped analysis-19/20 path.
for q = 1:4
    if istwofacet
        %converter argument order is (alpha, s_p, s_o, s_t, s_po, s_pt, s_ot,
        %log_sigma); the toolbox facet names map as documented in the header:
        %occ->o, trl->t, oid->po, tid->pt, to->ot
        vc = psyrat_gamma_varcomps_trt_ls(alpha(:,q), sd_p(:,q), ...
            sd_occ(:,q), sd_trl(:,q), sd_oid(:,q), sd_tid(:,q), ...
            sd_to(:,q), log_sigma(:,q));
        mats.p(:,q,q) = vc.sigma_p2;
        mats.i(:,q,q) = vc.sigma_t2;
        mats.o(:,q,q) = vc.sigma_o2;
        mats.pi(:,q,q) = vc.sigma_pt2;
        mats.po(:,q,q) = vc.sigma_po2;
        mats.io(:,q,q) = vc.sigma_ot2;
        mats.pio_e(:,q,q) = vc.sigma_res2;
        v_total = sd_p(:,q).^2 + sd_occ(:,q).^2 + sd_trl(:,q).^2 + ...
            sd_oid(:,q).^2 + sd_tid(:,q).^2 + sd_to(:,q).^2;
    else
        vc = psyrat_gamma_varcomps_ls(alpha(:,q), sd_p(:,q), ...
            sd_trl(:,q), log_sigma(:,q));
        mats.p(:,q,q) = vc.sigma_p2;
        mats.i(:,q,q) = vc.sigma_i2;
        mats.pi_e(:,q,q) = vc.sigma_pi_e2;
        v_total = sd_p(:,q).^2 + sd_trl(:,q).^2;
    end
    mats.expected(:,q) = vc.mu_bar;
    mats.v_total(:,q) = v_total;
end

%Off-diagonals: one concurrent-helper call per pair with cov_e_obs = 0 - the
%nonconcurrent boundary as a bit-identity (see header).
for a = 1:3
    for b = (a+1):4
        if istwofacet
            cv = psyrat_gamma_crosscov_trt_rescor( ...
                alpha(:,a), sd_p(:,a), sd_occ(:,a), sd_trl(:,a), ...
                sd_oid(:,a), sd_tid(:,a), sd_to(:,a), ...
                alpha(:,b), sd_p(:,b), sd_occ(:,b), sd_trl(:,b), ...
                sd_oid(:,b), sd_tid(:,b), sd_to(:,b), ...
                cor_p(:,a,b), cor_occ(:,a,b), cor_trl(:,a,b), ...
                cor_oid(:,a,b), cor_tid(:,a,b), cor_to(:,a,b), 0);
            vals = {cv.cov_p, cv.cov_t, cv.cov_o, cv.cov_pt, ...
                cv.cov_po, cv.cov_ot, cv.cov_res_e_obs};
        else
            cv = psyrat_gamma_crosscov_rescor( ...
                alpha(:,a), sd_p(:,a), sd_trl(:,a), ...
                alpha(:,b), sd_p(:,b), sd_trl(:,b), ...
                cor_p(:,a,b), cor_trl(:,a,b), 0);
            vals = {cv.cov_p_obs, cv.cov_i_obs, cv.cov_pi_e_obs};
        end
        for k = 1:numel(keys)
            mats.(keys{k})(:,a,b) = vals{k};
            mats.(keys{k})(:,b,a) = vals{k};
        end
    end
end

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
    error('psyrat_gamma_dod_components_ls:size',... %Error code and associated error
        '%s must be an ndraws x 4 numeric matrix (%d draws expected).',name,nd);
end
end

function local_check_corr(x,nd,name)
if ~isnumeric(x) || ndims(x) ~= 3 || size(x,1) ~= nd || ...
        size(x,2) ~= 4 || size(x,3) ~= 4
    error('psyrat_gamma_dod_components_ls:size',... %Error code and associated error
        '%s must be an ndraws x 4 x 4 numeric array (%d draws expected).',name,nd);
end
%A covariance matrix passed where a CORRELATION matrix belongs is the likeliest
%scale error on this signature; unit diagonals catch it cheaply.
for q = 1:4
    if any(abs(x(:,q,q) - 1) > 1e-8)
        error('psyrat_gamma_dod_components_ls:notcorrelation',... %Error code and associated error
            ['%s must be a CORRELATION matrix (unit diagonal); its [%d,%d] '...
            'diagonal is not 1. How to fix: pass the correlation block, not '...
            'the covariance block; SDs enter through the sd_* inputs.'],...
            name,q,q);
    end
end
end
