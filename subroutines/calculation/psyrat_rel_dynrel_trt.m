function out = psyrat_rel_dynrel_trt(varargin)
%Dynamic/conditional trial + occasion two-facet reliability as a function of the
%standardized dimension predictors (Rast & Clayson, analysis 14).
%
% out = psyrat_rel_dynrel_trt('sig_p',sig_p,'sig_log0',sig_log0,...
%   'b_sigma',b_sigma,'sig_occ',sig_occ,'sig_trl',sig_trl,...
%   'sig_trlxid',sig_trlxid,'sig_occxid',sig_occxid,'sig_trlxocc',sig_trlxocc,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',ntrl,'nocc',nocc,...
%   'reltype',reltype,'CI',.95)
%
%The two-facet dynamic-reliability model estimates a dimension-conditioned
%residual SD (the highest-order person x trial x occasion residual)
%
%   sigma_pto,e(z) = exp(Intercept_sigma + b_sigma . X(z))      (delta_p = 0)
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two. Reliability is then the two-facet coefficient evaluated at each z,
%reusing the validated test-retest psyrat_rel_trt with the dimension-conditioned
%residual feeding its 'err' argument and the five crossed facet SDs feeding the
%static facet arguments. psyrat_rel_trt's three reltypes give:
%
%   reltype 1 (coefficient of equivalence / internal consistency: occasion fixed)
%   reltype 2 (coefficient of stability: trial fixed)
%   reltype 3 (coefficient of trial equivalence and stability: both random)
%
%which is exactly the brms-notes G(z)/D(z) for reltype 3 (the dynamic G_t) and
%its occasion-fixed / trial-fixed analogues for reltypes 1/2. Generalizability
%(gcoeff 2) excludes, and dependability (gcoeff 1) includes, the facet main
%effects, exactly as the static test-retest path.
%
%Inputs (name-value)
% sig_p      - between-person (universe-score) SD posterior draws (draws x 1).
%              This is gro_sds(:,1); per the owner decision it is CONDITIONAL on
%              the dimensions (they enter the mean submodel).
% sig_log0   - log-residual intercept draws Intercept_sigma (draws x 1): the log
%              residual SD at z = 0. z = 0 is the OVERALL-sample mean of each
%              dimension (standardized once across the whole sample, pooled
%              across any groups), not a per-group mean.
% b_sigma    - log-residual (scale) dimension-slope draws (draws x KDIM).
% sig_occ    - occasion main-effect SD draws (draws x 1).
% sig_trl    - trial main-effect SD draws (draws x 1).
% sig_trlxid - trial x person SD draws (draws x 1).
% sig_occxid - occasion x person SD draws (draws x 1).
% sig_trlxocc- trial x occasion SD draws (draws x 1).
% ndim       - number of dimensions: 1 or 2.
% z1         - vector of standardized values for dimension 1 (the grid axis).
% obs        - n' (trial count) at which to evaluate the coefficients (scalar).
% CI         - credible-interval width in (0,1], e.g., .95.
%
%Optional Inputs
% z2       - vector of standardized values for dimension 2 (required if ndim=2).
% nocc     - n' (occasion count) at which to evaluate the coefficients (scalar);
%            default 1.
% reltype  - 1, 2, or 3 (see above); default 3 (both facets random).
% sig_dp   - person log-residual SD draws gro_sds(:,2). Required only when
%            'marginal' is 1.
% marginal - 0 (default): evaluate at the typical person (delta_p = 0).
%            1: lognormal population-average residual, multiplying the residual
%            SD by exp(sigma_delta_p^2) before forming the coefficient.
%
%Output (struct)
% out.ndim, out.obs, out.nocc, out.reltype, out.marginal, out.ci
% out.z1            - dimension-1 grid axis (1 x nz1)
% out.z2            - dimension-2 grid axis (1 x nz2; only when ndim=2)
% out.G, out.D      - structs with fields pt/ll/ul (global coefficients over the
%                     grid). ndim=1 -> 1 x nz1; ndim=2 -> nz1 x nz2 (rows index
%                     z1, columns index z2).
% out.ICCg, out.ICCd- relative/absolute single-observation coefficient (the same
%                     reltype evaluated at obs=1, nocc=1) over the grid.
% out.sige_pt       - posterior-mean residual SD sigma_pto,e(z) over the grid.

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

if mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_rel_dynrel_trt for more information about inputs'));
end

sig_p       = psyrat_req(varargin,'sig_p','Please input sig_p (between-person SD draws).');
sig_log0    = psyrat_req(varargin,'sig_log0','Please input sig_log0 (Intercept_sigma draws).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (scale dimension slopes).');
sig_occ     = psyrat_req(varargin,'sig_occ','Please input sig_occ (occasion main-effect SD draws).');
sig_trl     = psyrat_req(varargin,'sig_trl','Please input sig_trl (trial main-effect SD draws).');
sig_trlxid  = psyrat_req(varargin,'sig_trlxid','Please input sig_trlxid (trial x person SD draws).');
sig_occxid  = psyrat_req(varargin,'sig_occxid','Please input sig_occxid (occasion x person SD draws).');
sig_trlxocc = psyrat_req(varargin,'sig_trlxocc','Please input sig_trlxocc (trial x occasion SD draws).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1 (dimension-1 grid axis).');
obs         = psyrat_req(varargin,'obs','Please input obs (n'' trial count).');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

z2       = psyrat_opt(varargin,'z2',[]);
nocc     = psyrat_opt(varargin,'nocc',1);
reltype  = psyrat_opt(varargin,'reltype',3);
sig_dp   = psyrat_opt(varargin,'sig_dp',[]);
marginal = psyrat_opt(varargin,'marginal',0);

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end
if ~isscalar(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a positive scalar.\n');
end
if ~isscalar(nocc) || nocc < 1
    error('varargin:nocc','WARNING: nocc must be a positive scalar.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1, 2, or 3.\n');
end
if marginal == 1 && isempty(sig_dp)
    error('varargin:sig_dp',...
        'WARNING: sig_dp (person log-residual SD draws) is required when marginal = 1.\n');
end

%column vectors of draws
sig_p       = sig_p(:);
sig_log0    = sig_log0(:);
sig_occ     = sig_occ(:);
sig_trl     = sig_trl(:);
sig_trlxid  = sig_trlxid(:);
sig_occxid  = sig_occxid(:);
sig_trlxocc = sig_trlxocc(:);
if isempty(sig_dp); sig_dp = zeros(size(sig_log0)); else; sig_dp = sig_dp(:); end

%b_sigma must be draws x KDIM (KDIM = 1 for one dimension; 3 for two). The
%expectedk column-count check below catches a transposed/mis-shaped input.
kdim = size(b_sigma,2);
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if kdim ~= expectedk
    error('varargin:bsigma',...
        'WARNING: b_sigma has %d columns but ndim=%d expects %d.\n',...
        kdim,ndim,expectedk);
end

z1 = z1(:)';
nz1 = numel(z1);
if ndim == 2
    z2 = z2(:)';
    nz2 = numel(z2);
else
    nz2 = 1;
end

%preallocate the coefficient surfaces (rows index z1, columns index z2)
gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);
D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);
ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);

%shared static facet arguments for psyrat_rel_trt (SDs in; squared internally)
facetargs = {'bo',sig_occ,'bt',sig_trl,'txp',sig_trlxid,...
    'oxp',sig_occxid,'txo',sig_trlxocc};

for a = 1:nz1
    for c = 1:nz2
        %dimension design row matching the Stan model design
        if ndim == 1
            xrow = z1(a);
        else
            xrow = [z1(a) z2(c) z1(a)*z2(c)];
        end

        %residual SD at this z (typical person, delta_p = 0)
        sige_z = exp(sig_log0 + b_sigma * xrow(:));

        %optional lognormal population-average over the person scale random
        %effect: marginal residual SD = exp(alpha + b_sigma z) * exp(sigma_dp^2)
        if marginal == 1
            sige_z = sige_z .* exp(sig_dp.^2);
        end

        sige_pt(a,c) = mean(sige_z);

        %two-facet G(z)/D(z) at n' = (obs trials, nocc occasions) for the chosen
        %reltype; psyrat_rel_trt squares all SD inputs internally.
        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',2,'reltype',reltype,'bp',sig_p,facetargs{:},...
            'err',sige_z,'obs',obs,'nocc',nocc,'CI',ciperc);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',1,'reltype',reltype,'bp',sig_p,facetargs{:},...
            'err',sige_z,'obs',obs,'nocc',nocc,'CI',ciperc);
        %single-observation coefficient (one trial, one occasion) for the same
        %reltype, reported as the dynamic two-facet ICC analogue.
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',2,'reltype',reltype,'bp',sig_p,facetargs{:},...
            'err',sige_z,'obs',1,'nocc',1,'CI',ciperc);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = psyrat_rel_trt(...
            'gcoeff',1,'reltype',reltype,'bp',sig_p,facetargs{:},...
            'err',sige_z,'obs',1,'nocc',1,'CI',ciperc);
    end
end

out = struct();
out.ndim = ndim;
out.obs = obs;
out.nocc = nocc;
out.reltype = reltype;
out.marginal = marginal;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2
    out.z2 = z2;
end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.G    = gridops.squeeze(G,ndim);
out.D    = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
