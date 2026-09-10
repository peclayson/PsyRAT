function [ll,pt,ul] = psyrat_cutscore_sing(varargin)
%Calculate one-facet cut-score dependability from the Table 2 formula.
%
%[ll,pt,ul] = psyrat_cutscore_sing('gcoeff',1,'bp',sig_u,'wp',sig_e,...
%   'mu',mu,'cut',cutscore,'obs',20,'CI',.95)
%
%Cut-scores are ABSOLUTE-error only. A criterion-referenced (cut-score)
%coefficient is a dependability/absolute decision (Phi(lambda)) by definition
%(Rocha 2026, Table 2, p.7); the absolute error keeps the trial main effect
%sigma_i^2. A relative-error cut-score is not a valid quantity, so the only
%accepted decision type is gcoeff = 1; gcoeff = 2 (relative) is rejected.
%
%The one-facet absolute cut-score is
%   (sigma_p^2 + (mu-C)^2) / (sigma_p^2 + (mu-C)^2 + (sigma_pi,e^2 + sigma_i^2)/n')
%where the total within-person variance sigma_pi,e^2 + sigma_i^2 is supplied
%directly as wp^2 (callers build wp = sqrt(sigma_pi,e^2 + sigma_i^2)).
%
%Inputs
% gcoeff - decision type. Must be 1 (dependability/absolute). gcoeff = 2
%  (relative) is invalid for a cut-score and raises an error.
% bp - between-person standard deviation draws
% wp - total within-person standard deviation draws (sigma_pi,e^2 + sigma_i^2
%  in variance units)
% mu - universe score mean draws
% cut - criterion/cut-score
% obs - number of observations (scalar or [min max])
% CI - credible interval width (0,1]
%
%ARGUMENT CONTRACT -- READ BEFORE COMPARING THIS FUNCTION WITH
%psyrat_cutscore_dep (RC-03). The two take the SAME argument names with
%OPPOSITE meanings:
%
%  psyrat_cutscore_sing  wp = the TOTAL within-person SD (sigma_pi,e^2 +
%                             sigma_i^2 already inside); there is no i input,
%                             because the error term is complete as passed
%  psyrat_cutscore_dep   wp = the RELATIVE residual sigma_pi,e; sigma_i^2 is
%                             supplied separately as i and ADDED to the error
%
%Callers of this function must therefore build wp themselves. The shipped
%wiring does: psyrat_criterionfigures forms wp = sqrt(sig_e.^2 + sig_trl.^2)
%and passes no i, the same total-within construction as psyrat_relsummary's
%local_total_within (which additionally passes i on, because its consumer
%psyrat_rel_sing subtracts sigma_i^2 back out for the relative branch).
%Handing this function a relative residual, or handing psyrat_cutscore_dep the
%total built here, yields a plausible number rather than an error. This is the
%cut-score half of the psyrat_rel_sing / psyrat_ssrel collision (RC-19).
%
%Outputs
% ll - lower credible interval bound
% pt - posterior mean
% ul - upper credible interval bound

% Copyright (C) 2016-2025 Peter E. Clayson
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

if ~isempty(varargin)
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help psyrat_cutscore_sing for more information about inputs'));
    end
    
    gcoeff = psyrat_req(varargin,'gcoeff',...
        'Please input gcoeff as 1 (absolute/dependability).');
    if gcoeff ~= 1
        error('varargin:gcoeff',... %Error code and associated error
            strcat('WARNING: Cut-scores are absolute-error only.\n\n',...
            'A criterion-referenced cut-score is a dependability (absolute)\n',...
            'decision; a relative-error cut-score is not a valid quantity.\n',...
            'The only accepted value is gcoeff = 1. gcoeff = 2 (relative) is\n',...
            'not permitted. See help psyrat_cutscore_sing for more information.\n'));
    end

    bp = psyrat_req(varargin,'bp','Please input bp.');
    wp = psyrat_req(varargin,'wp','Please input wp.');
    mu = psyrat_req(varargin,'mu','Please input mu.');
    cut = psyrat_req(varargin,'cut','Please input cut.');
    obs = psyrat_req(varargin,'obs','Please input obs.');
    
    ciperc = psyrat_req(varargin,'CI','Please input CI.');
    if ciperc > 1 || ciperc < 0
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval should ',...
            'be a value between 0 and 1\n',...
            'A value of ',sprintf(' %2.2f',ciperc),...
            ' is invalid\n',...
            'See help psyrat_cutscore_sing for more information \n'));
    end
end

if ciperc < .5
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

ciedge = (1-ciperc)/2;

bp = bp(:).^2;
wp = wp(:).^2;
mu = mu(:);

ndraw = length(bp);
wp = psyrat_expanddraw(wp,ndraw,'wp');
mu = psyrat_expanddraw(mu,ndraw,'mu');

offsq = (mu - cut).^2;

if length(obs) == 1
    rel = psyrat_cutscore_rel(bp,wp,offsq,obs);

    ll = quantile(rel,ciedge);
    pt = mean(rel);
    ul = quantile(rel,1-ciedge);

elseif length(obs) == 2
    %one column per requested trial count in the range obs(1):obs(2)
    ll = zeros(1,obs(2)-obs(1)+1);
    pt = zeros(1,obs(2)-obs(1)+1);
    ul = zeros(1,obs(2)-obs(1)+1);

    for ii = obs(1):obs(2)
        rel = psyrat_cutscore_rel(bp,wp,offsq,ii);

        %index relative to the start of the requested range so the outputs
        %align to obs(1):obs(2) (no leading zeros when obs(1) > 1; identical
        %to the previous behavior when obs(1) == 1)
        idx = ii - obs(1) + 1;
        ll(idx) = quantile(rel,ciedge);
        pt(idx) = mean(rel);
        ul(idx) = quantile(rel,1-ciedge);
    end
else
    error('varargin:obs',... %Error code and associated error
        strcat('WARNING: obs should be a scalar or 2-element vector.\n'));
end

end

function rel = psyrat_cutscore_rel(bp,wp,offsq,obs)
%Absolute (dependability) cut-score only. wp is the total within-person
%variance sigma_pi,e^2 + sigma_i^2, so the absolute error keeps the trial
%main effect sigma_i^2 (Rocha 2026, Table 2).
err = wp ./ obs;
rel = (bp + offsq) ./ (bp + offsq + err);
end

function out = psyrat_expanddraw(v,ndraw,name)
if length(v) == 1
    out = repmat(v,ndraw,1);
elseif length(v) == ndraw
    out = v;
else
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: Posterior draw sizes do not align for ',name,'.\n'));
end
end
