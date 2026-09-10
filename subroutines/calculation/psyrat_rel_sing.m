function [ll,pt,ul] = psyrat_rel_sing(varargin)
%Calculate one-facet reliability/ICC formulas from Table 2.
%
%[ll,pt,ul] = psyrat_rel_sing('gcoeff',1,'metric','global',...
%   'bp',sig_u,'wp',sig_e,'obs',20,'CI',.95)
%
%Inputs
% gcoeff - decision type: 1 = dependability/absolute, 2 = generalizability/relative
% metric - 'global' or 'icc'
% bp - between-person standard deviation draws
% wp - total within-person standard deviation draws
% obs - number of observations (required for metric='global')
% CI - credible interval width (0,1]
%
%Optional Inputs
% i - item/trial standard deviation draws. Used to separate
%  sigma_i^2 from sigma_pi,e^2 for relative formulas. Default: 0.
%
%ARGUMENT CONTRACT -- READ BEFORE COMPARING THIS FUNCTION WITH psyrat_ssrel
%(RC-19). The two take the SAME argument names with OPPOSITE meanings:
%
%  psyrat_rel_sing  wp = the TOTAL within-person SD; i^2 is SUBTRACTED to
%                        recover the relative residual (pierr = wp^2 - i^2)
%  psyrat_ssrel     wp = the RELATIVE residual sigma_pi,e^2(s); i^2 is ADDED
%                        for the absolute branch
%
%wp and i must therefore be passed as a PAIR from the same source. Passing wp
%without its i, or a relative residual as wp, yields a plausible number rather
%than an error -- which is why the pierr < 0 branch below warns instead of
%clipping. psyrat_cutscore_sing and psyrat_cutscore_dep carry the same name
%collision. See psyrat_relsummary's local_total_within for the pairing rule.
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
            'See help psyrat_rel_sing for more information about inputs'));
    end
    
    gcoeff = psyrat_req(varargin,'gcoeff',...
        'Please input gcoeff as 1 (absolute) or 2 (relative).');
    if ~any(gcoeff == [1 2])
        error('varargin:gcoeff',... %Error code and associated error
            strcat('WARNING: gcoeff is invalid. Valid values are 1 or 2\n'));
    end
    
    metric = psyrat_req(varargin,'metric',...
        'Please input metric as ''global'' or ''icc''.');
    metric = lower(char(string(metric)));
    if ~any(strcmp(metric,{'global','icc'}))
        error('varargin:metric',... %Error code and associated error
            strcat('WARNING: metric is invalid. Valid values are ',...
            '''global'' or ''icc''\n'));
    end
    
    bp = psyrat_req(varargin,'bp','Please input bp.');
    wp = psyrat_req(varargin,'wp','Please input wp.');
    
    if strcmp(metric,'global')
        obs = psyrat_req(varargin,'obs','Please input obs for global coefficients.');
    else
        obs = psyrat_opt(varargin,'obs',1);
    end
    
    ciperc = psyrat_req(varargin,'CI','Please input CI.');
    if ciperc > 1 || ciperc < 0
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval should ',...
            'be a value between 0 and 1\n',...
            'A value of ',sprintf(' %2.2f',ciperc),...
            ' is invalid\n',...
            'See help psyrat_rel_sing for more information \n'));
    end
    
    i = psyrat_opt(varargin,'i',0);
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
i = i(:).^2;

ndraw = length(bp);
wp = psyrat_expanddraw(wp,ndraw,'wp');
i = psyrat_expanddraw(i,ndraw,'i');

pierr = wp - i;
if any(pierr < 0)
    %NOT clipped. On every correct call site this cannot happen: the caller
    %passes the TOTAL within-person SD as wp and the trial main effect as i, so
    %pierr = wp^2 - i^2 = sigma_(pi,e)^2, which is non-negative by
    %construction. A negative value therefore means the pair was not passed
    %together -- the contract violation described in psyrat_relsummary's
    %local_total_within -- and silently flooring it at zero would absorb that
    %mistake into a plausible-looking coefficient instead of exposing it.
    warning('varargin:i',... %Warning code and associated warning
        ['WARNING: i^2 exceeded wp^2 for %d of %d draws, so the relative ',...
        'error\nvariance came out negative and the generalizability ',...
        'coefficient below is not\ninterpretable. This means wp was not the ',...
        'TOTAL within-person SD while i was the\ntrial main effect; the two ',...
        'must be passed together. Nothing has been clipped --\nthe value ',...
        'shown is the formula evaluated exactly.\n'],...
        sum(pierr < 0),numel(pierr));
end

if strcmp(metric,'icc')
    if gcoeff == 1
        rel = bp ./ (bp + wp);
    else
        rel = bp ./ (bp + pierr);
    end
    
    ll = quantile(rel,ciedge);
    pt = mean(rel);
    ul = quantile(rel,1-ciedge);
    return;
end

if length(obs) == 1
    rel = psyrat_global_rel(bp,wp,pierr,obs,gcoeff);
    
    ll = quantile(rel,ciedge);
    pt = mean(rel);
    ul = quantile(rel,1-ciedge);
    
elseif length(obs) == 2
    %one column per requested trial count in the range obs(1):obs(2)
    nobs = obs(2)-obs(1)+1;
    ll = zeros(1,nobs);
    pt = zeros(1,nobs);
    ul = zeros(1,nobs);

    for ii = obs(1):obs(2)
        rel = psyrat_global_rel(bp,wp,pierr,ii,gcoeff);
        %index relative to the range start so outputs align to obs(1):obs(2)
        %(identical to the previous behavior when obs(1) == 1)
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

function rel = psyrat_global_rel(bp,wp,pierr,obs,gcoeff)
if gcoeff == 1
    err = wp ./ obs;
else
    err = pierr ./ obs;
end
rel = bp ./ (bp + err);
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
