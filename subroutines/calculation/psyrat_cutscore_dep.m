function [ll,pt,ul] = psyrat_cutscore_dep(varargin)
%Calculate cut-score-specific reliability for a one-facet design
%
%[ll,pt,ul] = psyrat_cutscore_dep('bp',var_u,'wp',var_e,'mu',mu,...
%  'cut',cutscore,'obs',[1 50],'CI',.95,'est','dep')
%
%LEGACY -- NOT ON ANY SHIPPED PATH (RC-03). This function has zero production
%call sites. The one-facet cut-score the toolbox actually computes is
%psyrat_cutscore_sing, called from psyrat_criterionfigures and
%psyrat_formula_prompt. The only remaining callers of this file are
%tests/TestCalculationFunctions.m and tests/TestCalculationAccuracyOracle.m,
%so agreement between this function and PsyRATAccuracyOracle says nothing about
%any shipped number (recorded at SCIENTIFIC_FORMULA_AUDIT.md, "Oracle agrees
%with Rocha", superseded 2026-07-25). Retained rather than deleted by owner
%decision. Note for any future deletion sweep: the equally dead sibling
%psyrat_dep is the target of eight production error messages that read "See
%help psyrat_dep", so removing the pair means retargeting those first.
%
%Inputs
% bp - between-person standard deviation components from CmdStan
% wp - RELATIVE within-person residual standard deviation draws (sigma_pi,e,
%  WITHOUT the trial main effect). See ARGUMENT CONTRACT below.
% mu - universe score mean draws from CmdStan
% cut - cut-score
% obs - number of observations. Can either be 1 number or range [min max]
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Optional Inputs
% i - trial/item standard deviation components (sigma_i). Default is 0. Added
%  to the error term, because a cut-score is absolute-error by definition
%  (Rocha 2026, Table 2). At the default of 0 the added term vanishes.
% est - type of estimate to be outputted. Cut-scores are absolute-error
%  (dependability) only, so 'dep' is the only accepted value (default 'dep').
%  'gen' (relative error) is not a valid cut-score coefficient and is rejected.
%
%ARGUMENT CONTRACT -- READ BEFORE COMPARING THIS FUNCTION WITH
%psyrat_cutscore_sing (RC-03). The two take the SAME argument names with
%OPPOSITE meanings:
%
%  psyrat_cutscore_dep   wp = the RELATIVE residual sigma_pi,e; sigma_i^2 is
%                             supplied separately as i and ADDED to the error
%  psyrat_cutscore_sing  wp = the TOTAL within-person SD (sigma_pi,e^2 +
%                             sigma_i^2 already inside); there is no i input
%
%Fed the total wp that psyrat_cutscore_sing expects together with that total's
%own i, this function counts sigma_i^2 twice. Nothing errors -- the result is
%simply a smaller, plausible-looking coefficient -- so wp and i must be passed
%as a PAIR from the same source. This is the cut-score half of the
%psyrat_rel_sing / psyrat_ssrel collision documented as RC-19.
%
%Outputs
% ll - lower limit of the credible interval specified by CI
% pt - the point estimate
% ul - upper limit of the credible interval specified by CI

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
    
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help psyrat_cutscore_dep for more information about inputs'));
    end
    
    %check if bp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person variance not specified \n\n',...
            'Please input bp. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if wp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('wp',varargin),1);
    if ~isempty(ind)
        wp = varargin{ind+1};
    else
        error('varargin:wp',... %Error code and associated error
            strcat('WARNING: Within-person variance not specified \n\n',...
            'Please input wp. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if mu was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('mu',varargin),1);
    if ~isempty(ind)
        mu = varargin{ind+1};
    else
        error('varargin:mu',... %Error code and associated error
            strcat('WARNING: Universe score means not specified \n\n',...
            'Please input mu. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if cut was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('cut',varargin),1);
    if ~isempty(ind)
        cut = varargin{ind+1};
    else
        error('varargin:cut',... %Error code and associated error
            strcat('WARNING: Cut-score not specified \n\n',...
            'Please input cut. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if obs was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1};
    else
        error('varargin:obs',... %Error code and associated error
            strcat('WARNING: Number of observations not specified \n\n',...
            'Please input obs. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if CI was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help psyrat_cutscore_dep for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_cutscore_dep for more information \n'));
    end
    
    %check if i was specified.
    %If it is not found, set as default: 0.
    i = psyrat_opt(varargin,'i',0);
    
    %check if est was specified. Cut-scores are absolute-error (dependability)
    %only: a criterion-referenced decision has no relative-error
    %(generalizability) coefficient (Rocha 2026, Table 2). 'dep' is the only
    %accepted value (the default); 'gen' is rejected.
    ind = find(strcmpi('est',varargin),1);
    if ~isempty(ind) && ~strcmpi(varargin{ind+1},'dep')
        error('varargin:est',... %Error code and associated error
            strcat('WARNING: Cut-scores are dependability (absolute error) ',...
            'only.\n\n',...
            'The only accepted value is est = ''dep''. est = ''gen'' ',...
            '(relative) is not a valid cut-score coefficient.\n',...
            'See help psyrat_cutscore_dep for more information \n'));
    end
    
end

%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
if ciperc < .5
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

%calculate edges for CI
ciedge = (1-ciperc)/2;

%inputs are in standard deviation units, but formulas need variances
bp = bp(:).^2;
wp = wp(:).^2;
i = i(:).^2;
mu = mu(:);

%expand scalar vectors to length of posterior draws
ndraw = length(bp);
[wp,mu,i] = psyrat_expanddraws(wp,mu,i,ndraw);

offsq = (mu - cut) .^ 2;

if length(obs) == 1
    
    %cut-scores use absolute error only (dependability): the trial main
    %effect i enters the error term
    err_term = (wp ./ obs) + (i ./ obs);

    cutrel = (bp + offsq) ./ (bp + offsq + err_term);
    
    %pull the lower limit, point estimate, and upper limit
    ll = quantile(cutrel,ciedge);
    pt = mean(cutrel);
    ul = quantile(cutrel,1-ciedge);
    
elseif length(obs) == 2
    
    %if obs contains two numbers calculate reliability for the range
    %specified, one column per requested trial count in obs(1):obs(2)
    ll = zeros(1,obs(2)-obs(1)+1);
    pt = zeros(1,obs(2)-obs(1)+1);
    ul = zeros(1,obs(2)-obs(1)+1);

    for ii = obs(1):obs(2)
        %cut-scores use absolute error only (dependability)
        err_term = (wp ./ ii) + (i ./ ii);

        cutrel = (bp + offsq) ./ (bp + offsq + err_term);

        %index relative to the start of the requested range so the outputs
        %align to obs(1):obs(2) (no leading zeros when obs(1) > 1; identical
        %to the previous behavior when obs(1) == 1)
        idx = ii - obs(1) + 1;

        %pull the lower limit, point estimate, and upper limit
        ll(idx) = quantile(cutrel,ciedge);
        pt(idx) = mean(cutrel);
        ul(idx) = quantile(cutrel,1-ciedge);
    end
    
end

end

function varargout = psyrat_expanddraws(varargin)
%Expand scalars to the posterior draw length and validate vector sizes.

ndraw = varargin{end};
nout = nargin - 1;
varargout = cell(1,nout);

for i = 1:nout
    v = varargin{i};
    if length(v) == 1
        varargout{i} = repmat(v,ndraw,1);
    elseif length(v) == ndraw
        varargout{i} = v;
    else
        error('varargin:drawsize',... %Error code and associated error
            strcat('WARNING: Posterior draw sizes do not align across inputs.\n',...
            'All vectors should be scalar or have matching lengths.\n'));
    end
end

end
