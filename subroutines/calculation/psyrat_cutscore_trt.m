function [ll,pt,ul] = psyrat_cutscore_trt(varargin)
%Calculate cut-score dependability for a Persons x Trials x Occasions design
%
%[ll,pt,ul] = psyrat_cutscore_trt('gcoeff',1,'reltype',1,...
% 'bp',sig_id,'bo',sig_occ,'bt',sig_trl,'txp',sig_trlxid,...
% 'oxp',sig_occxid,'txo',sig_trlxocc,'err',sig_err,'mu',mu,...
% 'cut',cutscore,'obs',[1 50],'nocc',1,'CI',.95)
%
%Cut-scores are ABSOLUTE-error only. A criterion-referenced (cut-score)
%coefficient is a dependability/absolute decision (Phi(lambda)) by definition
%(Rocha 2026, Table 3, p.8); the absolute error keeps the trial main effect
%sigma_i^2 and (for the two-facet design) sigma_o^2 and sigma_oi^2. A
%relative-error cut-score is not a valid quantity, so the only accepted
%decision type is gcoeff = 1; gcoeff = 2 (relative) is rejected.
%
%Inputs
% gcoeff - decision type. Must be 1 (dependability/absolute). gcoeff = 2
%  (relative) is invalid for a cut-score and raises an error.
% reltype - reliability coefficient to calculate
%  1 - equivalence, 2 - stability, 3 - equivalence and stability
% bp - between-person standard deviation components
% bo - between-occasion standard deviation components
% bt - between-trial standard deviation components
% txp - trial x person standard deviation components
% oxp - occasion x person standard deviation components
% txo - trial x occasion standard deviation components
% err - residual standard deviation components
% mu - universe score mean draws
% cut - cut-score
% obs - number of observations. Can either be 1 number or range [min max]
% nocc - number of occasions (default: 1)
% CI - size of the credible interval in decimal format: .95 = 95%
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
            'See help psyrat_cutscore_trt for more information about inputs'));
    end
    
    ind = find(strcmpi('gcoeff',varargin),1);
    if ~isempty(ind)
        gcoeff = varargin{ind+1};
        if gcoeff ~= 1
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: Cut-scores are absolute-error only.\n\n',...
                'A criterion-referenced cut-score is a dependability (absolute)\n',...
                'decision; a relative-error cut-score is not a valid quantity.\n',...
                'The only accepted value is gcoeff = 1. gcoeff = 2 (relative) is\n',...
                'not permitted. See help psyrat_cutscore_trt for more information.\n'));
        end
    else
        error('varargin:gcoeff',... %Error code and associated error
            strcat('WARNING: G-theory estimate not specified \n\n',...
            'Please input gcoeff. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('reltype',varargin),1);
    if ~isempty(ind)
        reltype = varargin{ind+1};
        if ~any(reltype == [1 2 3])
            error('varargin:reltype',... %Error code and associated error
                strcat('WARNING: reltype is invalid. Valid values are 1, 2, or 3\n'));
        end
    else
        error('varargin:reltype',... %Error code and associated error
            strcat('WARNING: reliability coefficient not specified \n\n',...
            'Please input reltype. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person variance not specified \n\n',...
            'Please input bp. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('bo',varargin),1);
    if ~isempty(ind)
        bo = varargin{ind+1};
    else
        error('varargin:bo',... %Error code and associated error
            strcat('WARNING: Between-occasion variance not specified \n\n',...
            'Please input bo. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('bt',varargin),1);
    if ~isempty(ind)
        bt = varargin{ind+1};
    else
        error('varargin:bt',... %Error code and associated error
            strcat('WARNING: Between-trial variance not specified \n\n',...
            'Please input bt. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('txp',varargin),1);
    if ~isempty(ind)
        txp = varargin{ind+1};
    else
        error('varargin:txp',... %Error code and associated error
            strcat('WARNING: trialxperson variance not specified \n\n',...
            'Please input txp. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('oxp',varargin),1);
    if ~isempty(ind)
        oxp = varargin{ind+1};
    else
        error('varargin:oxp',... %Error code and associated error
            strcat('WARNING: occasionxperson variance not specified \n\n',...
            'Please input oxp. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('txo',varargin),1);
    if ~isempty(ind)
        txo = varargin{ind+1};
    else
        error('varargin:txo',... %Error code and associated error
            strcat('WARNING: trialxoccasion variance not specified \n\n',...
            'Please input txo. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('err',varargin),1);
    if ~isempty(ind)
        err = varargin{ind+1};
    else
        error('varargin:err',... %Error code and associated error
            strcat('WARNING: error variance not specified \n\n',...
            'Please input err. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('mu',varargin),1);
    if ~isempty(ind)
        mu = varargin{ind+1};
    else
        error('varargin:mu',... %Error code and associated error
            strcat('WARNING: Universe score means not specified \n\n',...
            'Please input mu. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('cut',varargin),1);
    if ~isempty(ind)
        cut = varargin{ind+1};
    else
        error('varargin:cut',... %Error code and associated error
            strcat('WARNING: Cut-score not specified \n\n',...
            'Please input cut. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1};
    else
        error('varargin:obs',... %Error code and associated error
            strcat('WARNING: Number of observations not specified \n\n',...
            'Please input obs. See help psyrat_cutscore_trt for more information \n'));
    end
    
    ind = find(strcmpi('nocc',varargin),1);
    if ~isempty(ind)
        nocc = varargin{ind+1};
        if length(nocc) ~= 1 || nocc <= 0
            error('varargin:nocc',... %Error code and associated error
                strcat('WARNING: Number of occasions should be a positive scalar\n'));
        end
    else
        nocc = 1;
    end
    
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help psyrat_cutscore_trt for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_cutscore_trt for more information \n'));
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
bo = bo(:).^2;
bt = bt(:).^2;
txp = txp(:).^2;
oxp = oxp(:).^2;
txo = txo(:).^2;
err = err(:).^2;
mu = mu(:);

ndraw = length(bp);
[bo,bt,txp,oxp,txo,err,mu] = psyrat_expanddraws(bo,bt,txp,oxp,txo,err,mu,ndraw);
offsq = (mu - cut) .^ 2;

if length(obs) == 1
    cutrel = psyrat_cutscore_formula(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,...
        err,offsq,obs,nocc);
    
    %pull the lower limit, point estimate, and upper limit
    ll = quantile(cutrel,ciedge);
    pt = mean(cutrel);
    ul = quantile(cutrel,1-ciedge);
    
elseif length(obs) == 2
    %if obs contains two numbers calculate the reliability for the range
    %specified, one column per requested trial count in obs(1):obs(2)
    ll = zeros(1,obs(2)-obs(1)+1);
    pt = zeros(1,obs(2)-obs(1)+1);
    ul = zeros(1,obs(2)-obs(1)+1);

    for ii = obs(1):obs(2)
        cutrel = psyrat_cutscore_formula(gcoeff,reltype,bp,bo,bt,txp,oxp,...
            txo,err,offsq,ii,nocc);

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

function cutrel = psyrat_cutscore_formula(gcoeff,reltype,bp,bo,bt,txp,oxp,...
    txo,err,offsq,obs,nocc)

%Cut-scores are absolute-error only: gcoeff is guaranteed to be 1
%(dependability) by the caller's validation -- psyrat_cutscore_trt rejects
%gcoeff ~= 1 before this local helper is ever reached, so no second guard is
%needed here. The absolute error keeps the trial main effect and (two-facet)
%occasion / trial x occasion terms.

if reltype == 1 %coefficient of equivalence

    %universe score
    uni_scor = bp + (oxp ./ nocc);

    %absolute error
    err_term = (txp ./ obs) + (err ./ (obs*nocc)) +...
        (bt ./ obs) + (txo ./ (obs*nocc));

elseif reltype == 2 %coefficient of stability

    %universe score
    uni_scor = bp + (txp ./ obs);

    %absolute error
    err_term = (oxp ./ nocc) + (err ./ (obs*nocc)) +...
        (bo ./ nocc) + (txo ./ (obs*nocc));

elseif reltype == 3 %coefficient of trial equivalence and stability

    %universe score
    uni_scor = bp;

    %absolute error
    err_term = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc)) +...
        (bt ./ obs) + (bo ./ nocc) + (txo ./ (obs*nocc));

end

cutrel = (uni_scor + offsq) ./ (uni_scor + offsq + err_term);

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
