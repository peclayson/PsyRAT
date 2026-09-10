function [ll,pt,ul] = psyrat_icc(varargin)
%Calculate ICCs using standard deviation draws from CmdStan
%
%[ll,pt,ul] = psyrat_icc('bp',sig_u,'wp',sig_e,'CI',.95)
%
%
%Inputs
% bp - between-person standard deviation draws from CmdStan (sig_u);
%  squared internally
% wp - within-person standard deviation draws from CmdStan (sig_e);
%  squared internally
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Outputs
% ll - lower limit of the credible interval specified by CI for ICC
% pt - the point estimate for ICC
% ul - upper limit of the credible interval specified by CI for ICC

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
        'See help psyrat_dep for more information about inputs'));
    end
    
    %check if bp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1}; 
    else 
        error('varargin:bp',... %Error code and associated error
        strcat('WARNING: Between-person variance not specified \n\n',... 
        'Please input bp. See help psyrat_icc for more information \n'));
    end
    
    %check if wp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp',varargin),1);
    if ~isempty(ind)
        wp = varargin{ind+1}; 
    else 
        error('varargin:wp',... %Error code and associated error
        strcat('WARNING: Within-person variance not specified \n\n',... 
        'Please input wp. See help psyrat_icc for more information \n'));
    end
    
    %check if wp was specified. 
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
                'See help psyrat_icc for more information \n'));
        end
    else 
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',... 
            'Please input CI. See help psyrat_icc for more information \n'));
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

%compute ICC
icc = bp.^2 ./ (bp.^2 + wp.^2);

%calculate lower limit, point estimate, and upper limit
ll = quantile(icc,ciedge);
pt = mean(icc);
ul = quantile(icc,1-ciedge);

end
