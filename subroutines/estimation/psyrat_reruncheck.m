function psyrat_reruncheck
%Execute gui to ask user whether to rerun model with more iterations
%
%
%
%Input
% No inputs required
%
%Output
% No output generated. The user's choice will be checked in psyrat_startproc
%

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

figwidth = 500;
figheight = 200;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2.5;

%initialize gui
%CloseRequestFcn: the window's close box answers "Do Not Rerun" (see
%psyrat_rerunclose below). Without it MATLAB's default closereq deleted the
%figure, uiwait returned, and psyrat_computevarcompwarp's guidata call on the
%empty handle errored AFTER the estimation had finished (council G13).
psyrat_gui= psyrat_newfigure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off',...
  'CloseRequestFcn',@psyrat_rerunclose);

movegui(psyrat_gui,'center');

str = {'Chains did not converge';...
    'Would you like to rerun with more iterations?'};

%Write text
uicontrol(psyrat_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 50]);          

%Create a button that will not rerun the model
uicontrol(psyrat_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String',{'Do Not';'Rerun'},...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_setrerun,0,psyrat_gui}); 

%Create button that will rerun the model with more iterations
uicontrol(psyrat_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Rerun',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@psyrat_setrerun,1,psyrat_gui}); 

%display gui
set(psyrat_gui,'Visible','on');

%tag the gui
psyrat_gui.Tag = 'psyrat_gui';

%wait so the user can respond
uiwait(psyrat_gui);


end

function psyrat_setrerun(varargin)
%to rerun or not to rerun the model... that is the question :)
psyrat_gui = varargin{4};
guidata(psyrat_gui,varargin{3});
uiresume(psyrat_gui);

end

function psyrat_rerunclose(varargin)
%the close box is the third way out of this dialog and it means "Do Not
%Rerun": store 0 and release uiwait, exactly as the button does, and leave
%the figure alive for psyrat_computevarcompwarp to close and delete (it does
%so after either button too). Deleting here would hand the caller an empty
%handle again.
psyrat_gui = varargin{1};
guidata(psyrat_gui,0);
uiresume(psyrat_gui);

end
