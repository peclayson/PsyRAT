function psyrat_startview(varargin)
%Loads data for viewing and determines whether to load gui for examining
% single occasion or multiple occasion data
%
%psyrat_startview('file','/path/to/SomePsyRATData.psyrat')
%
%
%Required Inputs:
% No inputs are required.
%
%Optional Inputs:
% file - file of data processed using psyrat_computevarcomp. This optional input
%  is used by psyrat_startproc to provide an easy transition from processing
%  to viewing without the user having to re-select a file.
%
%Output:
% No data are outputted to the Matlab command window. However, the user
%  will have the option of saving various figures and plots that will be 
%  created by psyrat_relfigures, which is executed by this gui 

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

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%see if the file for the figures and tables has been specified in
%varargin
if ~isempty(varargin) && (isempty(psyrat_data) && isempty(psyrat_prefs))
    
    %check if data file has been provided
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = varargin{ind+1};
        
        %if file has been provided, load it
        load(file,'-mat');
    end
 
end %if ~isempty(varargin)

%check if psyrat_gui is open. If the user executes psyrat_startproc and skips
%psyrat_start then there will be no gui to close.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%if the file was not specified, prompt the user to indicate where the file
%is located.
if ~exist('file','var') && isempty(psyrat_data)
    [filepart, pathpart] = uigetfile({'*.psyrat',...
        'PsyRAT Toolbox files (*.psyrat)'},'Data');

    if filepart == 0
        %relaunch FIRST, dialog LAST (B20 creation-order lesson): modality
        %binds only the windows that exist when the dialog is created, so a
        %dialog fired before psyrat_start would be covered by the relaunched
        %home screen and escape its own modality
        psyrat_start;
        errordlg({'No PsyRAT results file was selected.';...
            'How to fix: click View Results again and choose a .psyrat file.'},...
            'File Error','modal');
        return;
    end

    fprintf('\n\nLoading Data...\n\n');
    
    %load data
    load(fullfile(pathpart,filepart),'-mat');
   
end

%if psyrat_prefs does not exist, load the default preferences. If this window
%was not gotten to using psyrat_start, psyrat_prefs will need to be defined
if isempty(psyrat_prefs)
    psyrat_prefs = psyrat_defaults;
    psyrat_prefs.ver = psyrat_defineversion;
end

%legacy or externally produced results may lack the version stamp the
%export callbacks read unguarded; label it as unrecorded rather than
%stamping the viewer's own version onto a file it did not produce
if ~isfield(psyrat_data,'ver') || isempty(psyrat_data.ver)
    psyrat_data.ver = '(not recorded)';
end

%flag a non-CmdStan (native-engine) run before viewing, so the user sees that
%the results are approximate / not bit-identical to CmdStan. CmdStan runs and
%legacy files (no engine field) produce no notice, so the default workflow is
%unchanged. psyrat_engine_notice returns '' in those cases.
%uiwait blocks here until the user dismisses the notice. Without it, warndlg is
%non-modal and the view-setup figure built immediately below (psyrat_startview_*)
%is drawn on top, hiding the warning the user is meant to acknowledge.
if isfield(psyrat_data,'rel')
    engine_notice = psyrat_engine_notice(psyrat_data.rel);
    if ~isempty(engine_notice)
        uiwait(warndlg(engine_notice,'Native estimation engine'));
    end
end

%flag estimation-diagnostic problems (failed convergence and/or divergent
%transitions) before viewing, so a GUI-only user who never reads the MATLAB
%console is not shown reliability coefficients from an untrustworthy fit with no
%warning. psyrat_convergence_notice returns '' for clean runs and for legacy
%files lacking the diagnostic fields, so the default workflow is unchanged.
%uiwait blocks here for the same reason as the engine notice above: a non-modal
%warndlg would be hidden behind the view-setup figure built below.
if isfield(psyrat_data,'rel')
    conv_notice = psyrat_convergence_notice(psyrat_data.rel);
    if ~isempty(conv_notice)
        uiwait(warndlg(conv_notice,'Estimation diagnostics'));
    end
end

%create a gui to allow the user to specify what aspects of the data will be
%viewed
%pic between two possible guis: gui for single session data and a gui for
%multiple session data
%route to the viewer for this analysis family. psyrat_analysis_family is the
%single source of truth for the analysis-string -> family mapping (shared with
%the dynrel/splits summary guards and the headless report dispatch). A missing
%analysis field is the legacy one-facet default -> 'sing'; an unrecognized
%family opens no viewer, matching the prior if/elseif chain's fall-through.
if isfield(psyrat_data.rel,'analysis')
    analysis_str = psyrat_data.rel.analysis;
else
    analysis_str = 'ic';
end
switch psyrat_analysis_family(analysis_str)
    case 'sing'
        psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    case 'dod'
        psyrat_startview_dodiff('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    case 'trt'
        psyrat_startview_trt('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    case 'dynrel'
        psyrat_startview_dynrel('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    case 'splits'
        %nonparallel data splits (observed design, no projection/cutoff)
        psyrat_startview_splits('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

end
