function psyrat_startview_sing(varargin)
%Prepares the data from a single occasion for viewing and lets user 
% specify which tables and figures to present
%
%psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
%
%
%Required Input
% psyrat_data - PsyRAT Toolbox data structure array
% psyrat_prefs - PsyRAT Toolbox preferences structure array
%
%Output
% No variables will be outputted to the Matlab workspace. Based on the
%  inputs from this gui, psyrat_relfigures will be executed to display various
%  figures and tables

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
if ~isfield(psyrat_prefs.view,'tablessrel')
    psyrat_prefs.view.tablessrel = 1;
end

%which type of analysis is this (are we using subject-level reliability?)
if ~isfield(psyrat_data.rel,'analysis') ||... 
        (isfield(psyrat_data.rel,'analysis') && ...
        ~strcmp(psyrat_data.rel.analysis,'ic_sserrvar') && ...
        ~strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar'))
    sserrvar = 0;
elseif isfield(psyrat_data.rel,'analysis') && ...
        (strcmp(psyrat_data.rel.analysis,'ic_sserrvar') || ...
        strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar'))
    sserrvar = 1;
end

%analysis-unit label fragments (trial vs split). A parallel-splits run measures
%splits, so the trial-labeled controls switch to "split"; psyrat_unitlabels
%defaults to "trial" when no splits flag is present, so non-splits output is
%byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

isdiff = isfield(psyrat_data.rel,'analysis') && ...
    strcmp(psyrat_data.rel.analysis,'ic_diff');

%RC-32. Which COEFFICIENT CONTROL this screen shows must follow which
%coefficient the summary actually consumes, not which family the analysis
%belongs to. Both difference analyses are summarized from diffgcoeff --
%psyrat_relsummary.m:1690-1691 ('sing_diff') and :4499 ('sing_diff_sserr') both
%assign relsummary.gcoeff FROM diffgcoeff -- so both need the Difference-Score
%popup. Before this, ic_diff_sserrvar matched sserrvar first and got a control
%labeled "G-Theory Coefficient" writing prefs.view.gcoeff, which that route
%never reads: psyrat_relfigures passes prefs.view.diffgcoeff (:416-419) and
%psyrat_relsummary_sing_diff_sserr takes only diffgcoeff. The visible control
%was inert. Confirmed live 2026-07-29: with it set to Generalizability, every
%rendered output still read "dependability".
%
%A SEPARATE flag, deliberately: isdiff means "this is the ic_diff screen",
%isdiffcoeff means "this screen's coefficient is diffgcoeff". Widening isdiff
%instead would have silently changed its TWO other uses -- the single-event
%"Cell/Event:" guard just below (which ic_diff suppresses on purpose) and the
%diffgcoeff backfill after it. Cited by position rather than line number on
%purpose: this comment block itself moves every line under it.
%
%It would NOT have loosened the criterion gate, and an earlier version of this
%comment claimed it would. That claim was backwards and is recorded here so it
%is not reintroduced: the gate is a conjunction of negations,
%showCriterionButton = ~sserrvar && ~isdiff && hasEstimatedMu, so adding
%labels to isdiff can only drive it FALSE -- suppressing the button further,
%never exposing it. ic_diff_sserrvar is already stopped there by ~sserrvar.
%RC-05's actual hazard is the opposite operation: REMOVING a term (relaxing the
%gate on the assumption that ~isdiff alone catches every difference analysis)
%is what would put mean(bp)-as-mu live. Neither term may be deleted; widening
%is simply a different change with different consequences.
isdiffcoeff = isfield(psyrat_data.rel,'analysis') && ...
    any(strcmp(psyrat_data.rel.analysis,{'ic_diff','ic_diff_sserrvar'}));

show_single_event_label = false;
single_event_label = '';
if ~isdiff && isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'events') && ...
        ~(ischar(psyrat_data.rel.events) && strcmpi(psyrat_data.rel.events,'none'))
    event_labels = cellstr(string(psyrat_data.rel.events(:)));
    if isscalar(event_labels)
        show_single_event_label = true;
        single_event_label = event_labels{1};
    end
end

%RC-32: gated on isdiffcoeff, not isdiff. The backfill must cover every screen
%that READS view.diffgcoeff, and since RC-32 that is both difference analyses,
%not just ic_diff -- the popup below reads it for ic_diff_sserrvar too. Left at
%isdiff, a prefs struct without the field faulted on the ic_diff_sserrvar screen
%only, which is exactly the shape least likely to be hit in testing.
if isdiffcoeff
    if isempty(psyrat_prefs) || ~isstruct(psyrat_prefs) || ...
            ~isfield(psyrat_prefs,'view') || ~isfield(psyrat_prefs.view,'diffgcoeff')
        psyrat_prefs.view.diffgcoeff = 1;
    end
end
if ~isfield(psyrat_prefs,'view') || ~isfield(psyrat_prefs.view,'gcoeff')
    psyrat_prefs.view.gcoeff = 1;
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%define parameters for figure position
figwidth = 680;
figheight = 620;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

%check that a default fsize has been defined
if ~isfield(psyrat_prefs,'guis') || ~isfield(psyrat_prefs.guis,'fsize')
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end

%create the gui
psyrat_gui= psyrat_newfigure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','View Results Setup',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' psyrat_data.rel.filename],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*.45);

%Print the name of the measurement analyzed
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' psyrat_data.proc.measheader],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]); 

if show_single_event_label
    %next row
    row = row - (rowspace*.55);

    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','center',...
        'String',['Cell/Event:  ' single_event_label],...
        'Tooltip','Single-event reliability cell used for this view',...
        'Position',[0 row figwidth 25]);

    %next row
    row = row - (rowspace*.85);
else
    %next row
    row = row - (rowspace*1.4);
end

str = sprintf(['The reliability threshold to use for retaining data\n'...
    'Participants that do not have enough ' L.units ' to meet this reliability\n'...
    'threshold will be recommended for exclusion']);

%Print the text for reliability cutoff with a box for the user to specify
%the input
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Reliability Cutoff:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.h(1) = uicontrol(psyrat_gui,'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.view.depvalue,...
    'Position', [rcol+5 row+6 figwidth/4 25]);  

%error text for invalid reliability cutoff (shown only when needed)
inputs.dep_error = uicontrol(psyrat_gui,'Style','text',...
    'fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','left',...
    'ForegroundColor',[0.70 0.00 0.00],...
    'String','Reliability cutoff must be numeric and between 0 and 1 (inclusive)',...
    'Visible','off',...
    'Position',[lcol row-14 figwidth-lcol 16]);

if ~sserrvar && ~isdiffcoeff
    row = row - rowspace;
    str = sprintf(['Choose dependability to use the absolute error variance\n'...
        'Choose generalizability to use the relative error variance']);
    
    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','G-Theory Coefficient:',...
        'Tooltip',str,...
        'Position',[lcol row figwidth/3 25]);
    
    inputs.gcoeff = uicontrol(psyrat_gui,...
        'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',{'Dependability' 'Generalizability'},...
        'Value',psyrat_prefs.view.gcoeff,...
        'Position',[rcol+5 row figwidth/4 25]);
elseif sserrvar && ~isdiffcoeff
    %subject-level analyses expose the same dependability (absolute phi_s) vs
    %generalizability (relative G_s) choice as the population one-facet view;
    %psyrat_ssrel uses gcoeff to select phi_s (1) or G_s (2). RC-32: guarded on
    %~isdiffcoeff because ic_diff_sserrvar also sets sserrvar and used to match
    %HERE, taking a gcoeff control its summary route never reads.
    row = row - rowspace;
    str = sprintf(['Choose dependability for the absolute coefficient (phi_s)\n'...
        'Choose generalizability for the relative coefficient (G_s)']);

    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','G-Theory Coefficient:',...
        'Tooltip',str,...
        'Position',[lcol row figwidth/3 25]);

    inputs.gcoeff = uicontrol(psyrat_gui,...
        'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',{'Dependability' 'Generalizability'},...
        'Value',psyrat_prefs.view.gcoeff,...
        'Position',[rcol+5 row figwidth/4 25]);
elseif isdiffcoeff
    %RC-32: both difference analyses land here now. On ic_diff_sserrvar this
    %one control also governs the per-event subject-level rows, because
    %psyrat_relsummary_sing_diff_sserr passes diffgcoeff to psyrat_ssrel as
    %well (:4569) -- the same scope the label already under-describes on
    %ic_diff, kept identical rather than split (owner ruling, 2026-07-29).
    row = row - rowspace;
    str = sprintf(['Choose dependability to use absolute error for the\n',...
        'difference score and generalizability for relative error']);
    
    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','Difference-Score Coefficient:',...
        'Tooltip',str,...
        'Position',[lcol row figwidth/3 25]);
    
    inputs.diffgcoeff = uicontrol(psyrat_gui,...
        'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',{'Dependability' 'Generalizability'},...
        'Value',psyrat_prefs.view.diffgcoeff,...
        'Position',[rcol+5 row figwidth/4 25]);
end

%next row
row = row - rowspace-10;

%indicate that a checked box means yes
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Output Selection',...
    'Position', [rcol+5 row figwidth/4 24]);  

uicontrol(psyrat_gui,'Style','text',...
    'fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','center',...
    'String','(checked = include)',...
    'Position',[rcol+5 row-14 figwidth/4 20]);

%next row
row = row - (rowspace+14);

%increase distance between rows as some descriptions take up more than one
%line
rowspace = 50;
rcol = (figwidth/4)*3;

chckstr = 'Checked = Yes; Unchecked = No';

if ~sserrvar
    str_label = ['Plot: ' L.Units ' vs Reliability'];
    str_tool = sprintf(['Display a plot that shows the impact of the number of\n'...
        L.units ' retained for averaging on reliability estimates']);
elseif sserrvar
    str_label = 'Plot: Subject-Level Reliability';
    str_tool = sprintf(['Display a plot that shows each participant''s \n'...
        'reliability estimate and their credible interval']);
end

%dependability with increasing trials
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',str_label,...
    'Tooltip',str_tool,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(2) = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.plotdep,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 


%next row
row = row - rowspace;

if ~sserrvar
    str_label = 'Plot: ICC Estimates';
    str_tool = 'Display a plot of the intraclass correlation coefficients';
elseif sserrvar
    str_label = 'Plot: Subject-Level ICC Estimates';
    str_tool = sprintf(['Display a plot that shows subject-level \n'...
        'intraclass correlation coefficients']);
end

%plot ICCs
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',str_label,...
    'Tooltip',str_tool,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(3) = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.ploticc,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 


%next row
row = row - rowspace;

%plot between-person standard deviations
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Between-Person Standard Deviations',...
    'Tooltip','Display a plot showing the between-person standard deviations',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(7) = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.showstddevf,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);


%next row
row = row - rowspace;

if ~sserrvar
    str = sprintf(['Display a table showing the ' L.numberof ' needed\n',...
        'to obtain the specified reliability threshold and the\n',...
        'reliability point estimate and credible interval for the ' L.unit ' cutoff']);
    str_label = ['Table: ' L.Unit ' Cutoffs at Reliability Threshold'];
else
    str = sprintf(['Display/export a table with subject-level reliability,\n',...
        'ICC, SEM, and participant include/exclude flags']);
    str_label = 'Table: Subject-Level Reliability';
end

%dependability cutoff table
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',str_label,...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);

if ~sserrvar
    inputs.h(4) = uicontrol(psyrat_gui,'Style','checkbox',...
        'Value',psyrat_prefs.view.inctrltable,...
        'Tooltip',chckstr,...
        'Position', [rcol row+20 figwidth/2 25]);
else
    
    inputs.h(4) = uicontrol(psyrat_gui,'Style','checkbox',...
        'Value',psyrat_prefs.view.tablessrel,...
        'Tooltip',chckstr,...
        'Position', [rcol row+20 figwidth/2 25]);
    
end

%next row
row = row - rowspace;

str = sprintf(['Display a table that summarizes the number of participants\n'...
    'with data that satisfy the reliability threshold, the number of\n'...
    'participants without data the satisfy the reliability threshold\n',...
    'the overall reliability point estimate and credible interval,\n'...
    'and the ' L.unit ' summary information (min, mean, median, max)']);


%overall reliability table
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Overall Reliability Summary',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(5) = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.overalltable,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 


%next row
row = row - rowspace;

%between- and within-person standard deviation tables
str = sprintf(['Display a table that shows the between- and within-person\n',...
    'standard deviations, standard errors of measurement, and intraclass\n',...
    'correlation coefficients']);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Sources of Variance',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(6) = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.showstddevt,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);


%next row for buttons
row = row - rowspace*1.4;

%The criterion coefficient is built from (mu - cut)^2, so it is only meaningful
%when the universe-score mean was actually estimated. The DoD -> single-event
%view can reach here without one: the gamma DoD model parameterizes a per-cell
%LOG-mean rather than an observed-scale mean, and result structs saved before
%the mean was extracted carry none either. Those cases are flagged by the
%producer as mu_source ~= 'estimated'; suppress the button for them rather than
%offering a criterion coefficient computed against a placeholder mean. The flag
%is absent for every ordinary run, so this is a no-op outside the DoD route.
hasEstimatedMu = true;
if isfield(psyrat_data.rel,'out') && isstruct(psyrat_data.rel.out) && ...
        isfield(psyrat_data.rel.out,'mu_source') && ...
        ~isempty(psyrat_data.rel.out.mu_source)
    hasEstimatedMu = strcmpi(char(string(psyrat_data.rel.out.mu_source)),'estimated');
end

showCriterionButton = ~sserrvar && ~isdiff && hasEstimatedMu;
nbuttons = 5 + double(showCriterionButton);
btn_margin = 14;
btn_gap = 8;
btn_height = 50;
btn_width = floor((figwidth - (2*btn_margin) - ((nbuttons-1)*btn_gap)) / nbuttons);
btn_x = btn_margin;

%Create a back button that will take the user back to psyrat_start
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back to Home',...
    'Position', [btn_x row btn_width btn_height],...
    'Callback',{@psyrat_svb,psyrat_gui});
btn_x = btn_x + btn_width + btn_gap;

%Create a button to let the user load a different .psyrat file
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Open Another';'.psyrat File'},...
    'Position', [btn_x row btn_width btn_height],...
    'Tooltip','Open a different processed PsyRAT (.psyrat) file',...
    'Callback',{@psyrat_view_loadnewfile,psyrat_gui});
btn_x = btn_x + btn_width + btn_gap;

%Create button that will check the inputs and begin processing the data
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Generate';'Figures/Tables'},...
    'Position', [btn_x row btn_width btn_height],...
    'Tooltip','Generate the selected plots and tables',...
    'Callback',{@psyrat_svh,'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inputs',inputs}); 
btn_x = btn_x + btn_width + btn_gap;

%Create button that will display preferences
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Output Prefs',...
    'Position', [btn_x row btn_width btn_height],...
    'Tooltip','Set plotting and summary preferences',...
    'Callback',{@psyrat_viewprefs,'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inputs',inputs}); 
btn_x = btn_x + btn_width + btn_gap;

if showCriterionButton
    uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','center',...
        'String',{'Criterion Score';'Outputs'},...
        'Position', [btn_x row btn_width btn_height],...
        'Tooltip','Open criterion-score-specific setup and outputs',...
        'Callback',{@psyrat_viewcriterion,'psyrat_prefs',psyrat_prefs,...
        'psyrat_data',psyrat_data,'inputs',inputs});
    btn_x = btn_x + btn_width + btn_gap;
end

%Create button that will close any open figures other than this gui
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Close PsyRAT';'Outputs'},...
    'Position', [btn_x row btn_width btn_height],...
    'Tooltip','Close only PsyRAT-generated output figures and tables',...
    'Callback',{@psyrat_closefigs}); 

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

end

function psyrat_svb(varargin)
%back button. takes user back to psyrat_start

%close gui
close(varargin{3});

%go back to psyrat_start
psyrat_start;

end

function psyrat_view_loadnewfile(varargin)
%give user the chance to load a new file

%close gui
close(varargin{3});
psyrat_startview;

end

function psyrat_viewcriterion(varargin)
%Open criterion-score-specific setup for one-facet outputs.

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('inputs',varargin),1);
if ~isempty(ind)
    inputs = varargin{ind+1};
    
    depeval = depcheck(str2double(inputs.h(1).String));
    if depeval ~= 0
        if isfield(inputs,'dep_error') && ishandle(inputs.dep_error)
            inputs.dep_error.Visible = 'on';
        end
        if ishghandle(inputs.h(1))
            uicontrol(inputs.h(1));
        end
        return;
    end
    
    psyrat_prefs.view.depvalue = str2double(inputs.h(1).String);
    psyrat_prefs.view.plotdep = inputs.h(2).Value;
    psyrat_prefs.view.ploticc = inputs.h(3).Value;
    is_ss = isfield(psyrat_data.rel,'analysis') && ...
        any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}));
    if is_ss
        psyrat_prefs.view.tablessrel = inputs.h(4).Value;
    else
        psyrat_prefs.view.inctrltable = inputs.h(4).Value;
    end
    psyrat_prefs.view.overalltable = inputs.h(5).Value;
    psyrat_prefs.view.showstddevt = inputs.h(6).Value;
    psyrat_prefs.view.showstddevf = inputs.h(7).Value;
    if isfield(inputs,'gcoeff') && ishghandle(inputs.gcoeff)
        psyrat_prefs.view.gcoeff = inputs.gcoeff.Value;
    end
    if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
        psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
    end
end

psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

psyrat_startview_criterion('psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,...
    'analysis','sing');

end

function psyrat_closefigs(varargin)
%button to close PsyRAT-generated output windows while leaving other figures
%alone.
outfigs = findall(0,'Type','figure','Tag','psyrat_output');
if ~isempty(outfigs)
    close(outfigs);
end


end

function psyrat_svh(varargin)
%parses inputs to psyrat_relfigures for displaying figures
%
%Input
% psyrat_data - PsyRAT Toolbox data structure array
% psyrat_prefs - PsyRAT Toolbox preferences structure array
%

%somersault through varargin inputs to check for psyrat_prefs and psyrat_data
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%find inputs
ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

%subject-level analysis mode for this callback
is_ss = isfield(psyrat_data.rel,'analysis') && ...
    any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}));

%check whether the reliability estimate provided is numeric and between 0
%and 1
depeval = depcheck(str2double(inputs.h(1).String));

psyrat_prefs.view.depvalue = str2double(inputs.h(1).String);
psyrat_prefs.view.plotdep = inputs.h(2).Value;
psyrat_prefs.view.ploticc = inputs.h(3).Value;
if is_ss
    psyrat_prefs.view.tablessrel = inputs.h(4).Value;
else
    psyrat_prefs.view.inctrltable = inputs.h(4).Value;
end
psyrat_prefs.view.overalltable = inputs.h(5).Value;
psyrat_prefs.view.showstddevt = inputs.h(6).Value;
psyrat_prefs.view.showstddevf = inputs.h(7).Value;
if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
    psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
end
if isfield(inputs,'gcoeff') && ishghandle(inputs.gcoeff)
    psyrat_prefs.view.gcoeff = inputs.gcoeff.Value;
end

%if the reliability estimate is invalid, show inline feedback and keep the
%current view open.
if depeval ~= 0 
    if isfield(inputs,'dep_error') && ishandle(inputs.dep_error)
        inputs.dep_error.Visible = 'on';
    else
        %create error text
        errorstr = {};
        errorstr{end+1} = 'The reliability estimate must be numeric';
        errorstr{end+1} = 'and between 0 and 1 (inclusive)';
        
        %display error prompt
        errordlg(errorstr,'Invalid reliability cutoff','modal'); %B21 stage 1: modal + purpose title
    end
    if ishghandle(inputs.h(1))
        uicontrol(inputs.h(1));
    end
    return;
end

if isfield(inputs,'dep_error') && ishandle(inputs.dep_error)
    inputs.dep_error.Visible = 'off';
end

%which type of analysis is this (are we using subject-level reliability?)
if ~isfield(psyrat_data.rel,'analysis') ||... 
        (isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'ic'))
    analysis_type = 'sing';
elseif isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'ic_sserrvar')
    analysis_type = 'sing_sserr';
elseif isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar')
    analysis_type = 'sing_diff_sserr';
elseif isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'ic_diff')
    analysis_type = 'sing_diff';
end

%pass inputs from gui to psyrat_relfigures
psyrat_relfigures('psyrat_data',psyrat_data,...
    'psyrat_prefs',psyrat_prefs,...
    'analysis',analysis_type);

end


function psyrat_viewprefs(varargin)
%displays various preferences for plotting or summarizing data
%
%Input
% psyrat_data - PsyRAT Toolbox data structure array
% psyrat_prefs - PsyRAT Toolbox preferences structure array

%somersault through varargin inputs to check for psyrat_prefs and psyrat_data
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%analysis-unit label fragments (trial vs split) for the preference prompts;
%defaults to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%find inputs
ind = find(strcmp('inputs',varargin),1);

if ~isempty(ind)
    inputs = varargin{ind+1};

    %check whether the reliability estimate provided is numeric and between 0
    %and 1
    depeval = depcheck(str2double(inputs.h(1).String));

    %if the reliability estimate is invalid, keep this view open and show
    %inline feedback.
    if depeval ~= 0 
        if isfield(inputs,'dep_error') && ishandle(inputs.dep_error)
            inputs.dep_error.Visible = 'on';
        else
            %create error text
            errorstr = {};
            errorstr{end+1} = 'The reliability estimate must be numeric';
            errorstr{end+1} = 'and between 0 and 1 (inclusive)';
            
            %display error prompt
            errordlg(errorstr,'Invalid reliability cutoff','modal'); %B21 stage 1: modal + purpose title
        end
        if ishghandle(inputs.h(1))
            uicontrol(inputs.h(1));
        end
        return;
    end

    if isfield(inputs,'dep_error') && ishandle(inputs.dep_error)
        inputs.dep_error.Visible = 'off';
    end

    psyrat_prefs.view.depvalue = str2double(inputs.h(1).String);
    psyrat_prefs.view.plotdep = inputs.h(2).Value;
    psyrat_prefs.view.ploticc = inputs.h(3).Value;
    is_ss = isfield(psyrat_data.rel,'analysis') && ...
        any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}));

    if is_ss
        psyrat_prefs.view.tablessrel = inputs.h(4).Value;
    else
        psyrat_prefs.view.inctrltable = inputs.h(4).Value;
    end
    psyrat_prefs.view.overalltable = inputs.h(5).Value;
    psyrat_prefs.view.showstddevt = inputs.h(6).Value;
    psyrat_prefs.view.showstddevf = inputs.h(7).Value;
    if isfield(inputs,'gcoeff') && ishghandle(inputs.gcoeff)
        psyrat_prefs.view.gcoeff = inputs.gcoeff.Value;
    end
    if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
        psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
    end
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    pos=[400 400 550 550];
end

%define list for plotting reliability against number of trials
deplist = {'Lower Limit' 'Point Estimate' 'Upper Limit'};

%define list for central tendency measures
centlist = {'Mean' 'Median'};

%define space between rows and first row location
rowspace = 35;
row = pos(4) - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (pos(3)/2+20);

%create the basic psyrat_prefs
psyrat_gui = psyrat_newfigure('unit','pix',...
  'position',pos,...
  'menub','no',...
  'name','View Preferences',...
  'numbertitle','off',...
  'resize','off');    

%print the gui headers
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [pos(4)/8 row pos(4)/3 25]);  

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[4.4*pos(4)/8 row pos(4)/3 25]);

%next row
row = row - rowspace*2;

str = sprintf(['Indicate the estimate that should be plotted in the\n',...
    'figure that shows the relationship between the ' L.numberof '\n',...
    'retained for averaging and reliability']);

%which lines should be plotted on reliability plot
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Line to plot for reliability',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

str = sprintf(['Lower Limit of Credible Interval\n',...
    'Reliability Point Estimate\n',...
    'Upper Limit of Credible Interval']);

newprefs.plotdepline = uicontrol(psyrat_gui,'Style','listbox',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',deplist,'Min',1,'Max',1,'Value',psyrat_prefs.view.plotdepline,...
    'Tooltip',str,... 
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2;

str = sprintf(['Indicate the ' L.numberof ' that should be plotted in the\n',...
    'figure that shows the relationship between the ' L.numberof '\n',...
    'retained for averaging and reliability']);

%number of trials to plot
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',['Number of ' L.units ' to plot'],...
    'Tooltip',str,...
    'Position', [lcol row+5 pos(4)/2 35]);  

newprefs.ntrials = uicontrol(psyrat_gui,...
    'Style','edit','fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.view.ntrials,... 
    'Position', [rcol row+21 pos(4)/3 25]);  

%next row
row = row - rowspace*2.2;

str = sprintf(['Indicate which reliability estimate should be used\n',...
    'for the reliability threshold that deems data as reliable']);

%how to determine cutoff
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',['Estimate to use for ' L.unit ' cutoffs'],...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

str = sprintf(['Lower Limit of Credible Interval\n',...
    'Reliability Point Estimate\n',...
    'Upper Limit of Credible Interval']);

newprefs.meascutoff = uicontrol(psyrat_gui,...
    'Style','listbox','fontsize',psyrat_prefs.guis.fsize,...
    'String',deplist,'Min',1,'Max',1,...
    'Value',psyrat_prefs.view.meascutoff,...
    'Tooltip',str,...
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2.2;

str = sprintf(['Indicate the measure of central tendency to use for the\n',...
    'estimation of the overall reliability of the dataset']);

%measure of central tendendcy for overall reliability
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Measure of central tendency for overall reliability calculations',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.depcentmeas = uicontrol(psyrat_gui,'Style','listbox',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',centlist,'Min',1,'Max',1,...
    'Value',psyrat_prefs.view.depcentmeas,...
    'Position', [rcol row pos(4)/3 40]);  

%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will save inputs for preferences
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [pos(4)/8 row pos(4)/3 40],...
    'Callback',{@psyrat_prefs_save,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data,'newprefs',newprefs}); 

%Create button that will go back to psyrat_gui without saving
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [4.4*pos(4)/8 row pos(4)/3 40],...
    'Callback',{@psyrat_prefs_back,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data});

%tag gui
psyrat_gui.Tag = 'psyrat_gui';

end

function psyrat_prefs_back(varargin)
%if the back button was pressed the inputs will not be saved

%somersault through varargin inputs to check for psyrat_prefs and psyrat_data
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%execute psyrat_startview_sing with the old preferences
psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function psyrat_prefs_save(varargin)
%if the save button was pressed use new inputs

%somersault through varargin inputs to check for psyrat_prefs and psyrat_data
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%analysis-unit label fragments (trial vs split) for the validation message;
%defaults to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%find newprefs
ind = find(strcmp('newprefs',varargin),1);
newprefs = varargin{ind+1};

%check the number of trials before leaving this gui.
%~isreal rejects complex literals like '50-4i': isnan is false for them and
%ntrials <= 1 compares only the real part, so one reached prefs.view.ntrials and
%then psyrat_relfigures, which passes it as [1 plotntrials] into
%psyrat_depvtrialsplot's "x = trials(1):trials(2)" -- a colon expression, which
%rejects complex operands. The criterion path's max(2,round(...)) does not
%sanitize it: round acts on both parts and max compares by magnitude, so a
%complex value with modulus above 2 survives both.
ntrials = str2double(newprefs.ntrials.String);
if ~isreal(ntrials) || isnan(ntrials) || ntrials <= 1
    %B21 stage 1 (2026-08-16): modal + purpose title (see the ledger's B21
    %record; the single-argument form was invisible behind the parent figure)
    errordlg(['Please specify at least two ' L.units ' to plot ', ...
        'for reliability estimates'],...
        'Invalid plot trial count','modal');
    if ishghandle(newprefs.ntrials)
        uicontrol(newprefs.ntrials);
    end
    return;
end

%pull new preferences
psyrat_prefs.view.plotdepline = newprefs.plotdepline.Value;
psyrat_prefs.view.ntrials = ntrials;
psyrat_prefs.view.meascutoff = newprefs.meascutoff.Value;
psyrat_prefs.view.depcentmeas = newprefs.depcentmeas.Value;

%check if psyrat_gui is open
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%execute psyrat_startview_sing with the new preferences
psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function checkout = depcheck(depvalue)
%ensure that the provided dependability estimate is numeric and between 0
%and 1
%
%Input
% depvalue - dependability threshold estimate from psyrat_startview_sing
%
%Output
% checkout
%   0: dependability estimate is numeric and between 0 and 1
%   1: dependability is string
%   2: dependability is not between 0 and 1

%check whether depvalue is numeric and real. ~isreal is needed because
%str2double parses complex literals ('0.8+1i'), isnan is false for them, and
%the range test below compares only the real part -- so such a cutoff was
%accepted as valid and silently behaved as its real part everywhere it is
%used as a threshold (find(llrel >= depcutoff,1) and the subject-level
%ind2include/ind2exclude split in psyrat_relsummary), while printing as that
%real part under %0.2f. Nothing downstream re-checks it.
if ~isreal(depvalue) || isnan(depvalue)
    checkout = 1;
else
    %check whether depvalue is between 0 and 1 (inclusive)
    if depvalue >= 0 && depvalue <= 1
        checkout = 0;
    else
        checkout = 2;
    end
end

end
