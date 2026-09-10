function psyrat_startview_trt(varargin)
%Prepares the data from multiple occasions for viewing and lets user 
% specify which tables and figures to present
%
%psyrat_startview_trt('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
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

if ~isfield(psyrat_prefs.view,'plotssrel')
    psyrat_prefs.view.plotssrel = 1;
end
if ~isfield(psyrat_prefs.view,'plotssicc')
    psyrat_prefs.view.plotssicc = 1;
end
if ~isfield(psyrat_prefs.view,'tablessrel')
    psyrat_prefs.view.tablessrel = 1;
end
%test-retest occasion estimand (single-occasion default vs multi-occasion
%composite); back-compat for prefs saved before this option existed
if ~isfield(psyrat_prefs.view,'noccmode') || isempty(psyrat_prefs.view.noccmode)
    psyrat_prefs.view.noccmode = 1;
end
if ~isfield(psyrat_prefs.view,'nocc')
    psyrat_prefs.view.nocc = [];
end
%difference-score coefficient (only used for two-facet difference scores,
%REL.analysis 'trt_diff'); default to dependability
if ~isfield(psyrat_prefs.view,'diffgcoeff') || isempty(psyrat_prefs.view.diffgcoeff)
    psyrat_prefs.view.diffgcoeff = 1;
end

%analysis-unit label fragments (trial vs split). A parallel-splits run measures
%splits, so the trial-labeled controls switch to "split"; psyrat_unitlabels
%defaults to "trial" when no splits flag is present, so non-splits output is
%byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%define parameters for figure position
figwidth = 680;
%figheight allows two extra rows for the test-retest occasion-score control.
%
%RC-25. 850 was too short for the rows this screen accumulates: the bottom
%button row was laid out at y = -23.75 on 'trt'/'trt_sserrvar' (26 of its 50 px
%inside the window) and at y = -68.75 on 'trt_diff', where NONE of it was inside
%and the user could not press Generate at all. The arithmetic is pure fixed
%pixels -- every Position here is an expression in figheight/rowspace/lcol/rcol,
%and fsize only ever sets the 'fontsize' property -- so this is neither font-size
%nor display dependent, contrary to what the entry assumed. 935 puts the row at
%+61.25 and +16.25 respectively.
figheight = 935;

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

%create the gui.
%
%RC-25. The origin y drops 310 -> 100 in the SAME change as the figheight above,
%and it is not optional. The origin is hardcoded, so window top = y + figheight:
%at 310 the old 850-tall window already reached 1160, which overruns the
%1728x1117 display this defect was reported from. Raising figheight without
%lowering the origin would have pushed the top to 1245 and traded a clipped
%button row for a window whose top is off-screen on the smaller display. 100
%puts the top at 1035, which fits both that display and a 1440-tall one.
psyrat_gui= psyrat_newfigure('unit','pix',...
  'position',[400 100 figwidth figheight],...
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

%next row
row = row - (rowspace*1.4);

str = sprintf(['The reliability threshold to use for retaining data\n'...
    'Participants that do not have enough ' L.units ' to meet this reliability\n'...
    'threshold will be recommended for exclusion']);

%Print the text for reliability cutoff with a box for the user to specify
%the input
uicontrol(psyrat_gui,...
    'Style','text',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Reliability Cutoff:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.relcutoff = uicontrol(psyrat_gui,...
    'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',psyrat_prefs.view.relvalue,...
    'Position', [rcol+5 row+6 figwidth/4 25]);  

%error text for invalid reliability cutoff (shown only when needed)
inputs.rel_error = uicontrol(psyrat_gui,'Style','text',...
    'fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','left',...
    'ForegroundColor',[0.70 0.00 0.00],...
    'String','Reliability cutoff must be numeric and between 0 and 1 (inclusive)',...
    'Visible','off',...
    'Position',[lcol row-14 figwidth-lcol 16]);

%next row
row = row - rowspace-5;

str = sprintf(['Choose dependability to use the absolute error variance\n'...
    'Choose generalizability to use the relative error variance']);

%Provide the user with the option to choose between dependability and
%generalizability coefficients
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','G-Theory Coefficient:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.depgen = uicontrol(psyrat_gui,...
    'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Dependability' 'Generalizability'},...
    'Value',psyrat_prefs.view.gcoeff,...
    'Position', [rcol row figwidth/4 25]);  

%next row
row = row - rowspace-10;

str = sprintf(['For internal consistency, use a coefficient of equivalence\n'...
    'For test-retest reliability, use a coefficient of stability\n'...
    'For combined consistency and stability, use equivalence and stability']);

%Provide the user with the option to choose between coefficients of 
%equivalence and coefficients of stability
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Type of Reliability Coefficient:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.equistab = uicontrol(psyrat_gui,...
    'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Equivalence' 'Stability' 'Equivalence + Stability'},...
    'Value',psyrat_prefs.view.reltype,...
    'Position', [rcol row figwidth/4 25]);

%next row
row = row - rowspace-10;

str = sprintf(['Single-occasion: reliability of a score from one occasion\n'...
    '(n''_o = 1), generalizing over the occasion facet - the default\n'...
    'Multi-occasion composite: reliability of a score averaged across\n'...
    'k occasions (n''_o = k). The SEM matches the selected coefficient.\n'...
    'See the scientific audit (sec.14) for the occasion-estimand convention']);

%Provide the user with the option to choose between a single-occasion score
%(the default test-retest estimand) and a multi-occasion composite score
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Test-Retest Occasion Score:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);

inputs.noccmode = uicontrol(psyrat_gui,...
    'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Single-occasion (n''_o = 1)' 'Multi-occasion composite (n''_o = k)'},...
    'Value',psyrat_prefs.view.noccmode,...
    'Position', [rcol row figwidth/4 25]);

%next row
row = row - rowspace-10;

%default the composite size k to the number of observed occasions
nobsocc = 1;
if isfield(psyrat_data.rel,'time') && ~isempty(psyrat_data.rel.time)
    nobsocc = max(1,numel(psyrat_data.rel.time));
end
if ~isempty(psyrat_prefs.view.nocc)
    noccstr = num2str(psyrat_prefs.view.nocc);
else
    noccstr = num2str(nobsocc);
end

str = sprintf(['Number of occasions (k) to average for the multi-occasion\n'...
    'composite score. Defaults to the number of observed occasions (%d).\n'...
    'Only used when the multi-occasion composite score is selected.'],nobsocc);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','# Occasions in Composite (k):',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);

inputs.nocc = uicontrol(psyrat_gui,...
    'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',noccstr,...
    'Position', [rcol+5 row+6 figwidth/4 25]);

%keep the k box enabled only for the multi-occasion composite score
set(inputs.noccmode,'Callback',{@psyrat_toggle_nocc,inputs.nocc});
if psyrat_prefs.view.noccmode == 2
    set(inputs.nocc,'Enable','on');
else
    set(inputs.nocc,'Enable','off');
end

%next row
row = row - rowspace-10;

%For two-facet difference scores, also let the user choose the difference-score
%coefficient (dependability vs generalizability). This is independent of the
%per-condition G-Theory coefficient above and is only shown for 'trt_diff'.
if isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'trt_diff')

    str = sprintf(['Coefficient for the DIFFERENCE score (condition X - Y)\n'...
        'Choose dependability to use the absolute error variance\n'...
        'Choose generalizability to use the relative error variance\n'...
        'This is separate from the per-condition coefficient above']);

    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','Difference-Score Coefficient:',...
        'Tooltip',str,...
        'Position', [lcol row figwidth/4 25]);

    inputs.diffgcoeff = uicontrol(psyrat_gui,...
        'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',{'Dependability' 'Generalizability'},...
        'Value',psyrat_prefs.view.diffgcoeff,...
        'Position', [rcol row figwidth/4 25]);

    %next row
    row = row - rowspace-10;
end


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
str = sprintf(['Display a plot that shows the impact of the number of\n'...
    L.units ' retained for averaging on reliability estimates']);

%reliability with increasing trials
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',['Plot: ' L.Units ' vs Reliability'],...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.plotrel = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.plotrel,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%plot ICCs
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: ICC Estimates',...
    'Tooltip','Display a plot of the intraclass correlation coefficients',...
    'Position', [lcol row figwidth/2 40]);  

inputs.ploticc = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.ploticc,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

str = sprintf(['Display a plot with subject-level reliability estimates\n'...
    'and credible intervals for each group/event']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Subject-Level Reliability',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);
inputs.plotssrel = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.plotssrel,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row
row = row - rowspace;

str = sprintf(['Display a plot with subject-level ICC estimates\n'...
    'and credible intervals for each group/event']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Subject-Level ICC',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);
inputs.plotssicc = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.plotssicc,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row
row = row - rowspace;

str = sprintf(['Display/export a table with subject-level reliability,\n'...
    'ICC, SEM, and include/exclude flags']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Subject-Level Reliability',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);
inputs.tablessrel = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.tablessrel,...
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

inputs.plotbetsd = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.showstddevf,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row
row = row - rowspace;

str = sprintf(['Display a table showing the ' L.numberof ' needed\n',...
    'to obtain the specified reliability threshold and the\n',...
    'reliability point estimate and credible interval for the ' L.unit ' cutoff']);

%reliability cutoff table
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',['Table: ' L.Unit ' Cutoffs at Reliability Threshold'],...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.relcutt = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.inctrltable,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

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

inputs.overallt = uicontrol(psyrat_gui,'Style','checkbox',...
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

inputs.sdt = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.showstddevt,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%RC-24. The subject-level test-retest design (trt_sserrvar, analysis 25) reports
%per-participant reliability, not the group-level quantities the rest of the
%test-retest family shows. Its relsummary branch populates ssrel_table only -- no
%event.rel, event.icc, event.betsd or event.trlcutoff -- so these five outputs
%cannot be produced for it and are turned off and disabled rather than offered
%and then silently skipped. psyrat_relfigures gates them too; this is the half
%the user can see. The subject-level plots/table above stay on.
%
%relsummary.reltype was in that list until RC-26 stored it. The condition below
%keys off rel.analysis, not off that field, so the five stay disabled.
if isfield(psyrat_data.rel,'analysis') && ...
        strcmpi(psyrat_data.rel.analysis,'trt_sserrvar')
    ss_na = ['Not available for subject-level test-retest data: this design '...
        'estimates a per-participant residual, not the group-level '...
        'reliability, ICC and variance terms this output needs.'];
    for h = [inputs.plotrel inputs.ploticc inputs.plotbetsd ...
            inputs.relcutt inputs.overallt]
        set(h,'Value',0,'Enable','off','Tooltip',ss_na);
    end
end

%next row for buttons
row = row - rowspace*1.4;

%The criterion route needs the group-level components psyrat_cutscore_trt
%consumes (sig_id, sig_err, mu, ...). The subject-level design (trt_sserrvar)
%estimates a per-participant residual instead of those, and the two-facet
%difference D-study (trt_diff) stores its crossed components but no mu, so on
%both the route dies on a missing field with nothing it could compute. Suppress
%the button for them, the same kind of gate psyrat_startview_sing applies on
%its side (showCriterionButton there; the sing gate additionally requires an
%estimated mu, a condition with no trt analogue wired here). This is the
%criterion-route half of the RC-24 class: psyrat_relfigures learned a
%trt_sserr branch, this screen's button never got the matching guard.
showCriterionButton = ~(isfield(psyrat_data.rel,'analysis') && ...
    any(strcmpi(psyrat_data.rel.analysis,{'trt_sserrvar','trt_diff'})));

btn_margin = 14;
btn_gap = 8;
btn_height = 50;
nbuttons = 5 + double(showCriterionButton);
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

function psyrat_toggle_nocc(src,~,noccbox)
%enable the composite-size (k) box only for the multi-occasion composite
%score; for the single-occasion score n'_o is fixed at 1
if src.Value == 2
    set(noccbox,'Enable','on');
else
    set(noccbox,'Enable','off');
end

end

function [noccmode, nocc] = psyrat_get_nocc(inputs)
%pull the test-retest occasion-score selection from the gui controls.
%noccmode: 1 = single-occasion (n'_o = 1), 2 = multi-occasion composite.
%nocc: composite size k (positive integer) or [] to use observed occasions.
noccmode = 1;
nocc = [];
if isfield(inputs,'noccmode') && ishghandle(inputs.noccmode)
    noccmode = inputs.noccmode.Value;
end
if isfield(inputs,'nocc') && ishghandle(inputs.nocc)
    val = str2double(inputs.nocc.String);
    %isreal must come FIRST: str2double parses complex literals like '2+3i',
    %isfinite is true for them and val > 0 compares only the real part, so
    %without it mod() is reached with a complex argument and errors outright
    %("Argument must be real"). This function is called only from pushbutton
    %callbacks and nothing in this file is wrapped in try/catch, so that
    %surfaced as an unhandled callback error rather than the silent fallback
    %to [] (observed occasions) that every other invalid entry produces.
    if isreal(val) && isfinite(val) && val > 0 && mod(val,1) == 0
        nocc = val;
    end
end

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
%Open criterion-score-specific setup for TRT outputs.

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('inputs',varargin),1);
if ~isempty(ind)
    inputs = varargin{ind+1};
    
    releval = relcheck(str2double(inputs.relcutoff.String));
    if releval ~= 0
        if isfield(inputs,'rel_error') && ishandle(inputs.rel_error)
            inputs.rel_error.Visible = 'on';
        end
        if ishghandle(inputs.relcutoff)
            uicontrol(inputs.relcutoff);
        end
        return;
    end
    
    psyrat_prefs.view.relvalue = str2double(inputs.relcutoff.String);
    psyrat_prefs.view.gcoeff = inputs.depgen.Value;
    psyrat_prefs.view.reltype = inputs.equistab.Value;
    [psyrat_prefs.view.noccmode, psyrat_prefs.view.nocc] = psyrat_get_nocc(inputs);
    if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
        psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
    end
    psyrat_prefs.view.plotrel = inputs.plotrel.Value;
    psyrat_prefs.view.ploticc = inputs.ploticc.Value;
    psyrat_prefs.view.plotssrel = inputs.plotssrel.Value;
    psyrat_prefs.view.plotssicc = inputs.plotssicc.Value;
    psyrat_prefs.view.tablessrel = inputs.tablessrel.Value;
    psyrat_prefs.view.inctrltable = inputs.relcutt.Value;
    psyrat_prefs.view.overalltable = inputs.overallt.Value;
    psyrat_prefs.view.showstddevt = inputs.sdt.Value;
    psyrat_prefs.view.showstddevf = inputs.plotbetsd.Value;
end

psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

psyrat_startview_criterion('psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,...
    'analysis','trt');

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

%check whether the reliability estimate provided is numeric and between 0
%and 1
releval = relcheck(str2double(inputs.relcutoff.String));

psyrat_prefs.view.relvalue = str2double(inputs.relcutoff.String);
psyrat_prefs.view.gcoeff = inputs.depgen.Value;
psyrat_prefs.view.reltype = inputs.equistab.Value;
[psyrat_prefs.view.noccmode, psyrat_prefs.view.nocc] = psyrat_get_nocc(inputs);
if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
    psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
end
psyrat_prefs.view.plotrel = inputs.plotrel.Value;
psyrat_prefs.view.ploticc = inputs.ploticc.Value;
psyrat_prefs.view.plotssrel = inputs.plotssrel.Value;
psyrat_prefs.view.plotssicc = inputs.plotssicc.Value;
psyrat_prefs.view.tablessrel = inputs.tablessrel.Value;
psyrat_prefs.view.inctrltable = inputs.relcutt.Value;
psyrat_prefs.view.overalltable = inputs.overallt.Value;
psyrat_prefs.view.showstddevt = inputs.sdt.Value;
psyrat_prefs.view.showstddevf = inputs.plotbetsd.Value;

%if the reliability estimate is invalid, show inline feedback and keep the
%current view open.
if releval ~= 0 
    if isfield(inputs,'rel_error') && ishandle(inputs.rel_error)
        inputs.rel_error.Visible = 'on';
    else
        %create error text
        errorstr = {};
        errorstr{end+1} = 'The reliability estimate must be numeric';
        errorstr{end+1} = 'and between 0 and 1 (inclusive)';
        
        %display error prompt
        errordlg(errorstr,'Invalid reliability cutoff','modal'); %B21 stage 1: modal + purpose title
    end
    if ishghandle(inputs.relcutoff)
        uicontrol(inputs.relcutoff);
    end
    return;
end

if isfield(inputs,'rel_error') && ishandle(inputs.rel_error)
    inputs.rel_error.Visible = 'off';
end

%pass inputs from gui to psyrat_relfigures. Two-facet difference-score data
%(REL.analysis 'trt_diff') route through the difference-score figures/tables;
%the difference-score coefficient (diffgcoeff) defaults to dependability unless
%set in psyrat_prefs.view.
if isfield(psyrat_data.rel,'analysis') && ...
        strcmp(psyrat_data.rel.analysis,'trt_diff')
    relfig_analysis = 'trt_diff';
else
    relfig_analysis = 'trt';
end
psyrat_relfigures('psyrat_data',psyrat_data,...
    'psyrat_prefs',psyrat_prefs,...
    'analysis',relfig_analysis);

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
    releval = relcheck(str2double(inputs.relcutoff.String));

    %if the reliability estimate is invalid, keep this view open and show
    %inline feedback.
    if releval ~= 0 
        if isfield(inputs,'rel_error') && ishandle(inputs.rel_error)
            inputs.rel_error.Visible = 'on';
        else
            %create error text
            errorstr = {};
            errorstr{end+1} = 'The reliability estimate must be numeric';
            errorstr{end+1} = 'and between 0 and 1 (inclusive)';

            %display error prompt
            errordlg(errorstr,'Invalid reliability cutoff','modal'); %B21 stage 1: modal + purpose title
        end
        if ishghandle(inputs.relcutoff)
            uicontrol(inputs.relcutoff);
        end
        return;
    end

    if isfield(inputs,'rel_error') && ishandle(inputs.rel_error)
        inputs.rel_error.Visible = 'off';
    end
    
    psyrat_prefs.view.relvalue = str2double(inputs.relcutoff.String);
    psyrat_prefs.view.gcoeff = inputs.depgen.Value;
    psyrat_prefs.view.reltype = inputs.equistab.Value;
    [psyrat_prefs.view.noccmode, psyrat_prefs.view.nocc] = psyrat_get_nocc(inputs);
    if isfield(inputs,'diffgcoeff') && ishghandle(inputs.diffgcoeff)
        psyrat_prefs.view.diffgcoeff = inputs.diffgcoeff.Value;
    end
    psyrat_prefs.view.plotrel = inputs.plotrel.Value;
    psyrat_prefs.view.ploticc = inputs.ploticc.Value;
    psyrat_prefs.view.plotssrel = inputs.plotssrel.Value;
    psyrat_prefs.view.plotssicc = inputs.plotssicc.Value;
    psyrat_prefs.view.tablessrel = inputs.tablessrel.Value;
    psyrat_prefs.view.inctrltable = inputs.relcutt.Value;
    psyrat_prefs.view.overalltable = inputs.overallt.Value;
    psyrat_prefs.view.showstddevt = inputs.sdt.Value;
    psyrat_prefs.view.showstddevf = inputs.plotbetsd.Value;
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    pos = psyrat_gui.Position;
    close(psyrat_gui);
else
    pos=[400 400 550 760];
end

%define list for plotting reliability against number of trials
rellist = {'Lower Limit' 'Point Estimate' 'Upper Limit'};

%define list for central tendency measures
centlist = {'Mean' 'Median'};

%define space between rows and first row location
rowspace = 35;
row = pos(4) - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (pos(3)/1.5);

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
    'Position', [lcol row pos(4)/4 25]);  

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[rcol row pos(4)/4 25]);

%next row
row = row - rowspace*2;

str = sprintf(['Indicate the estimate that should be plotted in the\n',...
    'figure that shows the relationship between the ' L.numberof '\n',...
    'retained for averaging and reliability']);

%which lines should be plotted on relplot
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Line to plot for reliability',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/4 35]);  

str = sprintf(['Lower Limit of Credible Interval\n',...
    'Reliability Point Estimate\n',...
    'Upper Limit of Credible Interval']);

newprefs.plotrelline = uicontrol(psyrat_gui,'Style','listbox',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',rellist,'Min',1,'Max',1,'Value',psyrat_prefs.view.plotrelline,...
    'Tooltip',str,... 
    'Position', [rcol row pos(4)/4 50]);  

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
    'Position', [rcol row+21 pos(4)/4 25]);  

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
    'String',rellist,'Min',1,'Max',1,...
    'Value',psyrat_prefs.view.meascutoff,...
    'Tooltip',str,...
    'Position', [rcol row pos(4)/4 50]);  

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

newprefs.relcentmeas = uicontrol(psyrat_gui,'Style','listbox',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',centlist,'Min',1,'Max',1,...
    'Value',psyrat_prefs.view.relcentmeas,...
    'Position', [rcol row pos(4)/4 40]);  

%Create a back button that will save inputs for preferences
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [pos(4)/8 50 pos(4)/4 40],...
    'Callback',{@psyrat_prefs_save,'psyrat_prefs',psyrat_prefs,'psyrat_data',...
    psyrat_data,'newprefs',newprefs}); 

%Create button that will go back to psyrat_gui without saving
uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [4.4*pos(4)/8 50 pos(4)/4 40],...
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

%execute psyrat_startview_trt with the old preferences
psyrat_startview_trt('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

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
%psyrat_trt_relvtrialsplot's "x = trials(1):trials(2)" -- a colon expression,
%which rejects complex operands. The criterion path's max(2,round(...)) does not
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
psyrat_prefs.view.plotrelline = newprefs.plotrelline.Value;
psyrat_prefs.view.ntrials = ntrials;
psyrat_prefs.view.meascutoff = newprefs.meascutoff.Value;
psyrat_prefs.view.relcentmeas = newprefs.relcentmeas.Value;

%check if psyrat_gui is open
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

%execute psyrat_startview_trt with the new preferences
psyrat_startview_trt('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);

end

function checkout = relcheck(relvalue)
%ensure that the provided reliability estimate is numeric and between 0
%and 1
%
%Input
% relvalue - reliability threshold estimate from psyrat_startview_sing
%
%Output
% checkout
%   0: reliability estimate is numeric and between 0 and 1
%   1: reliability is string
%   2: reliability is not between 0 and 1

%check whether relvalue is numeric and real. ~isreal is needed because
%str2double parses complex literals ('0.8+1i'), isnan is false for them, and
%the range test below compares only the real part -- so such a cutoff was
%accepted as valid and silently behaved as its real part everywhere it is
%used as a threshold (find(llrel >= relcutoff,1) in psyrat_relsummary), while
%printing as that real part under %0.2f. Nothing downstream re-checks it.
if ~isreal(relvalue) || isnan(relvalue)
    checkout = 1;
else
    %check whether relvalue is between 0 and 1 (inclusive)
    if relvalue >= 0 && relvalue <= 1
        checkout = 0;
    else
        checkout = 2;
    end
end

end
