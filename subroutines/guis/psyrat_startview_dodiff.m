function psyrat_startview_dodiff(varargin)
%Prepare/view task-free 4-event difference-of-differences reliability
%
%psyrat_startview_dodiff('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
%
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
%

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

if ~isfield(psyrat_data,'rel') || ~isfield(psyrat_data.rel,'analysis') || ...
        ~strcmpi(psyrat_data.rel.analysis,'ic_dodiff')
    psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
    return;
end

if ~isfield(psyrat_prefs,'view')
    psyrat_prefs.view = struct;
end

[groupnames,eventnames,trialstats] = local_dod_stats(psyrat_data.rel);
psyrat_prefs = local_dod_initprefs(psyrat_prefs,psyrat_data.rel,...
    groupnames,eventnames,trialstats);

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

if ~isfield(psyrat_prefs,'guis') || ~isfield(psyrat_prefs.guis,'fsize')
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end

figwidth = 820;
figheight = 610;
rowspace = 34;
row = figheight - 60;
lcol = 30;
rcol = 470;

psyrat_gui = psyrat_newfigure('unit','pix',...
    'position',[380 320 figwidth figheight],...
    'menub','no',...
    'name','DoD Reliability Setup',...
    'numbertitle','off',...
    'resize','off');

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' psyrat_data.rel.filename],...
    'Position',[0 row figwidth 24]);
row = row - rowspace;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' psyrat_data.proc.measheader],...
    'Position',[0 row figwidth 24]);
row = row - rowspace*1.2;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Contrast form: (ERP_1 - ERP_2) - (ERP_3 - ERP_4)',...
    'Position',[lcol row figwidth-2*lcol 24]);
row = row - rowspace;

if numel(groupnames) > 1
    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','Group:',...
        'Position',[lcol row 130 24]);

    inputs.group = uicontrol(psyrat_gui,'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',groupnames(:),...
        'Value',psyrat_prefs.view.dodgroup,...
        'Position',[lcol+130 row 220 26]);
else
    inputs.group = [];
    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','Group: none',...
        'Position',[lcol row 260 24]);
end

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Heatmap metric:',...
    'Position',[rcol row 130 24]);
inputs.metric = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',{'Dependability';'Generalizability'},...
    'Value',psyrat_prefs.view.dodheatmetric,...
    'Position',[rcol+130 row 180 26]);
row = row - rowspace;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','ERP_1:',...
    'Position',[lcol row 130 24]);
inputs.erp1 = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',eventnames(:),...
    'Value',psyrat_prefs.view.dodmap(1),...
    'Position',[lcol+130 row 220 26]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Trials for ERP_1:',...
    'Position',[rcol row 130 24]);
inputs.trials1 = uicontrol(psyrat_gui,'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.dodtrials(1)),...
    'Position',[rcol+130 row 90 26]);
row = row - rowspace;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','ERP_2:',...
    'Position',[lcol row 130 24]);
inputs.erp2 = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',eventnames(:),...
    'Value',psyrat_prefs.view.dodmap(2),...
    'Position',[lcol+130 row 220 26]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Trials for ERP_2:',...
    'Position',[rcol row 130 24]);
inputs.trials2 = uicontrol(psyrat_gui,'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.dodtrials(2)),...
    'Position',[rcol+130 row 90 26]);
row = row - rowspace;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','ERP_3:',...
    'Position',[lcol row 130 24]);
inputs.erp3 = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',eventnames(:),...
    'Value',psyrat_prefs.view.dodmap(3),...
    'Position',[lcol+130 row 220 26]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Trials for ERP_3:',...
    'Position',[rcol row 130 24]);
inputs.trials3 = uicontrol(psyrat_gui,'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.dodtrials(3)),...
    'Position',[rcol+130 row 90 26]);
row = row - rowspace;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','ERP_4:',...
    'Position',[lcol row 130 24]);
inputs.erp4 = uicontrol(psyrat_gui,'Style','pop',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',eventnames(:),...
    'Value',psyrat_prefs.view.dodmap(4),...
    'Position',[lcol+130 row 220 26]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Trials for ERP_4:',...
    'Position',[rcol row 130 24]);
inputs.trials4 = uicontrol(psyrat_gui,'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.dodtrials(4)),...
    'Position',[rcol+130 row 90 26]);
row = row - rowspace*1.2;

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Heatmap trial grid min:',...
    'Position',[lcol row 190 24]);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','0',...
    'Position',[lcol+190 row 40 24]);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Heatmap trial grid max:',...
    'Position',[rcol row 190 24]);
inputs.gridmax = uicontrol(psyrat_gui,'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.dodgridmax),...
    'Position',[rcol+190 row 90 26]);
row = row - rowspace*1.1;

uicontrol(psyrat_gui,'Style','text','fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','left',...
    'String','Heatmaps use posterior mean variance components for speed. Summary table uses posterior draws.',...
    'Position',[lcol row figwidth-(2*lcol) 22]);
row = row - rowspace*1.8;

btn_margin = 12;
btn_gap = 6;
btn_height = 50;
btn_count = 8;
btn_width = floor((figwidth - (2*btn_margin) - ((btn_count-1)*btn_gap)) / btn_count);
btn_x = btn_margin;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back to Home',...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_back,psyrat_gui});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Open Another';'.psyrat File'},...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_opennew,psyrat_gui});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Generate DoD';'Outputs'},...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_generate,...
    'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inputs',inputs,'groupnames',groupnames,'eventnames',eventnames});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','center',...
    'String',{'Other Single-';'Session Metrics...'},...
    'TooltipString',['Open a compact window to route selected group DoD ',...
    'results into standard single-event or 2-event difference outputs.'],...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_othermetrics,...
    'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inputs',inputs,'groupnames',groupnames,'eventnames',eventnames});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','center',...
    'String',{'Observed';'DoD Table'},...
    'TooltipString',['Build observed DoD reliability table using the mean ',...
    'recorded trial counts for the selected group and contrast.'],...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_observedtable,...
    'psyrat_data',psyrat_data,'inputs',inputs,'eventnames',eventnames});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','center',...
    'String',{'Export DoD';'Var Components'},...
    'TooltipString',['Export DoD contrast variance components (point + CrI) ',...
    'using observed mean trials by default.'],...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_exportvarcomp,...
    'psyrat_data',psyrat_data,'inputs',inputs,'eventnames',eventnames});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Refresh',...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_refresh,...
    'psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data,...
    'inputs',inputs,'groupnames',groupnames,'eventnames',eventnames});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Close PsyRAT';'Outputs'},...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_dod_closefigs});

psyrat_gui.Tag = 'psyrat_gui';

end

function psyrat_dod_back(varargin)
close(varargin{3});
psyrat_start;
end

function psyrat_dod_opennew(varargin)
close(varargin{3});
psyrat_startview;
end

function psyrat_dod_refresh(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
%Read the CURRENT selections before rebuilding (pre-beta batch, 2026-09-05,
%owner ruling). This used to rebuild from the preferences captured when the
%button was built and read none of the controls, so a group change followed
%by Refresh came back on the original group, while chapter 9 says Refresh is
%how the trial defaults follow a group change. The reader is Generate's
%(local_read_main_inputs); the trial fields are NOT required here (false),
%and the stored trial counts are cleared so local_dod_initprefs recomputes
%the selected group's defaults. The ERP mapping and the metric are kept as
%selected, and the typed grid max is kept when it parses. If the selections
%cannot be read (an invalid mapping), the rebuild falls back to the captured
%preferences exactly as before; Generate reports the problem when pressed.
ind = find(strcmp('inputs',varargin),1);
if ~isempty(ind)
    inputs = varargin{ind+1};
    ind = find(strcmp('eventnames',varargin),1);
    eventnames = varargin{ind+1};
    [ok,state] = local_read_main_inputs(inputs,eventnames,false);
    if ok
        if ~isfield(psyrat_prefs,'view'); psyrat_prefs.view = struct; end
        psyrat_prefs.view.dodgroup = state.gidx;
        psyrat_prefs.view.dodmap = state.map;
        psyrat_prefs.view.dodheatmetric = state.metric;
        psyrat_prefs.view.dodtrials = [];   %recomputed for the selected group
        gridmax = local_nonnegint(inputs.gridmax.String);
        if ~isnan(gridmax)
            psyrat_prefs.view.dodgridmax = gridmax;
        end
    end
end
psyrat_startview_dodiff('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

function psyrat_dod_generate(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

ind = find(strcmp('groupnames',varargin),1);
groupnames = varargin{ind+1};

ind = find(strcmp('eventnames',varargin),1);
eventnames = varargin{ind+1};

[ok,state,errmsg,errtitle] = local_read_main_inputs(inputs,eventnames,true);
if ~ok
    errordlg(errmsg,errtitle,'modal');
    return;
end

if ~isfield(psyrat_prefs,'view')
    psyrat_prefs.view = struct;
end
psyrat_prefs.view.dodgroup = state.gidx;
psyrat_prefs.view.dodmap = state.map;
psyrat_prefs.view.dodtrials = state.obs;
psyrat_prefs.view.dodgridmax = state.gridmax;
psyrat_prefs.view.dodheatmetric = state.metric;

[bp,bt,erlog] = local_groupdraws(psyrat_data.rel,state.gidx);
bp = bp(:,state.map,state.map);
bt = bt(:,state.map,state.map);
erlog = erlog(:,state.map);

sel_events = eventnames(state.map);
obs = state.obs;

dep = local_dod_summary(bp,bt,erlog,obs,'dep');
gen = local_dod_summary(bp,bt,erlog,obs,'gen');

local_plot_dod_outputs(bp,bt,erlog,obs,sel_events,groupnames{state.gidx},...
    state.gridmax,state.metric,dep,gen);

%Explain the numbers just drawn when one of the quantities that has to be a
%variance was not (RC-43). The DoD route builds no relsummary, so
%psyrat_admissibility_collect -- which walks relsummary.group(g).diffscore --
%cannot reach this record; before this it was written by psyrat_dodiffrel and
%read by nothing, leaving a GUI user with an uninterpretable coefficient and no
%console to see the warning in. Raised AFTER the figure, like
%psyrat_relfigures.m does, because it explains what is on screen rather than
%replacing it. Nothing was clipped to produce those numbers.
%
%dep only, not dep and gen: the two records are bit-identical, and the note sums
%counts across whatever it is given. See psyrat_dod_admissibility_note.
adm_note = psyrat_dod_admissibility_note(dep);
if ~isempty(adm_note)
    warndlg(adm_note,'Coefficient is not interpretable as a reliability');
end

end

function psyrat_dod_othermetrics(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

ind = find(strcmp('groupnames',varargin),1);
groupnames = varargin{ind+1};

ind = find(strcmp('eventnames',varargin),1);
eventnames = varargin{ind+1};

[gidx,ok] = local_current_group(inputs,numel(groupnames));
if ~ok
    errordlg({'Selected group index is invalid.';...
        'How to fix: refresh the DoD window and re-select a group.'},...
        'Invalid group selection','modal');
    return;
end

local_open_metrics_window(psyrat_prefs,psyrat_data,gidx,groupnames{gidx},eventnames);
end

function psyrat_dod_observedtable(varargin)
ind = find(strcmp('psyrat_data',varargin),1);
psyrat_data = varargin{ind+1};

ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

ind = find(strcmp('eventnames',varargin),1);
eventnames = varargin{ind+1};

[ok,state,errmsg,errtitle] = local_read_main_inputs(inputs,eventnames,false);
if ~ok
    errordlg(errmsg,errtitle,'modal');
    return;
end

try
    psyrat_dod_observedt('psyrat_data',psyrat_data,...
        'groupidx',state.gidx,...
        'map',state.map,...
        'gui',1);
catch ME
    errordlg({ME.message},'Unable to build observed DoD table','modal');
end
end

function psyrat_dod_exportvarcomp(varargin)
ind = find(strcmp('psyrat_data',varargin),1);
psyrat_data = varargin{ind+1};

ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

ind = find(strcmp('eventnames',varargin),1);
eventnames = varargin{ind+1};

[ok,state,errmsg,errtitle] = local_read_main_inputs(inputs,eventnames,false);
if ~ok
    errordlg(errmsg,errtitle,'modal');
    return;
end

try
    psyrat_dod_varcompt('psyrat_data',psyrat_data,...
        'groupidx',state.gidx,...
        'map',state.map,...
        'savefile',1,...
        'gui',0);
catch ME
    errordlg({ME.message},'Unable to export DoD variance components','modal');
end
end

function [groupnames,eventnames,trialstats] = local_dod_stats(rel)
eventnames = cellstr(string(rel.events(:)));

if ischar(rel.groups) && strcmpi(rel.groups,'none')
    groupnames = {'none'};
    hasgroup = false;
else
    groupnames = cellstr(string(rel.groups(:)));
    hasgroup = true;
end

ne = numel(eventnames);
ng = numel(groupnames);
trialstats.mean = zeros(ng,ne);
trialstats.max = zeros(ng,ne);

datatable = rel.data;
for g = 1:ng
    for e = 1:ne
        if hasgroup
            gmask = string(datatable.group) == string(groupnames{g});
            emask = string(datatable.event) == string(eventnames{e});
            subset = datatable(gmask & emask,:);
        else
            emask = string(datatable.event) == string(eventnames{e});
            subset = datatable(emask,:);
        end

        if isempty(subset)
            trialstats.mean(g,e) = 0;
            trialstats.max(g,e) = 0;
        else
            trltable = varfun(@length,subset,'GroupingVariables',{'id'});
            trialstats.mean(g,e) = mean(trltable.GroupCount);
            trialstats.max(g,e) = max(trltable.GroupCount);
        end
    end
end
end

function psyrat_prefs = local_dod_initprefs(psyrat_prefs,rel,groupnames,eventnames,trialstats)
if ~isfield(psyrat_prefs,'view')
    psyrat_prefs.view = struct;
end

ng = numel(groupnames);
ne = numel(eventnames);
relmap = local_dod_relmap(rel,eventnames);

if ~isfield(psyrat_prefs.view,'dodgroup') || ...
        isempty(psyrat_prefs.view.dodgroup) || ...
        psyrat_prefs.view.dodgroup < 1 || ...
        psyrat_prefs.view.dodgroup > ng
    psyrat_prefs.view.dodgroup = 1;
end

if ~isfield(psyrat_prefs.view,'dodmap') || ...
        isempty(psyrat_prefs.view.dodmap) || ...
        numel(psyrat_prefs.view.dodmap) ~= 4 || ...
        any(psyrat_prefs.view.dodmap < 1) || ...
        any(psyrat_prefs.view.dodmap > ne)
    if ~isempty(relmap)
        psyrat_prefs.view.dodmap = relmap;
    else
        psyrat_prefs.view.dodmap = 1:4;
    end
elseif ~isempty(relmap) && ...
        isequal(psyrat_prefs.view.dodmap,1:4) && ...
        ~isequal(relmap,1:4)
    %Use process-stage DoD mapping as the initial default when available.
    psyrat_prefs.view.dodmap = relmap;
end

gidx = psyrat_prefs.view.dodgroup;
defaults = max(round(trialstats.mean(gidx,psyrat_prefs.view.dodmap)),0);

if ~isfield(psyrat_prefs.view,'dodtrials') || ...
        isempty(psyrat_prefs.view.dodtrials) || ...
        numel(psyrat_prefs.view.dodtrials) ~= 4
    psyrat_prefs.view.dodtrials = defaults;
else
    t = round(psyrat_prefs.view.dodtrials(:)');
    t(t < 0) = 0;
    psyrat_prefs.view.dodtrials = t;
end

gridmax_default = max(trialstats.max(gidx,:)) + 25;
if ~isfield(psyrat_prefs.view,'dodgridmax') || ...
        isempty(psyrat_prefs.view.dodgridmax) || ...
        psyrat_prefs.view.dodgridmax < 0
    psyrat_prefs.view.dodgridmax = gridmax_default;
end

if ~isfield(psyrat_prefs.view,'dodheatmetric') || ...
        ~any(psyrat_prefs.view.dodheatmetric == [1 2])
    psyrat_prefs.view.dodheatmetric = 1;
end
end

function relmap = local_dod_relmap(rel,eventnames)
%Return DoD map indices aligned to rel.events order, if available.

relmap = [];
if ~isfield(rel,'dod_map') || numel(rel.dod_map) ~= 4
    return;
end

maplabels = cellstr(string(rel.dod_map(:)));
if numel(unique(maplabels)) ~= 4
    return;
end

relmap = zeros(1,4);
for i = 1:4
    idx = find(strcmp(maplabels{i},eventnames),1);
    if isempty(idx)
        relmap = [];
        return;
    end
    relmap(i) = idx;
end

if numel(unique(relmap)) ~= 4
    relmap = [];
end
end

function [bp,bt,erlog] = local_groupdraws(rel,gidx)
bp = rel.out.id_varcov{1,gidx};
bt = rel.out.trl_varcov{1,gidx};
erlog = rel.out.b_sigma{1,gidx};

if iscell(bp)
    bp = cell2mat(bp);
end
if iscell(bt)
    bt = cell2mat(bt);
end
if iscell(erlog)
    erlog = cell2mat(erlog);
end
end

function summary = local_dod_summary(bp,bt,erlog,obs,metric)
if any(obs <= 0)
    summary = struct(...
        'pt',0,'ll',0,'ul',0,...
        'icc_pt',NaN,'icc_ll',NaN,'icc_ul',NaN);
    return;
end

summary = psyrat_dodiffrel(...
    'bp',bp,...
    'bt',bt,...
    'er_var',erlog,...
    'obs',obs,...
    'est',metric,...
    'CI',.95);
end

function local_plot_dod_outputs(bp,bt,erlog,obs,sel_events,gname,gridmax,metric,dep,gen)
outfigs = findall(0,'Type','figure','Name','DoD Reliability Outputs');
if ~isempty(outfigs)
    close(outfigs);
end

fig = figure;
fig.Tag = 'psyrat_output';
set(fig,'NumberTitle','Off');
set(fig,'Name','DoD Reliability Outputs');
set(fig,'Position',[120 80 1180 760]);

metric_name = 'Dependability';
if metric == 2
    metric_name = 'Generalizability';
end

contrast_text = sprintf('Contrast: (%s - %s) - (%s - %s)',...
    sel_events{1},sel_events{2},sel_events{3},sel_events{4});
trial_text = sprintf('Entered trials: [%d %d %d %d]',obs(1),obs(2),obs(3),obs(4));

uicontrol(fig,'Style','text','FontSize',12,'HorizontalAlignment','left',...
    'String',['Group: ' gname],...
    'Position',[20 720 360 24]);
uicontrol(fig,'Style','text','FontSize',12,'HorizontalAlignment','left',...
    'String',contrast_text,...
    'Position',[20 694 760 24]);
%The ICC columns of the table below are evaluated with the UNSCALED contrast
%variances (psyrat_dodiffrel :151-152), i.e. at a single observation, so say so
%here rather than leave the reader to assume the entered trial counts apply to
%them. Same qualifier the formula catalog uses for the conditional ICCs.
uicontrol(fig,'Style','text','FontSize',11,'HorizontalAlignment','left',...
    'String',[trial_text ' | Heatmap metric: ' metric_name ...
    ' | ICC columns are single-observation (n''=1)'],...
    'Position',[20 670 900 24]);

if any(obs == 0)
    uicontrol(fig,'Style','text','FontSize',10,'ForegroundColor',[0.70 0.00 0.00],...
        'HorizontalAlignment','left',...
        'String','Note: one or more entered trial counts are 0; summary reliability is reported as 0.',...
        'Position',[20 646 760 22]);
end

tab = cell(2,7);
tab(1,:) = {'Dependability',dep.pt,dep.ll,dep.ul,dep.icc_pt,dep.icc_ll,dep.icc_ul};
tab(2,:) = {'Generalizability',gen.pt,gen.ll,gen.ul,gen.icc_pt,gen.icc_ll,gen.icc_ul};
tab = local_round_cells(tab,4);

uitable('Parent',fig,...
    'Data',tab,...
    'ColumnName',{'Coefficient','Point','Low CrI','High CrI','ICC Point','ICC Low','ICC High'},...
    'ColumnWidth',{140 90 90 90 90 90 90},...
    'RowName',[],...
    'Position',[20 500 700 130]);

grid = 0:gridmax;
bp_mean = squeeze(mean(bp,1));
bt_mean = squeeze(mean(bt,1));
er_mean = mean(exp(erlog).^2,1);

z12 = zeros(length(grid),length(grid));
z34 = zeros(length(grid),length(grid));

for xi = 1:length(grid)
    for yi = 1:length(grid)
        obs12 = [grid(xi) grid(yi) obs(3) obs(4)];
        obs34 = [obs(1) obs(2) grid(xi) grid(yi)];
        z12(yi,xi) = local_dod_point_metric(bp_mean,bt_mean,er_mean,obs12,metric);
        z34(yi,xi) = local_dod_point_metric(bp_mean,bt_mean,er_mean,obs34,metric);
    end
end

z12 = max(min(z12,1),0);
z34 = max(min(z34,1),0);

%Event names print literally: the default TeX interpreter would turn the
%underscore in a name such as inc_err into a subscript (owner ruling
%2026-09-04, extending the 2026-09-03 plotting-function ruling to this window).
ax1 = axes('Parent',fig,'Position',[0.08 0.10 0.38 0.32]);
imagesc(ax1,grid,grid,z12);
axis(ax1,'xy');
colormap(ax1,parula(256));
colorbar(ax1);
xlabel(ax1,['Trials ' sel_events{1}],'Interpreter','none');
ylabel(ax1,['Trials ' sel_events{2}],'Interpreter','none');
title(ax1,sprintf('%s Heatmap: %s vs %s',metric_name,sel_events{1},sel_events{2}),...
    'Interpreter','none');

ax2 = axes('Parent',fig,'Position',[0.55 0.10 0.38 0.32]);
imagesc(ax2,grid,grid,z34);
axis(ax2,'xy');
colormap(ax2,parula(256));
colorbar(ax2);
xlabel(ax2,['Trials ' sel_events{3}],'Interpreter','none');
ylabel(ax2,['Trials ' sel_events{4}],'Interpreter','none');
title(ax2,sprintf('%s Heatmap: %s vs %s',metric_name,sel_events{3},sel_events{4}),...
    'Interpreter','none');

uicontrol(fig,'Style','text','FontSize',10,'HorizontalAlignment','left',...
    'String',sprintf('%s/%s heatmap fixes %s=%d and %s=%d.',...
    sel_events{1},sel_events{2},sel_events{3},obs(3),sel_events{4},obs(4)),...
    'Position',[70 52 430 20]);
uicontrol(fig,'Style','text','FontSize',10,'HorizontalAlignment','left',...
    'String',sprintf('%s/%s heatmap fixes %s=%d and %s=%d.',...
    sel_events{3},sel_events{4},sel_events{1},obs(1),sel_events{2},obs(2)),...
    'Position',[630 52 430 20]);
end

function val = local_dod_point_metric(bp,bt,er,obs,metric)
if any(obs <= 0)
    val = 0;
    return;
end

cvec = [1; -1; -1; 1];
u = cvec' * bp * cvec;

bt_scaled = 0;
for i = 1:4
    bt_scaled = bt_scaled + (cvec(i)^2) * (bt(i,i) / obs(i));
end

for i = 1:3
    for j = i+1:4
        nharm = 2 / ((1/obs(i)) + (1/obs(j)));
        bt_scaled = bt_scaled + (2 * cvec(i) * cvec(j) * bt(i,j) / nharm);
    end
end

er_scaled = 0;
for i = 1:4
    er_scaled = er_scaled + (cvec(i)^2) * (er(i) / obs(i));
end

if metric == 1
    denom = u + bt_scaled + er_scaled;
else
    denom = u + er_scaled;
end

if denom <= 0 || ~isfinite(denom)
    val = NaN;
else
    val = u / denom;
end
end

function arr = local_round_cells(arr,n)
scale = 10^n;
for r = 1:size(arr,1)
    for c = 1:size(arr,2)
        if isnumeric(arr{r,c}) && isscalar(arr{r,c}) && ~isnan(arr{r,c})
            arr{r,c} = round(arr{r,c} * scale) / scale;
        end
    end
end
end

function val = local_nonnegint(str)
%isreal rejects complex literals like '50-4i': str2double parses them, isnan is
%false for them, val < 0 compares only the real part, and round() acts on both
%parts so abs(round(val)-val) is 0. A complex trial count therefore passed every
%check here and at the any(isnan([t1 t2 t3 t4 gridmax])) gate in
%local_read_main_inputs below, then reached psyrat_dodiffrel and threw
%"First argument must be an array of real values" from quantile; a complex
%gridmax threw "Colon operands must be real scalars" at 0:gridmax instead.
val = str2double(str);
if ~isreal(val) || isnan(val) || val < 0 || abs(round(val) - val) > 1e-9
    val = NaN;
else
    val = round(val);
end
end

function [ok,state,errmsg,errtitle] = local_read_main_inputs(inputs,eventnames,need_trial_inputs)
ok = false;
state = struct;
errmsg = {};
errtitle = 'Invalid DoD inputs';

[gidx,good_group] = local_current_group(inputs,inf);
if ~good_group
    errmsg = {'Selected group index is invalid.';...
        'How to fix: refresh the DoD window and re-select a group.'};
    errtitle = 'Invalid group selection';
    return;
end

map = [inputs.erp1.Value inputs.erp2.Value inputs.erp3.Value inputs.erp4.Value];
if any(map < 1) || any(map > numel(eventnames))
    errmsg = {'ERP mapping contains an invalid event index.';...
        'How to fix: refresh the DoD window and reselect ERP_1 to ERP_4.'};
    errtitle = 'Invalid DoD mapping';
    return;
end

if numel(unique(map)) ~= 4
    errmsg = {'ERP_1 to ERP_4 must map to four unique events.';...
        'How to fix: select each event only once across ERP_1, ERP_2, ERP_3, and ERP_4.'};
    errtitle = 'DoD mapping requires four unique events';
    return;
end

state.gidx = gidx;
state.map = map;
state.metric = inputs.metric.Value;

if need_trial_inputs
    t1 = local_nonnegint(inputs.trials1.String);
    t2 = local_nonnegint(inputs.trials2.String);
    t3 = local_nonnegint(inputs.trials3.String);
    t4 = local_nonnegint(inputs.trials4.String);
    gridmax = local_nonnegint(inputs.gridmax.String);

    if any(isnan([t1 t2 t3 t4 gridmax]))
        errmsg = {'Trial inputs and heatmap max must be whole numbers >= 0.';...
            'How to fix: enter non-negative integers only.'};
        errtitle = 'Invalid trial input';
        return;
    end

    if gridmax > 300
        errmsg = {'Heatmap grid max is too large for interactive plotting (>300).';...
            'How to fix: choose a grid max of 300 or lower.'};
        errtitle = 'Heatmap grid too large';
        return;
    end

    state.obs = [t1 t2 t3 t4];
    state.gridmax = gridmax;
end

ok = true;
end

function [gidx,ok] = local_current_group(inputs,max_groups)
ok = true;
if isempty(inputs.group)
    gidx = 1;
    inferred_max = 1;
else
    gidx = inputs.group.Value;
    inferred_max = numel(inputs.group.String);
end
if ~isscalar(gidx) || ~isfinite(gidx) || gidx < 1 || abs(gidx-round(gidx)) > 0
    ok = false;
end
gidx = round(gidx);
if nargin < 2 || isempty(max_groups) || ~isfinite(max_groups)
    max_groups = inferred_max;
end
if gidx > max_groups
    ok = false;
end
end

function local_open_metrics_window(psyrat_prefs,psyrat_data,gidx,gname,eventnames)
aux_gui = findobj('Tag','psyrat_dod_metrics');
if ~isempty(aux_gui)
    close(aux_gui);
end

fsize = psyrat_prefs.guis.fsize;
figwidth = 520;
figheight = 300;
rowspace = 40;
row = figheight - 55;
lcol = 30;
rcol = 220;

aux_gui = psyrat_newfigure('unit','pix',...
    'position',[430 410 figwidth figheight],...
    'menub','no',...
    'name','DoD: Other Single-Session Metrics',...
    'numbertitle','off',...
    'resize','off');

uicontrol(aux_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Other Single-Session Metrics from DoD Model',...
    'Position',[0 row figwidth 24]);
row = row - rowspace*0.9;

uicontrol(aux_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String',['Selected group: ' gname],...
    'Position',[lcol row figwidth-2*lcol 24]);
row = row - rowspace;

uicontrol(aux_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Mode:',...
    'Position',[lcol row 170 24]);

aux.mode = uicontrol(aux_gui,'Style','pop',...
    'fontsize',fsize,...
    'String',{'Single-event reliability';'2-event difference reliability'},...
    'Value',1,...
    'Position',[rcol row 250 26],...
    'Callback',{@local_dod_modechange,'aux_gui',aux_gui});
row = row - rowspace;

aux.event1_label = uicontrol(aux_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Event 1:',...
    'Position',[lcol row 170 24]);

aux.event1 = uicontrol(aux_gui,'Style','pop',...
    'fontsize',fsize,...
    'String',eventnames(:),...
    'Value',1,...
    'Position',[rcol row 250 26]);
row = row - rowspace;

aux.event2_label = uicontrol(aux_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Event 2:',...
    'Position',[lcol row 170 24],...
    'Visible','off');

aux.event2 = uicontrol(aux_gui,'Style','pop',...
    'fontsize',fsize,...
    'String',eventnames(:),...
    'Value',min(2,numel(eventnames)),...
    'Position',[rcol row 250 26],...
    'Visible','off');

uicontrol(aux_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Open Standard Output GUI',...
    'Position',[figwidth/8 18 3*figwidth/4 46],...
    'Callback',{@local_dod_open_standard,...
    'aux_gui',aux_gui,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'groupidx',gidx});

setappdata(aux_gui,'dod_aux_inputs',aux);
aux_gui.Tag = 'psyrat_dod_metrics';

end

function local_dod_modechange(varargin)
ind = find(strcmp('aux_gui',varargin),1);
aux_gui = varargin{ind+1};
aux = getappdata(aux_gui,'dod_aux_inputs');

if aux.mode.Value == 1
    aux.event2_label.Visible = 'off';
    aux.event2.Visible = 'off';
else
    aux.event2_label.Visible = 'on';
    aux.event2.Visible = 'on';
end
end

function local_dod_open_standard(varargin)
ind = find(strcmp('aux_gui',varargin),1);
aux_gui = varargin{ind+1};
aux = getappdata(aux_gui,'dod_aux_inputs');

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

ind = find(strcmp('groupidx',varargin),1);
gidx = varargin{ind+1};

modeval = aux.mode.Value;

try
    if modeval == 1
        eventidx = aux.event1.Value;
        psyrat_data_out = psyrat_dod_build_virtual_rel(...
            'psyrat_data',psyrat_data,...
            'mode','single',...
            'groupidx',gidx,...
            'eventidx',eventidx);
    else
        eventpair = [aux.event1.Value aux.event2.Value];
        if numel(unique(eventpair)) ~= 2
            errordlg({'2-event difference reliability requires two unique events.';...
                'How to fix: pick different events for Event 1 and Event 2.'},...
                'Difference score requires unique events','modal');
            return;
        end

        psyrat_data_out = psyrat_dod_build_virtual_rel(...
            'psyrat_data',psyrat_data,...
            'mode','diff',...
            'groupidx',gidx,...
            'eventpair',eventpair);
    end
catch ME
    errordlg({ME.message},'Unable to build standard-output view','modal');
    return;
end

if ishandle(aux_gui)
    close(aux_gui);
end

psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data_out);
end

function psyrat_dod_closefigs(varargin)
outfigs = findall(0,'Type','figure','Tag','psyrat_output');
if ~isempty(outfigs)
    close(outfigs);
end
end
