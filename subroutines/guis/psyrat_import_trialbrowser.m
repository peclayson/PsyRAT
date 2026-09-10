function psyrat_import_trialbrowser(varargin)
%Browse imported single-trial ERP waveforms trial-by-trial.
%
% psyrat_import_trialbrowser('project',project,'spec_id',specId)

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

opts = struct('project',[],'spec_id','');
if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end
for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'project'
            opts.project = val;
        case 'spec_id'
            opts.spec_id = char(string(val));
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(opts.project)
    error('psyrat_browser:missingProject','Input ''project'' is required.');
end
project = opts.project;

if ~isfield(project,'canonical') || ~project.canonical.has_single_trial
    error('psyrat_browser:noSingleTrial', ...
        'Project does not contain single-trial data to browse.');
end

signals = project.canonical.single_trial.signals;
time_ms = project.canonical.single_trial.time_ms(:);
labels = project.canonical.single_trial.channel_labels;
trial_table = project.canonical.trial_table;

nTrials = size(signals,1);
nChans = size(signals,2);
if nTrials < 1
    error('psyrat_browser:emptyData','No trials were found in canonical single-trial data.');
end

spec = [];
if ~isempty(opts.spec_id) && ~isempty(project.scoring_specs)
    idx = find(strcmp({project.scoring_specs.spec_id},opts.spec_id),1,'first');
    if ~isempty(idx)
        spec = project.scoring_specs(idx);
    end
end

if isempty(spec)
    displayCh = 1:min(8,nChans);
    specName = '(none)';
else
    displayCh = localResolveChannels(spec.channel_selector,labels);
    if numel(displayCh) > 16
        displayCh = displayCh(1:16);
    end
    specName = spec.name;
end

figW = 980;
figH = 620;
hFig = psyrat_newfigure('Visible','off','unit','pix', ...
    'position',[300 120 figW figH],'menub','none','numbertitle','off', ...
    'resize','off','Name','PsyRAT Trial Browser');

%TAG REQUIRED, not decoration (B48). psyrat_start's teardown closes PsyRAT
%windows by 'psyrat_' Tag prefix instead of the bare `close all` it used to
%run, so an untagged window is no longer swept up. This one is the concrete
%case: psyrat_startimportscore's localBackToStart closes only
%psyrat_import_gui and then calls psyrat_start, and no PsyRAT code closes
%the browser explicitly -- this file has no close, delete or CloseRequestFcn
%anywhere -- so untagged it would be stranded for the rest of the session.
hFig.Tag = 'psyrat_gui_trialbrowser';

movegui(hFig,'center');

hAxes = axes('Parent',hFig,'Units','pixels','Position',[60 200 900 380]);

hTrialTxt = uicontrol(hFig,'Style','text','HorizontalAlignment','left', ...
    'Position',[60 165 900 24],'String','');
hMetaTxt = uicontrol(hFig,'Style','text','HorizontalAlignment','left', ...
    'Position',[60 140 900 24],'String','');
hScoreTxt = uicontrol(hFig,'Style','text','HorizontalAlignment','left', ...
    'Position',[60 115 900 24],'String','');

uicontrol(hFig,'Style','text','String','Trial #','HorizontalAlignment','right', ...
    'Position',[60 70 70 24]);

hTrialEdit = uicontrol(hFig,'Style','edit','String','1', ...
    'Position',[135 70 80 26],'Callback',@onJumpToTrial, ...
    'TooltipString','Type a trial number and press Enter to jump directly.');

uicontrol(hFig,'Style','push','String','Prev', ...
    'Position',[230 70 90 28], ...
    'TooltipString','Go to previous trial.', ...
    'Callback',@(~,~) onStepTrial(-1));

uicontrol(hFig,'Style','push','String','Next', ...
    'Position',[330 70 90 28], ...
    'TooltipString','Go to next trial.', ...
    'Callback',@(~,~) onStepTrial(1));

hSlider = uicontrol(hFig,'Style','slider','Min',1,'Max',nTrials,'Value',1, ...
    'Position',[440 70 520 28], ...
    'TooltipString','Drag to cycle quickly through trials.', ...
    'Callback',@onSlider);

if nTrials > 1
    set(hSlider,'SliderStep',[1/(nTrials-1) min(50/(nTrials-1),1)]);
else
    set(hSlider,'Enable','off');
end

state = struct();
state.project = project;
state.spec = spec;
state.specName = specName;
state.displayCh = displayCh;
state.labels = labels;
state.signals = signals;
state.time_ms = time_ms;
state.trial_table = trial_table;
state.nTrials = nTrials;
state.currTrial = 1;
state.axes = hAxes;
state.hTrialTxt = hTrialTxt;
state.hMetaTxt = hMetaTxt;
state.hScoreTxt = hScoreTxt;
state.hTrialEdit = hTrialEdit;
state.hSlider = hSlider;

guidata(hFig,state);
set(hFig,'KeyPressFcn',@onKey);
set(hFig,'Visible','on');
localRefresh(hFig);

    function onStepTrial(delta)
        st = guidata(hFig);
        st.currTrial = max(1,min(st.nTrials,st.currTrial + delta));
        guidata(hFig,st);
        localRefresh(hFig);
    end

    function onJumpToTrial(~,~)
        st = guidata(hFig);
        v = str2double(get(st.hTrialEdit,'String'));
        %~isreal is needed alongside ~isfinite: str2double parses complex
        %literals like '3-4i' and isfinite is true for them, while round keeps
        %both parts and min/max compare by magnitude, so the clamp below does
        %not restore a real value. A complex currTrial then reached
        %localRefresh, where set(hSlider,'Value',...) throws first ("Complex
        %inputs are not supported"); that aborts the callback before
        %st.signals(t,...) a couple of lines below, which would otherwise
        %throw in its own right ("Array indices must be positive integers").
        %Nothing in this file is wrapped in try/catch.
        if ~isreal(v) || ~isfinite(v)
            return;
        end
        st.currTrial = max(1,min(st.nTrials,round(v)));
        guidata(hFig,st);
        localRefresh(hFig);
    end

    function onSlider(~,~)
        st = guidata(hFig);
        st.currTrial = max(1,min(st.nTrials,round(get(st.hSlider,'Value'))));
        guidata(hFig,st);
        localRefresh(hFig);
    end

    function onKey(~,evt)
        switch evt.Key
            case {'rightarrow','downarrow'}
                onStepTrial(1);
            case {'leftarrow','uparrow'}
                onStepTrial(-1);
        end
    end

end

function localRefresh(hFig)
st = guidata(hFig);
t = st.currTrial;

set(st.hTrialEdit,'String',num2str(t));
set(st.hSlider,'Value',t);

wave = squeeze(st.signals(t,st.displayCh,:));
if numel(st.displayCh) == 1
    wave = wave(:)';
else
    wave = squeeze(wave);
end

axes(st.axes);
cla(st.axes);
hold(st.axes,'on');

if numel(st.displayCh) == 1
    plot(st.axes,st.time_ms,wave,'Color',[0.1 0.4 0.85],'LineWidth',1.8);
else
    cmap = lines(size(wave,1));
    for i = 1:size(wave,1)
        plot(st.axes,st.time_ms,wave(i,:),'Color',cmap(i,:),'LineWidth',1.0);
    end
    meanWave = mean(wave,1,'omitnan');
    plot(st.axes,st.time_ms,meanWave,'k','LineWidth',2.0);
end

if ~isempty(st.spec) && isfield(st.spec,'window_start_ms') && isfield(st.spec,'window_end_ms')
    y = ylim(st.axes);
    xs = [st.spec.window_start_ms st.spec.window_end_ms st.spec.window_end_ms st.spec.window_start_ms];
    ys = [y(1) y(1) y(2) y(2)];
    patch('Parent',st.axes,'XData',xs,'YData',ys,'FaceColor',[0.95 0.88 0.55], ...
        'FaceAlpha',0.25,'EdgeColor','none');
    uistack(findobj(st.axes,'Type','line'),'top');
end

xlabel(st.axes,'Time (ms)');
ylabel(st.axes,'Amplitude (uV)');
grid(st.axes,'on');

chNames = st.labels(st.displayCh);
if iscell(chNames)
    cstr = strjoin(chNames,', ');
else
    cstr = strjoin(cellstr(chNames),', ');
end

title(st.axes,sprintf('Trial %d/%d | Channels: %s | Spec: %s', ...
    t,st.nTrials,cstr,st.specName),'Interpreter','none');

meta = st.trial_table(t,:);
set(st.hTrialTxt,'String',sprintf('Trial %d of %d',t,st.nTrials));
set(st.hMetaTxt,'String',sprintf('subject=%s | event(cell)=%s | occasion=%s | group=%s | uid=%s', ...
    localCell(meta.subject_id),localCell(meta.event_id),localCell(meta.occasion_id), ...
    localCell(meta.group_id),localCell(meta.trial_uid)));

scoreMsg = localScoreMessage(st.project,st.spec,t,meta.trial_uid{1});
set(st.hScoreTxt,'String',scoreMsg);

drawnow;
end

function s = localScoreMessage(project,spec,trialIdx,trialUid)
if isempty(spec) || isempty(project.scoring_results)
    s = sprintf('Score: not available (run scoring to see trial-level results). [trial index %d]',trialIdx);
    return;
end

tbl = project.scoring_results;
if any(strcmp(tbl.Properties.VariableNames,'is_stale'))
    tbl = tbl(~tbl.is_stale,:);
end

ind = strcmp(tbl.spec_id,spec.spec_id) & strcmp(tbl.trial_uid,trialUid);
if ~any(ind)
    s = sprintf('Score: not available for spec "%s" on this trial.',spec.name);
    return;
end

r = tbl(find(ind,1,'last'),:);
if r.excluded_from_handoff
    s = sprintf('Score: EXCLUDED | status=%s | reason=%s',r.status{1},r.status_reason{1});
else
    if isfinite(r.latency_ms)
        s = sprintf('Score: value=%.6g | latency=%.3f ms | status=%s',r.value,r.latency_ms,r.status{1});
    else
        s = sprintf('Score: value=%.6g | status=%s',r.value,r.status{1});
    end
end
end

function out = localCell(in)
if iscell(in)
    if isempty(in)
        out = '';
    else
        out = char(string(in{1}));
    end
else
    out = char(string(in));
end
end

function ch = localResolveChannels(selector,labels)
if isnumeric(selector)
    ch = selector(:)';
    ch = unique(ch,'stable');
    ch = ch(ch >= 1 & ch <= numel(labels));
    if isempty(ch)
        ch = 1;
    end
    return;
end

if ischar(selector) || (isstring(selector) && isscalar(selector))
    selector = {char(string(selector))};
elseif isstring(selector)
    selector = cellstr(selector(:));
elseif ~iscell(selector)
    ch = 1;
    return;
end

labelsLower = lower(string(labels(:)));
ch = zeros(1,numel(selector));
k = 0;
for i = 1:numel(selector)
    idx = find(labelsLower == lower(string(selector{i})),1,'first');
    if ~isempty(idx)
        k = k + 1;
        ch(k) = idx;
    end
end
if k == 0
    ch = 1;
else
    ch = unique(ch(1:k),'stable');
end
end
