function psyrat_data_out = psyrat_dod_build_virtual_rel(varargin)
%Build a virtual single-session REL structure from DoD model outputs
%
%psyrat_data_out = psyrat_dod_build_virtual_rel('psyrat_data',psyrat_data,...
%  'mode','single','groupidx',1,'eventidx',2)
%
%psyrat_data_out = psyrat_dod_build_virtual_rel('psyrat_data',psyrat_data,...
%  'mode','diff','groupidx',1,'eventpair',[1 3])
%
%Required Inputs:
% psyrat_data - PsyRAT data structure containing rel.analysis='ic_dodiff'
% mode - 'single' for single-event outputs, 'diff' for 2-event difference outputs
% groupidx - index of group in rel.groups (use 1 when rel.groups='none')
%
%For mode='single':
% eventidx - index of selected event in rel.events
%
%For mode='diff':
% eventpair - [event1 event2] unique event indices in rel.events
%
%Output:
% psyrat_data_out - copy of psyrat_data with rel adapted for psyrat_startview_sing
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

if mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        ['Inputs are incomplete. Provide name/value pairs, for example ',...
        '''psyrat_dod_build_virtual_rel(''psyrat_data'',psyrat_data,', ...
        '''mode'',''single'',''groupidx'',1,''eventidx'',1)''.']);
end

psyrat_data = local_req(varargin,'psyrat_data',...
    'psyrat_data was not provided.');
mode = char(lower(string(local_req(varargin,'mode',...
    'mode was not provided. Use ''single'' or ''diff''.'))));
groupidx = local_req(varargin,'groupidx',...
    'groupidx was not provided.');

if ~isstruct(psyrat_data) || ~isfield(psyrat_data,'rel') || ...
        ~isfield(psyrat_data.rel,'analysis') || ...
        ~strcmpi(psyrat_data.rel.analysis,'ic_dodiff')
    error('varargin:analysis',... %Error code and associated error
        ['Virtual REL building requires DoD model outputs ',...
        '(rel.analysis = ''ic_dodiff'').']);
end

if ~isscalar(groupidx) || ~isnumeric(groupidx) || isnan(groupidx) || ...
        groupidx < 1 || abs(groupidx-round(groupidx)) > 0
    error('varargin:groupidx',... %Error code and associated error
        'groupidx must be a positive integer scalar.');
end
groupidx = round(groupidx);

rel_in = psyrat_data.rel;
[groupnames,hasgroup] = local_groupnames(rel_in);
eventnames = cellstr(string(rel_in.events(:)));

if groupidx > numel(groupnames)
    error('varargin:groupidx',... %Error code and associated error
        'groupidx exceeds the number of groups in the DoD outputs.');
end

if numel(eventnames) < 2
    error('varargin:eventcount',... %Error code and associated error
        'At least 2 events are required in DoD outputs.');
end

[id_varcov,trl_varcov,b_sigma,b_cell] = local_groupdraws(rel_in,groupidx);
ndraws = size(b_sigma,1);

switch mode
    case {'single','sing'}
        eventidx = local_req(varargin,'eventidx',...
            'eventidx is required for mode=''single''.');
        if ~isscalar(eventidx) || ~isnumeric(eventidx) || isnan(eventidx) || ...
                eventidx < 1 || eventidx > numel(eventnames) || ...
                abs(eventidx-round(eventidx)) > 0
            error('varargin:eventidx',... %Error code and associated error
                'eventidx must be an integer from 1 to the number of events.');
        end
        eventidx = round(eventidx);

        event_label = eventnames{eventidx};
        data_sub = local_filter_data(rel_in,hasgroup,groupnames{groupidx},event_label);
        if isempty(data_sub)
            error('varargin:data',... %Error code and associated error
                ['No rows were found for selected group/event in DoD data. ',...
                'How to fix: choose an event that exists in the selected group.']);
        end

        rel_out = rel_in;
        rel_out.analysis = 'ic';
        rel_out.groups = 'none';
        rel_out.events = {event_label};
        rel_out.data = data_sub;
        rel_out.dod_virtual_source = 'ic_dodiff';
        rel_out.dod_virtual_group = groupnames{groupidx};
        rel_out.dod_virtual_events = {event_label};

        sig_u = sqrt(max(id_varcov(:,eventidx,eventidx),0));
        sig_e = exp(b_sigma(:,eventidx));
        %Trial main effect (sigma_i) for the SELECTED cell only: the one-facet
        %model this virtual REL feeds is a single event, so only the diagonal
        %element of the DoD trial covariance belongs to it (mode='diff' below
        %needs the full 2x2 because it forms a difference across two cells).
        %Without this, psyrat_relsummary backfills sig_trl with zeros and every
        %downstream psyrat_rel_sing call returns the relative coefficient under
        %an absolute label -- dependability would equal generalizability exactly
        %across the whole view (Rocha 2026 Table 2; RC-15).
        sig_trl = sqrt(max(trl_varcov(:,eventidx,eventidx),0));

        rel_out.out = struct;
        %The universe-score mean. b_cell is the per-cell mean intercept of the
        %Gaussian DoD model; psyrat_cutscore_sing forms (mu - cut)^2 from it, so
        %a placeholder zero would evaluate every criterion coefficient as though
        %the mean were exactly 0 (RC-18). Older result structs (saved before the
        %extraction was added) and the gamma family (whose DoD model parameterizes
        %a per-cell LOG-mean, not an observed-scale mean) carry no b_cell; those
        %keep the legacy zeros and are flagged so the criterion button can be
        %suppressed rather than shown with a fabricated mean.
        if isempty(b_cell)
            rel_out.out.mu = zeros(ndraws,1);
            rel_out.out.mu_source = 'unavailable';
        else
            rel_out.out.mu = b_cell(:,eventidx);
            rel_out.out.mu_source = 'estimated';
        end
        rel_out.out.sig_u = sig_u(:);
        rel_out.out.sig_e = sig_e(:);
        rel_out.out.sig_trl = sig_trl(:);
        rel_out.out.sig_trl_source = 'estimated';
        rel_out.out.labels = event_label;
        rel_out.out.conv.data = local_group_conv(rel_in,groupidx);

    case {'diff','difference'}
        eventpair = local_req(varargin,'eventpair',...
            'eventpair is required for mode=''diff''.');
        if ~isnumeric(eventpair) || numel(eventpair) ~= 2
            error('varargin:eventpair',... %Error code and associated error
                'eventpair must be a numeric vector [event1 event2].');
        end
        eventpair = round(eventpair(:)');
        if any(eventpair < 1) || any(eventpair > numel(eventnames))
            error('varargin:eventpair',... %Error code and associated error
                'eventpair indices must be within the available event range.');
        end
        if numel(unique(eventpair)) ~= 2
            error('varargin:eventpair',... %Error code and associated error
                'eventpair must contain two unique event indices.');
        end

        sel_events = eventnames(eventpair);
        data_sub = local_filter_data(rel_in,hasgroup,groupnames{groupidx},sel_events);
        if isempty(data_sub)
            error('varargin:data',... %Error code and associated error
                ['No rows were found for selected group/event pair in DoD data. ',...
                'How to fix: choose events that exist in the selected group.']);
        end

        id_sub = id_varcov(:,eventpair,eventpair);
        trl_sub = trl_varcov(:,eventpair,eventpair);
        b_sigma_sub = b_sigma(:,eventpair);
        sd_id_sub = zeros(ndraws,2);
        sd_trl_sub = zeros(ndraws,2);
        for k = 1:2
            sd_id_sub(:,k) = sqrt(max(id_sub(:,k,k),0));
            sd_trl_sub(:,k) = sqrt(max(trl_sub(:,k,k),0));
        end

        b_draws = zeros(ndraws,2);
        rescor = zeros(ndraws,1);
        err_varcov = zeros(ndraws,2,2);
        err_varcov(:,1,1) = exp(b_sigma_sub(:,1)).^2;
        err_varcov(:,2,2) = exp(b_sigma_sub(:,2)).^2;

        rel_out = rel_in;
        rel_out.analysis = 'ic_diff';
        rel_out.groups = 'none';
        rel_out.events = sel_events(:);
        rel_out.data = data_sub;
        rel_out.diff_names = sel_events(:);
        rel_out.diffwpcov = 1;
        rel_out.diffrescor = 1;
        rel_out.dod_virtual_source = 'ic_dodiff';
        rel_out.dod_virtual_group = groupnames{groupidx};
        rel_out.dod_virtual_events = sel_events(:);

        rel_out.out = struct;
        rel_out.out.b = {num2cell(b_draws)};
        rel_out.out.b_sigma = {num2cell(b_sigma_sub)};
        rel_out.out.sd_id = {num2cell(sd_id_sub)};
        rel_out.out.sd_trl = {num2cell(sd_trl_sub)};
        rel_out.out.id_varcov = {num2cell(id_sub)};
        rel_out.out.trl_varcov = {num2cell(trl_sub)};
        rel_out.out.rescor = {num2cell(rescor)};
        rel_out.out.err_varcov = {num2cell(err_varcov)};
        rel_out.out.elabels = sel_events(:);
        rel_out.out.glabels = {'none'};
        conv_entry = local_group_conv(rel_in,groupidx);
        rel_out.out.conv.data = {conv_entry};

    otherwise
        error('varargin:mode',... %Error code and associated error
            'mode must be ''single'' or ''diff''.');
end

psyrat_data_out = psyrat_data;
psyrat_data_out.rel = rel_out;
if isfield(psyrat_data_out,'relsummary')
    psyrat_data_out = rmfield(psyrat_data_out,'relsummary');
end

end

function val = local_req(args,name,msg)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    error(['varargin:' lower(name)],... %Error code and associated error
        msg);
end
val = args{ind+1};
end

function [groupnames,hasgroup] = local_groupnames(rel)
if ischar(rel.groups) && strcmpi(rel.groups,'none')
    groupnames = {'none'};
    hasgroup = false;
else
    groupnames = cellstr(string(rel.groups(:)));
    hasgroup = true;
end
end

function [id_varcov,trl_varcov,b_sigma,b_cell] = local_groupdraws(rel,groupidx)
id_varcov = rel.out.id_varcov{1,groupidx};
trl_varcov = rel.out.trl_varcov{1,groupidx};
b_sigma = rel.out.b_sigma{1,groupidx};

%b_cell is optional: it is absent from result structs saved before the DoD
%extraction stored it, and from the gamma family (see the mu comment above).
%Return [] in that case so the caller can flag the mean as unavailable.
b_cell = [];
if isfield(rel.out,'b_cell') && size(rel.out.b_cell,2) >= groupidx && ...
        ~isempty(rel.out.b_cell{1,groupidx})
    b_cell = rel.out.b_cell{1,groupidx};
    if iscell(b_cell)
        b_cell = cell2mat(b_cell);
    end
end

if iscell(id_varcov)
    id_varcov = cell2mat(id_varcov);
end
if iscell(trl_varcov)
    trl_varcov = cell2mat(trl_varcov);
end
if iscell(b_sigma)
    b_sigma = cell2mat(b_sigma);
end

if ndims(id_varcov) ~= 3 || ndims(trl_varcov) ~= 3 || size(b_sigma,2) < 2
    error('varargin:drawshape',... %Error code and associated error
        'DoD draws are malformed. Expected id/trial covariance matrices and residual log-SDs.');
end
end

function conv_data = local_group_conv(rel,groupidx)
conv_data = [];
if isfield(rel,'out') && isfield(rel.out,'conv') && ...
        isfield(rel.out.conv,'data') && numel(rel.out.conv.data) >= groupidx
    conv_data = rel.out.conv.data{groupidx};
end
end

function data_sub = local_filter_data(rel,hasgroup,glabel,event_labels)
datatable = rel.data;
if ischar(event_labels) || isstring(event_labels)
    event_labels = cellstr(string(event_labels));
else
    event_labels = cellstr(string(event_labels(:)));
end

emask = ismember(string(datatable.event),string(event_labels));
if hasgroup
    gmask = string(datatable.group) == string(glabel);
else
    gmask = true(height(datatable),1);
end

data_sub = datatable(gmask & emask,:);
end
