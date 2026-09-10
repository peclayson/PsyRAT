function out = psyrat_adapter_eeglab_read(filepath,varargin)
%Read EEGLAB dataset into canonical contribution struct.

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

% varargin is accepted for interface parity with the other import adapters.
if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

out = localEmptyOut('eeglab');

if exist('pop_loadset','file') ~= 2
    error('psyrat_import:eeglabUnavailable', ...
        ['EEGLAB loader function pop_loadset was not found in the MATLAB path. ',...
        'How to fix: add EEGLAB to the path or use another supported format.']);
end

EEG = pop_loadset(filepath);
if ~isstruct(EEG) || ~isfield(EEG,'data')
    error('psyrat_import:eeglabLoadFailed', ...
        'EEGLAB file did not load to a valid EEG struct: %s', filepath);
end

data = double(EEG.data);

if ndims(data) == 3 && size(data,3) > 1
    signals = permute(data,[3 1 2]);
    nTrials = size(signals,1);
    nChan = size(signals,2);
    nSamp = size(signals,3);

    if isfield(EEG,'times') && ~isempty(EEG.times)
        time_ms = double(EEG.times(:).');
    else
        time_ms = localTimeFromSrate(nSamp,EEG);
    end

    channel_labels = localChannelLabels(EEG,nChan);
    fs_hz = localFs(EEG,time_ms);

    [subject_id,event_id,occasion_id,group_id,session_id,nUnresolvedEvents] = ...
        localTrialMeta(EEG,nTrials,filepath);

    if nUnresolvedEvents > 0
        out.warnings{end+1} = sprintf(['%d of %d EEGLAB epochs contained multiple ',...
            'events but no usable eventlatency to identify the time-locking event; ',...
            'the first listed event was used for those trials. Verify the condition ',...
            '(event) labels.'],nUnresolvedEvents,nTrials);
    end

    trial_table = psyrat_build_trial_table(nTrials, ...
        'subject_id',subject_id, ...
        'event_id',event_id, ...
        'occasion_id',occasion_id, ...
        'group_id',group_id, ...
        'session_id',session_id, ...
        'source_ref',filepath);

    out.single_trial.available = true;
    out.single_trial.signals = signals;
    out.single_trial.time_ms = time_ms;
    out.single_trial.fs_hz = fs_hz;
    out.single_trial.amplitude_unit = 'uV';
    out.single_trial.channel_labels = channel_labels;
    out.single_trial.trial_table = trial_table;
    out.single_trial.data_hash = psyrat_hash_struct(struct( ...
        'signals_size',size(signals), ...
        'time_ms',time_ms, ...
        'channel_labels',{channel_labels}));

    avg = squeeze(mean(signals,1,'omitnan'));
    out.average.available = true;
    out.average.signals = avg;
    out.average.time_ms = time_ms;
    out.average.fs_hz = fs_hz;
    out.average.channel_labels = channel_labels;
    out.average.source = 'derived_from_single_trial';
    out.average.data_hash = psyrat_hash_struct(struct('signals',avg,'time_ms',time_ms));
else
    nChan = size(data,1);
    nSamp = size(data,2);

    if isfield(EEG,'times') && ~isempty(EEG.times)
        time_ms = double(EEG.times(:).');
    else
        time_ms = localTimeFromSrate(nSamp,EEG);
    end

    out.average.available = true;
    out.average.signals = data(:,1:nSamp);
    out.average.time_ms = time_ms;
    out.average.fs_hz = localFs(EEG,time_ms);
    out.average.channel_labels = localChannelLabels(EEG,nChan);
    out.average.source = 'imported_average';
    out.average.data_hash = psyrat_hash_struct(struct('signals',out.average.signals,'time_ms',time_ms));

    out.warnings{end+1} = ['EEGLAB source did not provide multi-trial waveforms. ',...
        'This import contributes averaged/reference data only and cannot be used for scoring/handoff.'];
end

end

function out = localEmptyOut(adapter)
out = struct();
out.adapter = adapter;
out.single_trial = struct('available',false,'signals',[],'time_ms',[],'fs_hz',[], ...
    'amplitude_unit','','channel_labels',{{}},'trial_table',table(),'data_hash','');
out.average = struct('available',false,'signals',[],'time_ms',[],'fs_hz',[], ...
    'channel_labels',{{}},'source','','data_hash','');
out.warnings = {};
out.notes = '';
end

function time_ms = localTimeFromSrate(nSamp,EEG)
fs = localFs(EEG,[]);
if isfield(EEG,'xmin') && ~isempty(EEG.xmin)
    start_ms = double(EEG.xmin) * 1000;
else
    start_ms = 0;
end
time_ms = start_ms + ((0:(nSamp-1))/fs)*1000;
end

function fs = localFs(EEG,time_ms)
fs = [];
if isfield(EEG,'srate') && ~isempty(EEG.srate)
    fs = double(EEG.srate);
end
if (isempty(fs) || ~isfinite(fs) || fs <= 0) && ~isempty(time_ms) && numel(time_ms) > 1
    dt = median(diff(double(time_ms(:))));
    if dt > 0
        fs = 1000/dt;
    end
end
if isempty(fs) || ~isfinite(fs) || fs <= 0
    error('psyrat_import:eeglabInvalidSrate', ...
        'Could not determine a valid sampling rate from EEGLAB source.');
end
end

function labels = localChannelLabels(EEG,nChan)
labels = cell(nChan,1);
if isfield(EEG,'chanlocs') && numel(EEG.chanlocs) >= nChan
    for i = 1:nChan
        if isfield(EEG.chanlocs(i),'labels') && ~isempty(EEG.chanlocs(i).labels)
            labels{i} = char(string(EEG.chanlocs(i).labels));
        else
            labels{i} = sprintf('Ch%d',i);
        end
    end
else
    for i = 1:nChan
        labels{i} = sprintf('Ch%d',i);
    end
end
end

function [sub,event,occ,grp,sess,nUnresolved] = localTrialMeta(EEG,nTrials,filepath)
[~,fname] = fileparts(filepath);

if isfield(EEG,'subject') && ~isempty(EEG.subject)
    sid = char(string(EEG.subject));
else
    sid = fname;
end
sub = repmat({sid},nTrials,1);

event = repmat({''},nTrials,1);
nUnresolved = 0;
if isfield(EEG,'epoch') && numel(EEG.epoch) >= nTrials
    % Resolve the time-locking event per epoch (eventlatency closest to zero),
    % not the first listed event; see psyrat_eeglab_epoch_events.
    [event,nUnresolved] = psyrat_eeglab_epoch_events(EEG.epoch,nTrials);
elseif isfield(EEG,'event') && ~isempty(EEG.event)
    for i = 1:min(nTrials,numel(EEG.event))
        if isfield(EEG.event(i),'type')
            event{i} = char(string(EEG.event(i).type));
        end
    end
end

if isfield(EEG,'group') && ~isempty(EEG.group)
    gid = char(string(EEG.group));
else
    gid = '';
end
grp = repmat({gid},nTrials,1);

if isfield(EEG,'session') && ~isempty(EEG.session)
    oid = char(string(EEG.session));
else
    oid = '';
end
occ = repmat({oid},nTrials,1);
sess = occ;
end
