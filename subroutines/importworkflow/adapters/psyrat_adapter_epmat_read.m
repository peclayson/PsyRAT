function out = psyrat_adapter_epmat_read(filepath,varargin)
%Read EP Toolkit exported MATLAB struct into canonical contribution struct.

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

% varargin accepted for interface parity.
if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

out = localEmptyOut('ep_toolkit_mat');
s = load('-mat',filepath);
EP = localFindEPStruct(s);

if isempty(EP)
    error('psyrat_import:epMatInvalid', ...
        'Could not locate an EP Toolkit-like struct in %s', filepath);
end

if ~isfield(EP,'data') || isempty(EP.data)
    error('psyrat_import:epMatNoData', ...
        'EP Toolkit struct did not contain data field in %s', filepath);
end

data = double(EP.data);
if ndims(data) < 2
    error('psyrat_import:epMatShape', ...
        'EP Toolkit data must include channel and time dimensions.');
end

sz = size(data);
while numel(sz) < 7
    sz(end+1) = 1; %#ok<AGROW>
end

nChan = sz(1);
nSamp = sz(2);
nWave = sz(3);
nSub = sz(4);

if sz(5) > 1 || sz(6) > 1 || sz(7) > 1
    out.warnings{end+1} = ['EP data include extra dimensions (factor/rel/frequency). ',...
        'Only the first factor/rel/frequency slice is used for v1 ERP scoring import.'];
end

% Use first factor/rel/frequency slice for ERP scoring workflow.
core = data(:,:,:,:,1,1,1);

% Determine time axis and fs.
if isfield(EP,'timeNames') && ~isempty(EP.timeNames)
    time_ms = double(EP.timeNames(:).');
else
    time_ms = localTimeFromFs(nSamp,EP);
end

fs_hz = localFs(EP,time_ms);
channel_labels = localChannelLabels(EP,nChan);

isSingleTrial = false;
if isfield(EP,'dataType') && ~isempty(EP.dataType)
    isSingleTrial = strcmpi(char(string(EP.dataType)),'single_trial');
end
if ~isSingleTrial
    isSingleTrial = (nWave*nSub) > 1;
end

if isSingleTrial
    signals = zeros(nWave*nSub,nChan,nSamp);

    subNames = localNames(EP,'subNames',nSub,'sub');
    cellNames = localNames(EP,'cellNames',nWave,'cell');

    subject_id = cell(nWave*nSub,1);
    event_id = cell(nWave*nSub,1);

    idx = 1;
    for iSub = 1:nSub
        for iWave = 1:nWave
            wave = squeeze(core(:,:,iWave,iSub));
            signals(idx,:,:) = wave;
            subject_id{idx} = subNames{iSub};
            event_id{idx} = cellNames{iWave};
            idx = idx + 1;
        end
    end

    trial_table = psyrat_build_trial_table(size(signals,1), ...
        'subject_id',subject_id, ...
        'event_id',event_id, ...
        'source_ref',filepath);

    out.single_trial.available = true;
    out.single_trial.signals = signals;
    out.single_trial.time_ms = time_ms;
    out.single_trial.fs_hz = fs_hz;
    out.single_trial.amplitude_unit = 'uV';
    out.single_trial.channel_labels = channel_labels;
    out.single_trial.trial_table = trial_table;
    out.single_trial.data_hash = psyrat_hash_struct(struct( ...
        'size',size(signals), ...
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
    out.average.available = true;
    out.average.signals = core;
    out.average.time_ms = time_ms;
    out.average.fs_hz = fs_hz;
    out.average.channel_labels = channel_labels;
    out.average.source = 'imported_average';
    out.average.data_hash = psyrat_hash_struct(struct('size',size(core),'time_ms',time_ms));

    out.warnings{end+1} = ['EP Toolkit source appears to be averaged/reference data only. ',...
        'Single-trial data are required for PsyRAT scoring and handoff.'];
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

function EP = localFindEPStruct(s)
EP = [];
if isfield(s,'EPdata') && isstruct(s.EPdata)
    EP = s.EPdata;
    return;
end

fn = fieldnames(s);
for i = 1:numel(fn)
    cand = s.(fn{i});
    if isstruct(cand) && isfield(cand,'data') && ...
            (isfield(cand,'chanNames') || isfield(cand,'timeNames'))
        EP = cand;
        return;
    end
end
end

function names = localNames(EP,fieldName,n,defaultPrefix)
if isfield(EP,fieldName) && ~isempty(EP.(fieldName))
    vals = EP.(fieldName);
    if ischar(vals)
        vals = cellstr(vals);
    elseif isstring(vals)
        vals = cellstr(vals(:));
    elseif ~iscell(vals)
        vals = cellstr(string(vals(:)));
    else
        vals = cellfun(@(x) char(string(x)),vals(:),'UniformOutput',false);
    end

    if numel(vals) >= n
        names = vals(1:n);
        return;
    end
end

names = cell(n,1);
for i = 1:n
    names{i} = sprintf('%s_%d',defaultPrefix,i);
end
end

function labels = localChannelLabels(EP,nChan)
if isfield(EP,'chanNames') && ~isempty(EP.chanNames)
    vals = EP.chanNames;
    if ischar(vals)
        vals = cellstr(vals);
    elseif isstring(vals)
        vals = cellstr(vals(:));
    elseif ~iscell(vals)
        vals = cellstr(string(vals(:)));
    else
        vals = cellfun(@(x) char(string(x)),vals(:),'UniformOutput',false);
    end
    if numel(vals) >= nChan
        labels = vals(1:nChan);
        return;
    end
end

labels = cell(nChan,1);
for i = 1:nChan
    labels{i} = sprintf('Ch%d',i);
end
end

function time_ms = localTimeFromFs(nSamp,EP)
fs = [];
if isfield(EP,'Fs') && ~isempty(EP.Fs)
    fs = double(EP.Fs);
elseif isfield(EP,'srate') && ~isempty(EP.srate)
    fs = double(EP.srate);
end
if isempty(fs) || fs <= 0
    error('psyrat_import:epNoTime', ...
        'Could not derive EP Toolkit time axis because sampling rate is missing.');
end

time_ms = ((0:(nSamp-1))/fs)*1000;
end

function fs = localFs(EP,time_ms)
if isfield(EP,'Fs') && ~isempty(EP.Fs)
    fs = double(EP.Fs);
elseif isfield(EP,'srate') && ~isempty(EP.srate)
    fs = double(EP.srate);
else
    dt = median(diff(time_ms));
    fs = 1000/dt;
end
end
