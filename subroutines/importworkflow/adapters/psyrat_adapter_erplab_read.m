function out = psyrat_adapter_erplab_read(filepath,varargin)
%Read ERPLAB ERP file into canonical contribution struct.

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

% varargin is accepted for interface parity.
if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

out = localEmptyOut('erplab');
ERP = [];

if exist('pop_loaderp','file') == 2
    try
        ERP = pop_loaderp('filename',filepath);
    catch
        ERP = [];
    end
end

if isempty(ERP)
    s = load('-mat',filepath);
    f = fieldnames(s);
    for i = 1:numel(f)
        cand = s.(f{i});
        if isstruct(cand) && isfield(cand,'bindata')
            ERP = cand;
            break;
        end
    end
end

if isempty(ERP) || ~isstruct(ERP) || ~isfield(ERP,'bindata')
    error('psyrat_import:erplabLoadFailed', ...
        'Could not locate a valid ERPLAB ERP struct with bindata in %s', filepath);
end

bindata = double(ERP.bindata);
if ndims(bindata) ~= 3
    error('psyrat_import:erplabShape', ...
        'Expected ERPLAB bindata shape channel x sample x bin in %s', filepath);
end

out.average.available = true;
out.average.signals = bindata;

if isfield(ERP,'times') && ~isempty(ERP.times)
    out.average.time_ms = double(ERP.times(:).');
else
    out.average.time_ms = localTimeFromSrate(size(bindata,2),ERP);
end

if isfield(ERP,'srate') && ~isempty(ERP.srate)
    out.average.fs_hz = double(ERP.srate);
else
    dt = median(diff(out.average.time_ms));
    out.average.fs_hz = 1000/dt;
end

out.average.channel_labels = localChannelLabels(ERP,size(bindata,1));
out.average.source = 'imported_average';
out.average.data_hash = psyrat_hash_struct(struct( ...
    'size',size(bindata), ...
    'time_ms',out.average.time_ms, ...
    'labels',{out.average.channel_labels}));

out.warnings{end+1} = ['ERPLAB ERP files store averaged ERP waveforms. ',...
    'Single-trial data are required for PsyRAT scoring and process-data handoff, ',...
    'so this source is reference-only unless matched single-trial data are imported.'];

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

function labels = localChannelLabels(ERP,nChan)
labels = cell(nChan,1);
if isfield(ERP,'chanlocs') && numel(ERP.chanlocs) >= nChan
    for i = 1:nChan
        if isfield(ERP.chanlocs(i),'labels') && ~isempty(ERP.chanlocs(i).labels)
            labels{i} = char(string(ERP.chanlocs(i).labels));
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

function time_ms = localTimeFromSrate(nSamp,ERP)
if ~isfield(ERP,'srate') || isempty(ERP.srate)
    error('psyrat_import:erplabNoTime', ...
        'Could not derive ERPLAB time axis because srate is missing.');
end
if isfield(ERP,'xmin') && ~isempty(ERP.xmin)
    start_ms = double(ERP.xmin)*1000;
else
    start_ms = 0;
end
time_ms = start_ms + ((0:(nSamp-1))/double(ERP.srate))*1000;
end
