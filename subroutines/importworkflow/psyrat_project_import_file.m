function [project, import_report] = psyrat_project_import_file(varargin)
%Import one ERP source file into a PsyRAT import/scoring project.
%
% [project, import_report] = psyrat_project_import_file('project',project,'file',path)

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

opts = struct('project',[],'file','','freeze_copy',false, ...
    'frozen_copy_dir','');

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'project'
            opts.project = val;
        case 'file'
            opts.file = char(string(val));
        case 'freeze_copy'
            opts.freeze_copy = logical(val);
        case 'frozen_copy_dir'
            opts.frozen_copy_dir = char(string(val));
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(opts.project)
    error('psyrat_import:missingProject','Input ''project'' is required.');
end
if isempty(opts.file)
    error('psyrat_import:missingFile','Input ''file'' is required.');
end

project = opts.project;
filepath = localAbsPath(opts.file);

if exist(filepath,'file') ~= 2
    error('psyrat_import:fileMissing','Source file was not found: %s',filepath);
end

det = psyrat_import_detect_format(filepath);
if ~det.ok
    error('psyrat_import:unsupportedFormat', ...
        ['No supported import adapter recognized this source. ',...
        'How to fix: provide EEGLAB (.set), ERPLAB (.erp), or EP Toolkit exported .mat/.ept.']);
end

adapters = psyrat_import_registry();
idx = find(strcmp({adapters.name},det.adapter),1,'first');
if isempty(idx)
    error('psyrat_import:adapterMissing', ...
        'Detected adapter ''%s'' is not registered.',det.adapter);
end
adapter = adapters(idx);

contrib = adapter.read_fn(filepath);
[source_hash, finfo] = psyrat_hash_file(filepath);
source_id = psyrat_uuid('src');

manifest = struct();
manifest.source_id = source_id;
manifest.source_path = filepath;
manifest.source_hash = source_hash;
manifest.source_size_bytes = finfo.bytes;
manifest.source_mtime = finfo.date;
manifest.detected_format = det.display_name;
manifest.adapter = adapter.name;
manifest.adapter_version = adapter.version;
manifest.detect_confidence = det.confidence;
manifest.detect_reason = det.reason;
manifest.import_utc = psyrat_iso8601_utc;
manifest.warnings = contrib.warnings;

ref = struct();
ref.source_id = source_id;
ref.source_path = filepath;
ref.source_hash = source_hash;
ref.frozen_copy_path = '';

if opts.freeze_copy
    if isempty(opts.frozen_copy_dir)
        contrib.warnings{end+1} = ['freeze_copy was requested but frozen_copy_dir was empty. ',...
            'Source is stored by reference only.'];
    else
        if exist(opts.frozen_copy_dir,'dir') ~= 7
            mkdir(opts.frozen_copy_dir);
        end
        [~,name,ext] = fileparts(filepath);
        freezeTarget = fullfile(opts.frozen_copy_dir,[name '_' source_hash(1:12) ext]);
        copyfile(filepath,freezeTarget,'f');
        ref.frozen_copy_path = freezeTarget;
    end
end

project = localAppendSourceManifest(project,manifest,ref);
project = localMergeCanonical(project,contrib,source_hash,source_id);
project.canonical.warning_flags.averaged_not_used_for_scoring_handoff = ...
    project.canonical.has_average;

if project.canonical.has_average
    warn = struct('id','psyrat_import:averageReferenceOnly', ...
        'what',['Averaged ERP data are present. Scoring/handoff use single-trial data only.'], ...
        'howto','Import single-trial data for scoring and process-data handoff.', ...
        'message',['What happened: Averages are retained only for reference. ',...
        'How to fix: use single-trial imports for scoring/handoff.']);
    project.validation_reports = localAppendStruct(project.validation_reports,warn);
end

project.updated_utc = psyrat_iso8601_utc;
project.audit_log = localAppendAudit(project.audit_log,'source_imported',struct( ...
    'source_id',source_id, ...
    'source_path',filepath, ...
    'adapter',adapter.name, ...
    'single_trial_available',contrib.single_trial.available, ...
    'average_available',contrib.average.available, ...
    'warnings',{contrib.warnings}));

import_report = struct();
import_report.source_id = source_id;
import_report.detect = det;
import_report.single_trial_available = contrib.single_trial.available;
import_report.average_available = contrib.average.available;
import_report.warnings = contrib.warnings;

end

function p = localAbsPath(in)
if ispc
    isAbs = ~isempty(regexp(in,'^[A-Za-z]:[\\/]','once'));
else
    isAbs = startsWith(in,'/');
end
if isAbs
    p = in;
else
    p = fullfile(pwd,in);
end
end

function project = localAppendSourceManifest(project,manifest,ref)
project.source_manifest = localAppendStruct(project.source_manifest,manifest);
project.raw_source_refs = localAppendStruct(project.raw_source_refs,ref);
end

function out = localAppendStruct(in,entry)
if isempty(in)
    out = entry;
else
    out = in;
    out(end+1) = entry;
end
end

function project = localMergeCanonical(project,contrib,source_hash,source_id)
canonical = project.canonical;

if contrib.single_trial.available
    st = contrib.single_trial;
    if ~canonical.has_single_trial
        canonical.single_trial.signals = st.signals;
        canonical.single_trial.time_ms = st.time_ms;
        canonical.single_trial.fs_hz = st.fs_hz;
        canonical.single_trial.amplitude_unit = st.amplitude_unit;
        canonical.single_trial.channel_labels = st.channel_labels;
        canonical.single_trial.data_hash = st.data_hash;
        canonical.trial_table = st.trial_table;
        canonical.has_single_trial = true;
    else
        localRequireCompatibleSingleTrial(canonical.single_trial,st);
        canonical.single_trial.signals = cat(1,canonical.single_trial.signals,st.signals);
        canonical.trial_table = [canonical.trial_table; st.trial_table];
        canonical.single_trial.data_hash = psyrat_hash_struct(struct( ...
            'prev',canonical.single_trial.data_hash, ...
            'source_hash',source_hash, ...
            'n_trials',size(canonical.single_trial.signals,1)));
    end
end

if contrib.average.available
    canonical.has_average = true;
    if isempty(canonical.averaged_reference.signals)
        canonical.averaged_reference.signals = contrib.average.signals;
        canonical.averaged_reference.time_ms = contrib.average.time_ms;
        canonical.averaged_reference.fs_hz = contrib.average.fs_hz;
        canonical.averaged_reference.channel_labels = contrib.average.channel_labels;
        canonical.averaged_reference.source = contrib.average.source;
        canonical.averaged_reference.data_hash = contrib.average.data_hash;
        canonical.averaged_reference_origin = contrib.average.source;
    else
        if ~localAverageCompatible(canonical.averaged_reference,contrib.average)
            project.validation_reports = localAppendStruct(project.validation_reports, ...
                struct('id','psyrat_import:averageMismatch', ...
                'what','Averaged reference data differ across imported sources.', ...
                'howto','Review source-level average differences in audit logs.', ...
                'message','What happened: imported averages are not shape-compatible; first average retained as reference.'));
        end
    end
end

project.canonical = canonical;
project.audit_log = localAppendAudit(project.audit_log,'canonical_updated',struct( ...
    'source_id',source_id, ...
    'has_single_trial',canonical.has_single_trial, ...
    'has_average',canonical.has_average, ...
    'n_trials',height(canonical.trial_table)));
end

function localRequireCompatibleSingleTrial(base,add)
if numel(base.time_ms) ~= numel(add.time_ms) || any(abs(base.time_ms(:)-add.time_ms(:)) > 1e-9)
    error('psyrat_import:timeAxisMismatch', ...
        ['Imported source time axis does not match existing canonical single-trial ',...
        'time axis. This import cannot be merged safely.']);
end

if abs(base.fs_hz - add.fs_hz) > 1e-9
    error('psyrat_import:srateMismatch', ...
        'Imported source sampling rate does not match existing canonical data.');
end

baseLabels = lower(string(base.channel_labels(:)));
addLabels = lower(string(add.channel_labels(:)));
if numel(baseLabels) ~= numel(addLabels) || any(baseLabels ~= addLabels)
    error('psyrat_import:channelMismatch', ...
        'Imported source channel labels do not match existing canonical data.');
end
end

function tf = localAverageCompatible(base,add)
tf = true;
if isempty(base.time_ms) || isempty(add.time_ms)
    return;
end

if numel(base.time_ms) ~= numel(add.time_ms)
    tf = false;
    return;
end

if any(abs(base.time_ms(:)-add.time_ms(:)) > 1e-9)
    tf = false;
    return;
end

if abs(base.fs_hz - add.fs_hz) > 1e-9
    tf = false;
end
end

function out = localAppendAudit(in,action,details)
entry = struct('timestamp_utc',psyrat_iso8601_utc, ...
    'action',action, ...
    'details',details);
if isempty(in)
    out = entry;
else
    out = in;
    out(end+1) = entry;
end
end
