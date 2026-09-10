function [project, run_report] = psyrat_project_score(varargin)
%Compute single-trial ERP component scores for one scoring spec.
%
% [project, run_report] = psyrat_project_score('project',project,'spec_id','spec_x')
%
%Required inputs (name/value pairs)
% project - import-workflow project struct (see psyrat_project_new and
%  psyrat_project_import_file).
% spec_id - id of the scoring spec to run (see
%  psyrat_project_add_scoring_spec).
%
%Optional inputs
% trial_filter - struct of trial-selection criteria applied before scoring;
%  trials failing the filter are excluded from the run (default: empty
%  struct, no filtering).
%
%Outputs
% project - project struct updated with the scored trials for this run.
% run_report - summary struct with fields run_id, spec_id, n_trials,
%  n_failed, n_included, and a human-readable message.

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

opts = struct('project',[],'spec_id','','trial_filter',struct());

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
        case 'trial_filter'
            opts.trial_filter = val;
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(opts.project)
    error('psyrat_scoring:missingProject','Input ''project'' is required.');
end

project = opts.project;

if ~project.canonical.has_single_trial
    error('psyrat_scoring:singleTrialRequired', ...
        ['Single-trial data are required for scoring. ',...
        'Import single-trial ERP sources before running psyrat_project_score.']);
end

if isempty(project.scoring_specs)
    error('psyrat_scoring:noSpecs','No scoring specs are defined in this project.');
end

if isempty(opts.spec_id)
    opts.spec_id = project.ui_state.active_spec_id;
end

specIdx = find(strcmp({project.scoring_specs.spec_id},opts.spec_id),1,'first');
if isempty(specIdx)
    error('psyrat_scoring:specNotFound', ...
        'Scoring spec was not found: %s', opts.spec_id);
end
spec = project.scoring_specs(specIdx);

signals = project.canonical.single_trial.signals;
time_ms = project.canonical.single_trial.time_ms(:);
trial_table = project.canonical.trial_table;
channel_labels = project.canonical.single_trial.channel_labels;

if height(trial_table) ~= size(signals,1)
    error('psyrat_scoring:trialMismatch', ...
        'trial_table row count does not match number of single-trial signals.');
end

mask = localBuildTrialMask(trial_table,opts.trial_filter);
mask = mask & trial_table.include_flag;
trialIdx = find(mask);
if isempty(trialIdx)
    error('psyrat_scoring:noTrialsSelected', ...
        'No trials matched the selection and include_flag filters.');
end

channelIdx = localResolveChannels(spec.channel_selector,project.canonical.single_trial.channel_labels);
winIdx = localResolveWindows(spec,time_ms);

run_id = psyrat_uuid('run');
now_utc = psyrat_iso8601_utc;
rows = repmat(localEmptyRow(),numel(trialIdx),1);
nEdgeTrials = 0;   % scored trials whose peak landed on the measurement-window boundary

for i = 1:numel(trialIdx)
    idx = trialIdx(i);
    chVals = NaN(1,numel(channelIdx));
    chLats = NaN(1,numel(channelIdx));
    chFailed = false(1,numel(channelIdx));
    chReasons = repmat({''},1,numel(channelIdx));
    chEdge = false(1,numel(channelIdx));

    for c = 1:numel(channelIdx)
        wave = squeeze(signals(idx,channelIdx(c),:));
        wave = wave(:);

        % Single trials are assumed to be already baseline-adjusted upstream;
        % the scorer performs no baseline correction.
        [ok,val,lat,reason,isEdge] = localScoreOne(wave,time_ms,winIdx,spec);
        if ok
            chVals(c) = val;
            chLats(c) = lat;
            chEdge(c) = isEdge;
        else
            chFailed(c) = true;
            chReasons{c} = reason;
        end
    end

    % Channel-failure policy (J.5, owner-kept 2026-06-26): if ANY channel in a
    % multi-channel ROI fails to score, the whole trial is excluded rather than
    % averaging over the surviving channels. This keeps a fixed, consistent ROI
    % across every scored trial; averaging a variable channel subset trial-to-
    % trial would introduce inconsistent spatial sampling (nuisance variance).
    if any(chFailed)
        ok = false;
        val = NaN;
        lat = NaN;
        reason = localChannelFailureReason(channelIdx,chFailed,chReasons,channel_labels);
    else
        ok = true;
        val = mean(chVals,'omitnan');
        lat = localMeanFinite(chLats);
        reason = '';
    end

    if ok && any(chEdge)
        nEdgeTrials = nEdgeTrials + 1;
    end

    r = localEmptyRow();
    r.run_id = run_id;
    r.spec_id = spec.spec_id;
    r.trial_uid = trial_table.trial_uid{idx};
    r.value = val;
    r.latency_ms = lat;
    if ok
        r.status = 'ok';
        r.status_reason = '';
        r.excluded_from_handoff = false;
    else
        r.status = 'failed';
        r.status_reason = reason;
        r.excluded_from_handoff = true;
    end
    r.spec_hash = spec.spec_hash;
    r.data_hash = project.canonical.single_trial.data_hash;
    r.computed_utc = now_utc;
    r.is_stale = false;
    rows(i) = r;
end

if nEdgeTrials > 0
    warning('psyrat_scoring:edgePeak', ...
        ['%d of %d trials had the peak at the measurement-window ' ...
         'boundary (method ''%s''). The scoring window may be too narrow; ' ...
         'the reported peak may not be a true local extremum.'], ...
        nEdgeTrials, numel(trialIdx), spec.method);
end

newRows = struct2table(rows,'AsArray',true);

if isempty(project.scoring_results)
    project.scoring_results = newRows;
else
    if ~any(strcmp(project.scoring_results.Properties.VariableNames,'is_stale'))
        project.scoring_results.is_stale = false(height(project.scoring_results),1);
    end
    oldSpec = strcmp(project.scoring_results.spec_id,spec.spec_id) & ...
        ~project.scoring_results.is_stale;
    project.scoring_results.is_stale(oldSpec) = true;
    project.scoring_results = [project.scoring_results; newRows];
end

project.ui_state.active_spec_id = spec.spec_id;
project.ui_state.active_run_id = run_id;
project.updated_utc = now_utc;

nFailed = sum(newRows.excluded_from_handoff);
project.audit_log = localAppendAudit(project.audit_log,'scoring_run',struct( ...
    'run_id',run_id, ...
    'spec_id',spec.spec_id, ...
    'method',spec.method, ...
    'n_trials',height(newRows), ...
    'n_failed',nFailed));

run_report = struct();
run_report.run_id = run_id;
run_report.spec_id = spec.spec_id;
run_report.n_trials = height(newRows);
run_report.n_failed = nFailed;
run_report.n_included = height(newRows) - nFailed;
run_report.message = sprintf(['Scoring run completed: %d included, %d excluded. ',...
    'Excluded trials are omitted from final handoff dataset.'], ...
    run_report.n_included, run_report.n_failed);

end

function mask = localBuildTrialMask(tbl,flt)
mask = true(height(tbl),1);
if isempty(flt) || ~isstruct(flt)
    return;
end

fields = {'subject_id','event_id','occasion_id','group_id','session_id'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(flt,f) && ~isempty(flt.(f))
        vals = string(flt.(f));
        mask = mask & ismember(string(tbl.(f)),vals(:));
    end
end
end

function ch = localResolveChannels(selector,labels)
if isnumeric(selector)
    ch = selector(:)';
    if any(ch < 1) || any(ch > numel(labels))
        error('psyrat_scoring:channelIndex', ...
            'channel_selector contained out-of-range channel indices.');
    end
    return;
end

if ischar(selector) || (isstring(selector) && isscalar(selector))
    selector = {char(string(selector))};
elseif isstring(selector)
    selector = cellstr(selector(:));
elseif ~iscell(selector)
    error('psyrat_scoring:channelSelectorType', ...
        'channel_selector must be numeric indices or channel labels.');
end

labelsLower = lower(string(labels(:)));
ch = zeros(1,numel(selector));
for i = 1:numel(selector)
    idx = find(labelsLower == lower(string(selector{i})),1,'first');
    if isempty(idx)
        error('psyrat_scoring:channelNotFound', ...
            'Requested channel label not found: %s', char(string(selector{i})));
    end
    ch(i) = idx;
end
ch = unique(ch,'stable');
end

function winIdx = localResolveWindows(spec,time_ms)
start_ms = spec.window_start_ms;
end_ms = spec.window_end_ms;

if strcmp(spec.boundary_policy,'strict')
    if start_ms < min(time_ms) || end_ms > max(time_ms)
        error('psyrat_scoring:windowOutOfBounds', ...
            'Scoring window [%.3f %.3f] ms falls outside available time axis [%.3f %.3f] ms.', ...
            start_ms,end_ms,min(time_ms),max(time_ms));
    end
end

winIdx = find(time_ms >= start_ms & time_ms <= end_ms);
if isempty(winIdx)
    error('psyrat_scoring:windowEmpty', ...
        'Scoring window did not select any samples.');
end
end

function [ok,val,lat,reason,isEdge] = localScoreOne(wave,time_ms,winIdx,spec)
ok = true;
val = NaN;
lat = NaN;
reason = '';
isEdge = false;   % true when a peak/latency extremum sits on a window boundary

seg = wave(winIdx);
segTimes = time_ms(winIdx);

switch spec.method
    case 'absolute_peak_amplitude'
        [ok,amp,idx,reason] = localResolvePeak(seg,spec.polarity,spec.tie_break_policy);
        if ~ok
            return;
        end
        val = amp;
        lat = segTimes(idx);
        isEdge = localIsEdgeIndex(idx,numel(seg));

    case 'local_peak_amplitude'
        n = spec.method_params.neighborhood_samples;
        [found,amp,idx] = localLocalPeak(seg,spec.polarity,n,spec.tie_break_policy);
        if ~found
            ok = false;
            reason = sprintf('No local peak found in window (neighborhood=%d).',n);
            return;
        end
        val = amp;
        lat = segTimes(idx);

    case 'time_window_mean_amplitude'
        val = mean(seg,'omitnan');
        lat = NaN;

    case 'peak_latency'
        [ok,~,idx,reason] = localResolvePeak(seg,spec.polarity,spec.tie_break_policy);
        if ~ok
            return;
        end
        lat = segTimes(idx);
        val = lat;
        isEdge = localIsEdgeIndex(idx,numel(seg));

    case 'centroid_latency'
        % Amplitude-weighted mean latency (center of mass of the window).
        % Nonstandard for ERP component timing; fractional_area_latency is the
        % more conventional measure and is preferred for new specs.
        mode = lower(char(string(spec.method_params.weight_mode)));
        w = localCentroidWeights(seg,mode);
        den = sum(w,'omitnan');
        if den <= 0 || ~isfinite(den)
            ok = false;
            reason = 'Centroid weights summed to zero in window.';
            return;
        end
        lat = sum(segTimes.*w,'omitnan')/den;
        val = lat;

    case 'integrated_area_amplitude'
        % Time integral (uV*ms) of the polarity-selected waveform over the
        % window (the "area amplitude"). Sign-preserving: positive/both yield a
        % positive integral, negative a negative integral. Distinct from
        % time_window_mean_amplitude, which is in uV. Baseline-sensitive; any
        % baseline correction is applied upstream before this point.
        contrib = localPolaritySegment(seg,spec.polarity);
        val = trapz(segTimes,contrib);
        lat = NaN;

    case 'fractional_area_latency'
        % Latency (ms) at which the cumulative area of the polarity-selected
        % waveform reaches a fraction (method_params.area_fraction) of the
        % total window area. Uses the area magnitude so the crossing is
        % polarity-agnostic in scale, with linear interpolation for sub-sample
        % latency. The conventional 50%-area latency measure for ERP timing.
        % For a sign-stable latency set polarity to the component's polarity
        % (positive/negative); 'both' rectifies and mixes lobes across a zero
        % crossing. A sign-aware biphasic variant is intentionally not provided
        % (J.6, owner-closed 2026-06-26: the standard use is covered by setting
        % polarity, and net-signed-area latency is ill-defined for biphasic).
        frac = spec.method_params.area_fraction;
        contrib = abs(localPolaritySegment(seg,spec.polarity));
        [found,lat] = localFractionalAreaLatency(contrib,segTimes,frac);
        if ~found
            ok = false;
            reason = 'Area was zero or non-finite in window for the selected polarity.';
            return;
        end
        val = lat;

    otherwise
        ok = false;
        reason = ['Unsupported method: ' spec.method];
end

if ~isfinite(val)
    ok = false;
    reason = 'Computed value was non-finite.';
end
end

function m = localMeanFinite(x)
xf = x(isfinite(x));
if isempty(xf)
    m = NaN;
else
    m = mean(xf,'omitnan');
end
end

function reason = localChannelFailureReason(channelIdx,chFailed,chReasons,channel_labels)
parts = {};
failIdx = find(chFailed);
for i = 1:numel(failIdx)
    c = failIdx(i);
    chNum = channelIdx(c);
    chLab = '';
    if iscell(channel_labels) && numel(channel_labels) >= chNum
        chLab = char(string(channel_labels{chNum}));
    end
    r = chReasons{c};
    if isempty(chLab)
        parts{end+1} = sprintf('ch%d: %s',chNum,r); %#ok<AGROW>
    else
        parts{end+1} = sprintf('ch%d (%s): %s',chNum,chLab,r); %#ok<AGROW>
    end
end
reason = ['One or more channels failed scoring; trial excluded. Details: ' strjoin(parts,' | ')];
end

function [amp,idx] = localPeak(seg,polarity,tiePolicy)
if strcmp(polarity,'negative')
    target = min(seg);
    idxAll = find(seg == target);
elseif strcmp(polarity,'both')
    % Largest absolute deviation. Build the full tie set (rather than taking
    % max()'s first index) so the 'both' path honors tie_break_policy the same
    % way the +/- paths do. Rectify once and reuse.
    a = abs(seg);
    target = max(a);
    idxAll = find(a == target);
else
    target = max(seg);
    idxAll = find(seg == target);
end
idx = localTie(idxAll,tiePolicy);
if isnan(idx)
    % An all-NaN (or all-non-finite) window leaves idxAll empty for every
    % polarity, so localTie returns NaN. Signal failure here instead of
    % indexing seg(NaN), and let the caller mark the trial excluded.
    amp = NaN;
    return;
end
amp = seg(idx);
end

function [ok,amp,idx,reason] = localResolvePeak(seg,polarity,tiePolicy)
% Shared peak/latency index resolution for absolute_peak_amplitude and
% peak_latency: locate the extremum, and report a clean failure when an all-NaN
% window leaves no index (so the caller never indexes seg(NaN)).
[amp,idx] = localPeak(seg,polarity,tiePolicy);
if isnan(idx)
    ok = false;
    reason = 'No finite samples in window for peak detection.';
else
    ok = true;
    reason = '';
end
end

function tf = localIsEdgeIndex(idx,n)
% True when the peak index sits on a window boundary (first or last sample),
% which usually signals a too-narrow measurement window.
tf = (idx == 1) || (idx == n);
end

function [found,amp,idx] = localLocalPeak(seg,polarity,n,tiePolicy)
found = false;
amp = NaN;
idx = NaN;

if numel(seg) < (2*n + 1)
    return;
end

cand = false(size(seg));
for i = (1+n):(numel(seg)-n)
    left = seg((i-n):(i-1));
    right = seg((i+1):(i+n));
    if strcmp(polarity,'negative')
        cand(i) = all(seg(i) <= left) && all(seg(i) <= right);
    elseif strcmp(polarity,'both')
        isPos = all(seg(i) >= left) && all(seg(i) >= right);
        isNeg = all(seg(i) <= left) && all(seg(i) <= right);
        cand(i) = isPos || isNeg;
    else
        cand(i) = all(seg(i) >= left) && all(seg(i) >= right);
    end
end

idxAll = find(cand);
if isempty(idxAll)
    return;
end

if strcmp(polarity,'negative')
    [~,k] = min(seg(idxAll));
elseif strcmp(polarity,'both')
    [~,k] = max(abs(seg(idxAll)));
else
    [~,k] = max(seg(idxAll));
end

bestVal = seg(idxAll(k));
equalIdx = idxAll(seg(idxAll)==bestVal);
idx = localTie(equalIdx,tiePolicy);
amp = seg(idx);
found = true;
end

function w = localCentroidWeights(seg,mode)
switch mode
    case 'positive_only'
        w = max(seg,0);
    case 'negative_only'
        w = abs(min(seg,0));
    otherwise
        w = abs(seg);
end
end

function out = localPolaritySegment(seg,polarity)
% Polarity-driven sign selection shared by the area-based methods, matching
% the convention the peak methods already use: positive keeps positive-going
% deflections, negative keeps negative-going deflections, both rectifies.
switch polarity
    case 'negative'
        out = min(seg,0);
    case 'both'
        out = abs(seg);
    otherwise
        out = max(seg,0);
end
end

function [found,lat] = localFractionalAreaLatency(mag,segTimes,frac)
% Latency at which the cumulative trapezoidal area of mag (a non-negative
% area density, uV) reaches frac of the total window area, with linear
% interpolation between the bracketing samples. segTimes is monotonic (ms).
found = false;
lat = NaN;

cum = cumtrapz(segTimes,mag);   % cumulative area at each sample (uV*ms)
total = cum(end);
if ~isfinite(total) || total <= 0
    return;                     % flat/empty window for the selected polarity
end

target = frac*total;

% First sample whose cumulative area reaches the target (cum is
% non-decreasing because mag >= 0, so the first crossing is well defined).
k = find(cum >= target,1,'first');
if isempty(k)
    return;
end

if k == 1
    lat = segTimes(1);
else
    c0 = cum(k-1);
    c1 = cum(k);
    if c1 > c0
        w = (target - c0)/(c1 - c0);
    else
        w = 0;              % degenerate flat step; snap to the bracket start
    end
    lat = segTimes(k-1) + w*(segTimes(k) - segTimes(k-1));
end
found = true;
end

function idx = localTie(idxAll,policy)
if isempty(idxAll)
    idx = NaN;
    return;
end
if strcmp(policy,'last')
    idx = idxAll(end);
else
    idx = idxAll(1);
end
end

function row = localEmptyRow()
row = struct('run_id','','spec_id','','trial_uid','', ...
    'value',NaN,'latency_ms',NaN,'status','','status_reason','', ...
    'excluded_from_handoff',false,'spec_hash','','data_hash','', ...
    'computed_utc','','is_stale',false);
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
