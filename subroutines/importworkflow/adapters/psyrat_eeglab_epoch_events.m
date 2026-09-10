function [event, nUnresolved] = psyrat_eeglab_epoch_events(epochs,nTrials)
%Resolve the time-locking event type for each EEGLAB epoch.
%
% [event, nUnresolved] = psyrat_eeglab_epoch_events(EEG.epoch, nTrials)
%
% EEGLAB stores, per epoch, the types (EEG.epoch(i).eventtype) and latencies
% (EEG.epoch(i).eventlatency) of EVERY event captured inside the epoch window,
% ordered by latency. The time-locking event is the one whose latency is ~0;
% any event that falls in the pre-stimulus baseline appears earlier in the
% list. Taking the first listed event is therefore wrong whenever an epoch
% contains more than the time-locking event (for example a preceding stimulus,
% or a response captured in the baseline), and that mislabeling would feed the
% wrong condition into the 'event' facet of every reliability estimate.
%
% This selects the eventtype whose eventlatency is closest to zero. When an
% epoch holds multiple events but no usable parallel eventlatency is available,
% it falls back to the first listed event (preserving the legacy behavior) and
% counts the epoch in nUnresolved so the caller can warn. Single-event epochs
% are unambiguous and are never counted as unresolved.

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

event = repmat({''},nTrials,1);
nUnresolved = 0;

for i = 1:nTrials
    if ~isfield(epochs(i),'eventtype') || isempty(epochs(i).eventtype)
        continue;
    end
    types = epochs(i).eventtype;

    % A single event may be stored as a plain value or a scalar cell; either
    % way there is only one choice, so no latency lookup is needed.
    if ~iscell(types)
        event{i} = char(string(types));
        continue;
    end
    if numel(types) == 1
        event{i} = char(string(types{1}));
        continue;
    end

    % Multiple events: use the parallel eventlatency to find the time-locking
    % event (latency closest to zero, robust to small non-zero latencies from
    % resampling or rounding).
    lat = localEventLatencies(epochs(i));

    if numel(lat) == numel(types) && any(isfinite(lat))
        [~,idx] = min(abs(lat));
        event{i} = char(string(types{idx}));
    else
        % No usable latency information: keep the legacy first-event behavior
        % and flag the epoch so the caller can surface a warning.
        event{i} = char(string(types{1}));
        nUnresolved = nUnresolved + 1;
    end
end
end

function lat = localEventLatencies(ep)
% Return the epoch's event latencies as a numeric row vector, or [] when the
% field is absent/empty. Non-numeric or empty entries become NaN so a partial
% latency list still disqualifies the epoch via the length/finite checks.
lat = [];
if ~isfield(ep,'eventlatency') || isempty(ep.eventlatency)
    return;
end
raw = ep.eventlatency;
if iscell(raw)
    lat = cellfun(@localScalarOrNaN,raw);
else
    lat = double(raw(:)).';
end
end

function v = localScalarOrNaN(x)
if isempty(x) || ~isnumeric(x)
    v = NaN;
else
    v = double(x(1));
end
end
