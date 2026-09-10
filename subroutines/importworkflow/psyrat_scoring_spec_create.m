function spec = psyrat_scoring_spec_create(varargin)
%Create a validated single-trial ERP scoring specification.
%
% spec = psyrat_scoring_spec_create('name','P3 local','method','local_peak_amplitude',...)

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

spec = struct();
spec.spec_id = psyrat_uuid('spec');
spec.name = '';
spec.method = '';
spec.target_level = 'single_trial';
spec.window_start_ms = [];
spec.window_end_ms = [];
spec.channel_selector = [];
spec.polarity = 'positive';
spec.boundary_policy = 'strict';
spec.tie_break_policy = 'first';
spec.method_params = struct();
spec.fallback_policy = struct('local_peak_not_found','exclude_trial');
spec.created_utc = psyrat_iso8601_utc;
spec.updated_utc = spec.created_utc;
spec.spec_hash = '';

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Provide name/value pairs.');
end

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'spec_id'
            spec.spec_id = char(string(val));
        case 'name'
            spec.name = char(string(val));
        case 'method'
            spec.method = lower(char(string(val)));
        case 'window_start_ms'
            spec.window_start_ms = double(val);
        case 'window_end_ms'
            spec.window_end_ms = double(val);
        case 'channel_selector'
            spec.channel_selector = val;
        case 'polarity'
            spec.polarity = lower(char(string(val)));
        case 'boundary_policy'
            spec.boundary_policy = lower(char(string(val)));
        case 'tie_break_policy'
            spec.tie_break_policy = lower(char(string(val)));
        case 'method_params'
            spec.method_params = val;
        case 'fallback_policy'
            spec.fallback_policy = val;
        otherwise
            error('varargin:unknown','Unknown option ''%s''.',key);
    end
end

if isempty(spec.name)
    spec.name = spec.spec_id;
end

% Apply method-specific parameter defaults before validation/hashing so the
% effective parameters are explicit on the spec and captured in spec_hash.
spec = localApplyMethodDefaults(spec);

localValidate(spec);
spec.spec_hash = localHashInput(spec);

end

function spec = localApplyMethodDefaults(spec)
% fractional_area_latency uses a configurable area fraction; default to the
% conventional 50% area latency when the caller does not specify one.
if strcmp(spec.method,'fractional_area_latency')
    if ~isstruct(spec.method_params)
        spec.method_params = struct();
    end
    if ~isfield(spec.method_params,'area_fraction') || ...
            isempty(spec.method_params.area_fraction)
        spec.method_params.area_fraction = 0.5;
    end
end
end

function localValidate(spec)
methods = psyrat_scoring_methods();   % shared with the GUI picker (single source of truth)

if ~any(strcmp(spec.method,methods))
    error('psyrat_scoring:unknownMethod', ...
        'Unknown scoring method ''%s''.',spec.method);
end

%COMPLEX-LITERAL GUARD (this one and the two below). ~isreal is added to each
%reject condition because nothing else here sees a complex: isfinite is TRUE for
%one, floor() acts on both parts so floor(n) == n holds for an integer-valued
%complex, and every ordering test -- including the window_end <= window_start
%check just below -- reads the REAL part alone.
%
%The window pair is the one that was SILENT. MEASURED against the unpatched file
%(2026-08-12): psyrat_scoring_spec_create BUILT a spec from window_start_ms =
%100+5i, stored the complex value, and folded it into spec_hash. READ, not
%measured: psyrat_project_score's resolver is
%find(time_ms >= start_ms & time_ms <= end_ms), which compares real parts, so
%scoring would proceed on the real-part window; and its out-of-bounds message
%formats the bound with %.3f, which prints the real part alone.
%
%These values arrive from GUI edit-field text via psyrat_scoring_spec_from_inputs
%(psyrat_startimportscore wires the form), so this is a user-reachable path, not
%a programmatic-only one.
if isempty(spec.window_start_ms) || isempty(spec.window_end_ms) || ...
        ~isreal(spec.window_start_ms) || ~isreal(spec.window_end_ms) || ...
        ~isfinite(spec.window_start_ms) || ~isfinite(spec.window_end_ms)
    error('psyrat_scoring:window', ...
        'Scoring window_start_ms and window_end_ms must be finite real numeric values.');
end

if spec.window_end_ms <= spec.window_start_ms
    error('psyrat_scoring:windowOrder', ...
        'Scoring window_end_ms must be greater than window_start_ms.');
end

if isempty(spec.channel_selector)
    error('psyrat_scoring:channels', ...
        'channel_selector is required.');
end

if ~any(strcmp(spec.polarity,{'positive','negative','both'}))
    error('psyrat_scoring:polarity', ...
        'polarity must be positive, negative, or both.');
end

if strcmp(spec.method,'local_peak_amplitude')
    if ~isstruct(spec.method_params) || ~isfield(spec.method_params,'neighborhood_samples')
        error('psyrat_scoring:localPeakParams', ...
            'local_peak_amplitude requires method_params.neighborhood_samples.');
    end
    n = spec.method_params.neighborhood_samples;
    if ~isscalar(n) || ~isreal(n) || ~isfinite(n) || n < 1 || floor(n) ~= n
        error('psyrat_scoring:localPeakParams', ...
            'method_params.neighborhood_samples must be a positive integer.');
    end
end

if strcmp(spec.method,'centroid_latency')
    if ~isstruct(spec.method_params) || ~isfield(spec.method_params,'weight_mode')
        error('psyrat_scoring:centroidParams', ...
            'centroid_latency requires method_params.weight_mode.');
    end
    mode = lower(char(string(spec.method_params.weight_mode)));
    if ~any(strcmp(mode,{'positive_only','negative_only','absolute'}))
        error('psyrat_scoring:centroidParams', ...
            'centroid weight_mode must be positive_only, negative_only, or absolute.');
    end
end

if strcmp(spec.method,'fractional_area_latency')
    % area_fraction is filled in by localApplyMethodDefaults when omitted, so
    % a missing field here indicates a malformed method_params struct.
    if ~isstruct(spec.method_params) || ~isfield(spec.method_params,'area_fraction')
        error('psyrat_scoring:areaLatencyParams', ...
            'fractional_area_latency requires method_params.area_fraction.');
    end
    f = spec.method_params.area_fraction;
    if ~isscalar(f) || ~isreal(f) || ~isfinite(f) || f <= 0 || f >= 1
        error('psyrat_scoring:areaLatencyParams', ...
            'method_params.area_fraction must be a scalar in the open interval (0,1).');
    end
end
end

function h = localHashInput(spec)
in = spec;
in.spec_hash = '';
in.created_utc = '';
in.updated_utc = '';
h = psyrat_hash_struct(in);
end
