function spec = psyrat_scoring_spec_from_inputs(method,in)
%Build a validated scoring spec from raw (string) GUI inputs.
%
% spec = psyrat_scoring_spec_from_inputs(method,in)
%
% method : scoring-method name (see psyrat_scoring_methods).
% in     : struct of char/string fields collected from the GUI form. Recognized
%          fields (only those relevant to the chosen method are used):
%            name, window_start_ms, window_end_ms, channel, polarity,
%            neighborhood, weight_mode, area_fraction
%
% Centralizing this parsing (instead of inlining it in the GUI callback) keeps
% the string-to-spec logic pure and unit-testable; the GUI only collects the
% strings and calls this. All scientific validation is delegated to
% psyrat_scoring_spec_create.

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

method = lower(char(string(method)));

name = localGet(in,'name','');
w1 = str2double(localGet(in,'window_start_ms',''));
w2 = str2double(localGet(in,'window_end_ms',''));
chsel = localParseChannelSelector(localGet(in,'channel',''));
pol = lower(strtrim(localGet(in,'polarity','positive')));

% Method-specific parameters. Only the parameter(s) the method actually uses
% are read; anything missing is left to psyrat_scoring_spec_create to default
% or to reject.
params = struct();
switch method
    case 'local_peak_amplitude'
        n = round(str2double(localGet(in,'neighborhood','')));
        %COMPLEX-LITERAL GUARD. ~isreal LEADS because round() acts on BOTH
        %parts of a complex (so it stays complex), isfinite is TRUE for one,
        %and n < 1 reads the real part alone. Complex text therefore falls back
        %to the conventional default, exactly as 'abc' already does.
        if ~isreal(n) || ~isfinite(n) || n < 1
            n = 3;   % conventional default neighborhood
        end
        params.neighborhood_samples = n;
    case 'centroid_latency'
        wm = lower(strtrim(localGet(in,'weight_mode','absolute')));
        if isempty(wm)
            wm = 'absolute';
        end
        params.weight_mode = wm;
    case 'fractional_area_latency'
        af = str2double(localGet(in,'area_fraction',''));
        %COMPLEX-LITERAL GUARD: isreal() leads the accept clause, so complex
        %text leaves the field unset and takes the 0.5 default below, exactly
        %as 'abc' already does.
        if isreal(af) && isfinite(af)
            params.area_fraction = af;
        end
        % If left blank, psyrat_scoring_spec_create fills the 0.5 default.
end

spec = psyrat_scoring_spec_create('name',name,'method',method, ...
    'window_start_ms',w1,'window_end_ms',w2, ...
    'channel_selector',chsel,'polarity',pol, ...
    'method_params',params);
end

function val = localGet(in,field,default)
if isstruct(in) && isfield(in,field) && ~isempty(in.(field))
    val = char(string(in.(field)));
else
    val = default;
end
end

function chsel = localParseChannelSelector(chraw)
%Parse channel selector text into a numeric index vector or cellstr labels.
%
% '1'        -> 1
% '1,2,3'    -> [1 2 3]
% 'Fz,Cz,Pz' -> {'Fz','Cz','Pz'} (each scored separately, then averaged)

chraw = strtrim(char(string(chraw)));
if isempty(chraw)
    error('psyrat_scoring:channels','Channel selector cannot be empty.');
end

parts = regexp(chraw,'[,;\s]+','split');
parts = parts(~cellfun(@isempty,parts));
if isempty(parts)
    error('psyrat_scoring:channels','Channel selector cannot be empty.');
end

nums = str2double(parts);
%COMPLEX-LITERAL GUARD. isreal() joins the numeric-branch test: isfinite is TRUE
%for a complex, round() acts on both parts, and any(nums < 1) reads the real
%part, so '1+2i' would have been carried through as a complex channel INDEX --
%and psyrat_scoring_spec_create only checks that channel_selector is non-empty.
%Falling through to the else-branch instead treats the text as channel LABELS,
%which is the same thing any other unparseable channel text does.
if all(isreal(nums)) && all(isfinite(nums))
    nums = round(nums(:)');
    if any(nums < 1)
        error('psyrat_scoring:channels', ...
            'Numeric channel indices must be >= 1.');
    end
    if numel(nums) == 1
        chsel = nums(1);
    else
        chsel = unique(nums,'stable');
    end
else
    chsel = parts(:)';
end
end
