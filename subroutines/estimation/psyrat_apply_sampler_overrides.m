function outargs = psyrat_apply_sampler_overrides(args)
%PSYRAT_APPLY_SAMPLER_OVERRIDES Fold optional NUTS overrides into CmdStan control.
%   OUTARGS = PSYRAT_APPLY_SAMPLER_OVERRIDES(ARGS) scans the name/value cell ARGS
%   for optional 'adapt_delta' and 'max_treedepth' overrides and folds them into
%   the CmdStan 'control' struct (creating one when the calling analysis passed
%   none), then removes the bare override pairs. 'adapt_delta' maps to
%   MatlabStan's control 'delta' (the NUTS target acceptance rate) and
%   'max_treedepth' maps to the NUTS 'max_depth'; both cannot be passed to stan()
%   as top-level arguments, which is why they must live inside 'control'.
%
%   If neither override is present, ARGS is returned unchanged (byte-identical),
%   so default estimation runs are unaffected. Empty override values are treated
%   as "not supplied". This is deliberately a standalone function (rather than a
%   local subfunction of psyrat_computevarcomp) so the merge logic can be unit
%   tested without CmdStan.
%
%   See also PSYRAT_COMPUTEVARCOMP.

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

outargs = args;

[has_ad, ad] = local_get_namevalue(outargs, 'adapt_delta');
[has_mt, mt] = local_get_namevalue(outargs, 'max_treedepth');

%no overrides supplied -> leave the argument list exactly as received
if (~has_ad || isempty(ad)) && (~has_mt || isempty(mt))
    outargs = local_remove_namevalue(outargs, 'adapt_delta');
    outargs = local_remove_namevalue(outargs, 'max_treedepth');
    return;
end

%start from the existing control struct when one was passed, else a fresh struct
[has_ctrl, control] = local_get_namevalue(outargs, 'control');
if ~has_ctrl || isempty(control) || ~isstruct(control)
    control = struct;
end

if has_ad && ~isempty(ad)
    control.adapt_delta = ad;
end
if has_mt && ~isempty(mt)
    control.max_treedepth = mt;
end

%drop the bare override pairs and any prior control, then re-append the merged
%control once so the result stays an even-length name/value list with a single
%'control' entry
outargs = local_remove_namevalue(outargs, 'adapt_delta');
outargs = local_remove_namevalue(outargs, 'max_treedepth');
outargs = local_remove_namevalue(outargs, 'control');
outargs = [outargs, {'control', control}];
end

function [found, value] = local_get_namevalue(args, name)
%Return the value paired with NAME (case-insensitive) in the name/value cell ARGS.
found = false;
value = [];
if isempty(args) || mod(numel(args), 2) ~= 0
    return;
end
for k = 1:2:numel(args)
    key = args{k};
    if (ischar(key) || (isa(key, 'string') && isscalar(key))) && ...
            strcmpi(strtrim(char(key)), name)
        found = true;
        value = args{k+1};
        return;
    end
end
end

function outargs = local_remove_namevalue(args, name)
%Return ARGS with every NAME (case-insensitive) name/value pair removed.
if isempty(args) || mod(numel(args), 2) ~= 0
    outargs = args;
    return;
end
outargs = {};
for k = 1:2:numel(args)
    key = args{k};
    drop = (ischar(key) || (isa(key, 'string') && isscalar(key))) && ...
        strcmpi(strtrim(char(key)), name);
    if ~drop
        outargs{end+1} = args{k};   %#ok<AGROW>
        outargs{end+1} = args{k+1}; %#ok<AGROW>
    end
end
end
