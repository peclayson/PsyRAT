function tbl = psyrat_build_trial_table(nTrials,varargin)
%Build canonical trial metadata table.
%
% Required output columns:
% trial_uid, subject_id, event_id, occasion_id, group_id, session_id,
% source_ref, include_flag, exclude_reason

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

meta = struct();
meta.subject_id = repmat({''},nTrials,1);
meta.event_id = repmat({''},nTrials,1);
meta.occasion_id = repmat({''},nTrials,1);
meta.group_id = repmat({''},nTrials,1);
meta.session_id = repmat({''},nTrials,1);
meta.source_ref = repmat({''},nTrials,1);
meta.include_flag = true(nTrials,1);
meta.exclude_reason = repmat({''},nTrials,1);

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete. Use name/value pairs.');
end

for i = 1:2:length(varargin)
    key = lower(char(string(varargin{i})));
    val = varargin{i+1};
    switch key
        case 'subject_id'
            meta.subject_id = localExpand(val,nTrials);
        case 'event_id'
            meta.event_id = localExpand(val,nTrials);
        case 'occasion_id'
            meta.occasion_id = localExpand(val,nTrials);
        case 'group_id'
            meta.group_id = localExpand(val,nTrials);
        case 'session_id'
            meta.session_id = localExpand(val,nTrials);
        case 'source_ref'
            meta.source_ref = localExpand(val,nTrials);
        case 'include_flag'
            if isscalar(val)
                meta.include_flag = repmat(logical(val),nTrials,1);
            else
                meta.include_flag = logical(val(:));
            end
        case 'exclude_reason'
            meta.exclude_reason = localExpand(val,nTrials);
        otherwise
            error('varargin:unknown','Unknown trial metadata key ''%s''.',key);
    end
end

trial_uid = cell(nTrials,1);
for i = 1:nTrials
    trial_uid{i} = psyrat_uuid('trial');
end

tbl = table(trial_uid,meta.subject_id,meta.event_id,meta.occasion_id, ...
    meta.group_id,meta.session_id,meta.source_ref,meta.include_flag, ...
    meta.exclude_reason, ...
    'VariableNames',{'trial_uid','subject_id','event_id','occasion_id', ...
    'group_id','session_id','source_ref','include_flag','exclude_reason'});

end

function out = localExpand(in,n)
if ischar(in) || (isstring(in) && isscalar(in))
    out = repmat({char(string(in))},n,1);
elseif isstring(in)
    out = cellstr(in(:));
elseif iscell(in)
    out = cellfun(@(x) char(string(x)),in(:),'UniformOutput',false);
else
    out = cellstr(string(in(:)));
end

if numel(out) == 1 && n > 1
    out = repmat(out,n,1);
elseif numel(out) ~= n
    error('psyrat_import:metaLengthMismatch', ...
        'Metadata length (%d) does not match number of trials (%d).', ...
        numel(out),n);
end
end
