function zvec = psyrat_zscore_subjectlevel(datatable,rawcol)
%Standardize a subject-level (between-person) covariate for the dynamic
%reliability model (Rast & Clayson, analysis 11).
%
% zvec = psyrat_zscore_subjectlevel(datatable,rawcol)
%
%The dimension predictors are continuous covariates that are constant within a
%participant (e.g., a questionnaire total such as pswq_total). This helper
%mirrors the brms standardize_subject_level step: it takes one value per id,
%z-scores those values across participants (mean 0, SD 1 between subjects), and
%broadcasts the standardized value back to every row of that participant. The
%returned vector is aligned to the rows of datatable.
%
%Standardizing once across the whole sample (before any group/event chunking)
%keeps z = 0 at the overall mean and puts every group on a comparable
%predictor scale, so the per-group reliability surfaces are directly
%comparable.
%
%Inputs
% datatable - long-format table with at least an 'id' column and the raw
%  covariate column named by rawcol. One row per single-trial score.
% rawcol - char/string name of the raw covariate column to standardize.
%
%Output
% zvec - column vector (height(datatable) x 1) of standardized values, one per
%  row, broadcast from the subject-level value.
%
%Errors if the covariate is not constant within a participant, contains
%missing/non-finite subject-level values, or has a non-positive between-subject
%standard deviation (no variance to standardize). For a multi-group run (a
%'group' column is present), additionally warns -- before the CmdStan run --
%when a group has no between-participant variation in the covariate or fewer
%than two participants, since that group's scale slope is then unidentified.

% Copyright (C) 2016-2026 Peter E. Clayson
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

rawcol = char(string(rawcol));

if ~ismember('id',datatable.Properties.VariableNames)
    error('dimprep:noid',... %Error code and associated error
        'The data table must contain an ''id'' column to standardize %s.',...
        rawcol);
end

if ~ismember(rawcol,datatable.Properties.VariableNames)
    error('dimprep:nocol',... %Error code and associated error
        'The dimension column ''%s'' was not found in the data table.',rawcol);
end

ids = datatable.id;
rawvals = datatable.(rawcol);

if ~isnumeric(rawvals)
    rawvals = str2double(string(rawvals));
end
rawvals = double(rawvals(:));

%collapse to one value per participant (preserving first appearance order).
%gidx maps each row to its participant index, so per-subject grouping is a
%cheap integer compare (gidx == s) rather than an O(N) ismember scan per
%participant, and the final broadcast is a single vectorized index.
[uid,ia,gidx] = unique(ids,'stable');
nsub = numel(uid);
subjvals = nan(nsub,1);

for s = 1:nsub
    vals = rawvals(gidx == s);
    finitevals = vals(isfinite(vals));

    if isempty(finitevals)
        error('dimprep:missing',... %Error code and associated error
            ['Participant ''%s'' has no finite value for dimension ''%s''. ',...
            'Subject-level dimensions must be present for every participant.'],...
            char(string(uid(s))),rawcol);
    end

    %constant-within-participant check (allow a tiny numeric tolerance)
    if (max(finitevals) - min(finitevals)) > ...
            (1e-9 * max(1,abs(mean(finitevals))))
        error('dimprep:notconstant',... %Error code and associated error
            ['Dimension ''%s'' is not constant within participant ''%s''. ',...
            'Dimension predictors must be between-person (one value per ',...
            'participant).'],rawcol,char(string(uid(s))));
    end

    subjvals(s) = finitevals(1);
end

m = mean(subjvals);
s = std(subjvals); %sample SD (N-1), matching the brms stats::sd convention

if ~isfinite(s) || s <= 0
    error('dimprep:novariance',... %Error code and associated error
        ['Dimension ''%s'' has a non-positive between-participant standard ',...
        'deviation, so it cannot be standardized (no between-person ',...
        'variance to model).'],rawcol);
end

%per-GROUP identifiability heads-up for multi-group runs. The covariate is
%standardized once across the whole sample (z = 0 is the overall mean; this is
%deliberate so the per-group surfaces share a predictor scale), but each group
%is fit as its own model. If a group has no between-participant variation in the
%covariate (or fewer than two participants), that group's scale slope is not
%identified. Warn now, before the (long) CmdStan run, rather than only at plot
%time. This is a warning, not an error: the other groups remain estimable and
%the whole-sample standardization itself is valid.
if ismember('group',datatable.Properties.VariableNames)
    subjgrp = datatable.group(ia);   % one group per subject, in uid order
    ugrp = unique(subjgrp,'stable');
    if numel(ugrp) > 1
        for g = 1:numel(ugrp)
            if isnumeric(subjgrp)
                ingrp = subjgrp == ugrp(g);
            else
                ingrp = strcmp(string(subjgrp),string(ugrp(g)));
            end
            gvals = subjvals(ingrp);
            gsd = std(gvals);
            if numel(gvals) < 2 || ~isfinite(gsd) || gsd <= 0
                warning('dimprep:groupnovariance',... %Warning code + message
                    ['Dimension ''%s'' has no between-participant variation '...
                    'within group ''%s'' (or fewer than two participants), so '...
                    'that group''s scale slope is not identified. Its dynamic '...
                    'surface over ''%s'' will be shown across a nominal range. '...
                    'Consider dropping that group or the dimension.'],...
                    rawcol,char(string(ugrp(g))),rawcol);
            end
        end
    end
end

subjz = (subjvals - m) ./ s;

%broadcast the per-subject z back to every row, aligned to datatable order, via
%the participant index computed above (no per-row map lookup).
zvec = subjz(gidx);

end
