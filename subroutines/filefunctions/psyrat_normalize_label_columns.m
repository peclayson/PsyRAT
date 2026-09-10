function datatable = psyrat_normalize_label_columns(datatable,convertnumeric)
%Convert non-cellstr label columns of a data table to cellstr
%
%datatable = psyrat_normalize_label_columns(datatable)
%datatable = psyrat_normalize_label_columns(datatable,convertnumeric)
%
%
%Input
% datatable - long-format data table. The label columns 'id', 'group',
%             'event' and 'time' are normalized when present; every other
%             column (including 'meas', 'dim1', 'dim2' and 'weight') is
%             left untouched.
% convertnumeric - OPTIONAL logical, default false. When false, numeric label
%             columns pass through untouched. When true, NUMERIC label columns
%             are converted with cellstr(string(...)).
%             psyrat_loadfile calls this function TWICE and the two calls
%             differ only in this flag. The FIRST, above its which* row
%             filters, takes the default false: those filters compare a
%             numeric column with == against a numeric selection, so
%             converting first breaks a working path. The SECOND, at the very
%             end of the loader, passes true and owns numeric conversion.
%             psyrat_computevarcomp likewise passes true, because by the time
%             the engine is reached the filters have run and a numeric label
%             column is nothing but a hazard there.
%             LOGICAL label columns are NEVER converted, at either setting --
%             see the numeric block for why.
%
%output
% datatable - the same table with any label column that is not already
%             cellstr converted to a cellstr column vector. A table whose
%             label columns are already cellstr is returned unchanged.
%
%Errors
% labels:missingvalues   - a converted label column contained a missing
%                          value (<missing>, <undefined>, NaT or NaN)
% labels:unsupportedtype - a label column held a type with no text form
% labels:numericlabelcollision - convertnumeric was true and two distinct
%                          numeric label values rendered as the SAME text,
%                          which would silently merge two levels
% labels:convertnumeric  - convertnumeric was not a logical scalar
%

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
%

%WHY THIS EXISTS (finding B12). The toolbox stores label columns as cellstr
%and appends them into cell arrays initialized as {} -- for example
%REL.out.elabels and REL.out.labels in psyrat_computevarcomp. Anything that is
%not a cell fails such an append with 'MATLAB:UnableToConvert'. Those appends
%run only AFTER CmdStan has finished sampling, so the failure arrives at the
%end of a run that may have taken hours: an observed case lost a completed
%four-chain fit after 4 h 20 m. Normalizing at the input boundary keeps the
%type from ever reaching the estimation core.
%
%The types get in easily. readtable(file,'TextType','string') yields string
%columns, and bare readtable AUTO-DETECTS a date-like column as datetime --
%so an occasion column holding session dates arrives as datetime through the
%ordinary file path, with no unusual options set by the user.

if nargin < 2 || isempty(convertnumeric)
    convertnumeric = false;
end

%Validate rather than relying on truthiness. Without this,
%psyrat_normalize_label_columns(tbl,'no') would CONVERT, because a non-empty
%char is truthy -- the opposite of what the caller wrote.
if ~(islogical(convertnumeric) && isscalar(convertnumeric))
    error('labels:convertnumeric',... %Error code and associated error
        ['Input ''convertnumeric'' must be a logical scalar (true or false). ' ...
        'How to fix: pass true to convert numeric label columns, or omit it ' ...
        'to leave them untouched.']);
end

labelcols = {'id','group','event','time'};
varnames  = datatable.Properties.VariableNames;

%Match the canonical lowercase names exactly. psyrat_loadfile builds this
%table with exactly these names, and psyrat_computevarcomp dereferences them
%as datatable.id / .group / .event / .time, which is case-sensitive. Matching
%case-insensitively here would normalize a column the estimation core cannot
%read anyway.
for ii = 1:numel(labelcols)

    thiscol = labelcols{ii};

    %group, event and time are optional, so skip any that this design does
    %not carry rather than dereferencing a column that is not there.
    if ~ismember(thiscol,varnames)
        continue;
    end

    coldata = datatable.(thiscol);

    %ALREADY CELLSTR: nothing to convert, and nothing to check either.
    %
    %The missing check below is deliberately NOT extended here, and the reason
    %is that it would be redundant rather than that it would be unsafe: an
    %empty char '' inside a cellstr column is ALREADY rejected, by
    %psyrat_preflight_validate's table-wide 'varargin:missingcells' (measured
    %on a cellstr group column with '' entries). Duplicating that check here
    %would change which error id the user sees, not whether the run stops.
    %
    %Do not read this as "empty labels are harmless". They are not:
    %psyrat_relsummary.m's psyrat_relsummary_groupdata treats '' as "no group"
    %and skips group subsetting entirely, so an empty group label that did slip
    %through would silently return every row for that group. That is what the
    %preflight guard is protecting, and nothing here should be taken as a
    %reason to weaken it.
    if iscellstr(coldata) %#ok<ISCLSTR> cellfun(@ischar) is not equivalent here
        continue;
    end

    %NUMERIC AND LOGICAL ARE CALL-SITE DEPENDENT, and this one is easy to get
    %wrong in either direction.
    %
    %WHY THE LOADER'S FIRST CALL MUST NOT CONVERT (convertnumeric = false, the
    %default). psyrat_loadfile calls this function twice, and that first call
    %sits ABOVE its which* row filters. That ordering is required, not
    %incidental: the filters pick a branch from the type of the user's
    %SELECTION and compare a numeric column with `dataout.group == spec`
    %(psyrat_loadfile.m, the whichgroups block). Converting a numeric column to
    %text before that point turns a working comparison into an error. The
    %loader's SECOND call, at the end of that function, passes true and does
    %the numeric conversion once the filters are done with it.
    %
    %WHY THE ENGINE -- AND THE LOADER'S SECOND CALL -- MUST CONVERT
    %(convertnumeric = true). By the time either of those call sites runs, the
    %which* filters are done with the numeric column. In the engine's case a
    %numeric label column there is a live defect rather than a working path:
    %the engine's own conversion for cases 2/3/4 used bare num2str, which
    %returns a CHAR MATRIX for a column vector, then linear-indexed it as if it
    %were a cellstr while the data column stayed numeric -- so ismember
    %compared label values against character CODE POINTS and every group mask
    %came back EMPTY (measured 2026-08-12: NG1 = NG2 = 0, for same-width labels
    %as well as mixed-width). Some later cases convert for themselves (case 9,
    %and case 8 via psyrat_compute_diff_sserr), but others store the raw
    %pre-switch unique() output, which for a numeric column is a numeric
    %REL.groups that the downstream gnames{i} consumers cannot brace-index
    %either. Converting here covers both, and the four bare-num2str branches
    %were deleted with this change because they are now unreachable.
    %
    %The spelling is cellstr(string(...)), NOT cellstr(num2str(...)), for ONE
    %verified reason: num2str pads a column vector to a common width (' 1'
    %against '10'), and that leading space changes sort order and breaks string
    %matches against the unpadded labels the rest of the pipeline produces.
    %string() does not pad, and is already the convention at
    %psyrat_computevarcomp's own cellstr(string(unique(...))) sites and in
    %psyrat_preflight_validate.
    %
    %DO NOT claim string() also buys precision -- it does not, and an earlier
    %version of this comment did. MEASURED 2026-08-12: string() uses the same
    %short format as num2str, so string(pi) is "3.1416" and BOTH idioms render
    %0.123456789 and 0.123456111 as '0.12346'. Distinct non-integer labels can
    %therefore collapse into one level under either spelling, which would
    %silently MERGE two groups, events or occasions. That is why the collision
    %check below exists: the round trip is not lossless, so the merge is
    %detected and refused rather than assumed away.
    %
    %THE LOADER NOW USES THIS SPELLING TOO (ruled 2026-08-12, taken 2026-08-14).
    %psyrat_loadfile's four hand-rolled cellstr(num2str(...)) blocks are gone,
    %replaced by a second call to this function with convertnumeric = true, so
    %there is no longer a padded spelling anywhere in the toolbox. That change
    %moves shipped label text and level ORDER, and it was only safe to take
    %after B23: before it, the difference-score cases read event order from a
    %SORTED unique while the Stan model matrix was keyed to a 'stable' one, so
    %order carried a VALUE and switching the spelling would have silently
    %reshuffled which numeric-label datasets were accidentally correct.
    %
    %The missing check has to run BEFORE the numeric conversion, for the same
    %reason it does for the other types below and with a sharper edge here:
    %string(NaN) is the literal "NaN", which preflight's missing-cell guard
    %CANNOT see. Both converting call sites run before preflight, so converting
    %first would turn a caught missing label into an accepted level named NaN.
    %The loader's old numeric blocks DID produce a literal 'NaN' label; routing
    %them through here replaced that long-standing behavior with a
    %labels:missingvalues error, which is an approved, user-visible change.
    %
    %LOGICAL IS DELIBERATELY NOT CONVERTED, even here. string(true) is "true",
    %so a logical column would become {'true';'false'} on this path, where the
    %loader's old num2str idiom would have given {'1';'0'} -- and that idiom
    %never ran on logical anyway, since isnumeric(true) is false. Inventing a
    %third label spelling for a type no shipped path converts is not something
    %this fix is entitled to do, so logical keeps passing through untouched.
    if isnumeric(coldata) && convertnumeric
        local_reject_missing(coldata,thiscol);
        converted = cellstr(string(coldata(:)));
        local_reject_label_collision(coldata,converted,thiscol);
        datatable.(thiscol) = converted;
        continue;
    end
    if isnumeric(coldata) || islogical(coldata)
        continue;
    end

    %Reject missing BEFORE converting, because converting first would disguise
    %what went wrong: cellstr() maps a <missing> string to '' and an
    %<undefined> categorical to the literal '<undefined>', so the problem
    %reappears downstream as an odd-looking level rather than as a missing one.
    %
    %BE PRECISE ABOUT WHAT THIS BUYS, because it is easy to overstate. It is
    %NOT the only thing standing between a missing label and a wrong number:
    %psyrat_preflight_validate already rejects any missing cell table-wide with
    %'varargin:missingcells', and it catches both the raw <missing> and the ''
    %that conversion would produce (measured, both cases). So a missing label
    %never reaches estimation either way.
    %
    %What this adds is a BETTER FAILURE: it fires in psyrat_loadfile, before
    %preflight runs, and it names the offending column and row and says why a
    %label cannot be imputed -- where the preflight message says only that the
    %data include missing cells, without locating them. A usability
    %improvement on an existing guard, not a new correctness guarantee.
    local_reject_missing(coldata,thiscol);

    datatable.(thiscol) = local_to_cellstr(coldata,thiscol);

end

end

function local_reject_label_collision(coldata,converted,colname)
%Error if numeric-to-text conversion merged two distinct label values.
%
%MATLAB's default numeric-to-text format keeps about five significant digits,
%for string() exactly as for num2str (measured 2026-08-12: string(pi) is
%"3.1416"). So a label column holding, say, 0.123456789 and 0.123456111 renders
%both as '0.12346'. unique() would then see ONE level where the data have two,
%silently pooling two groups, events or occasions into a single stratum and
%changing every variance component computed from them.
%
%That is a wrong number rather than a crash, which is the failure mode this
%toolbox least wants to ship, so it is refused here. Integer labels -- the
%overwhelmingly common case -- are unaffected: they round-trip exactly.
%The relation is >=, not ==. Only a DROP in level count is a merge, and that is
%the only thing this guard is entitled to refuse. An == test would also fire if
%the text somehow carried MORE levels than the numbers -- for example a column
%holding both +0 and -0, which unique() collapses to one numeric value while
%the text forms can differ -- and would then report "levels would be merged"
%about valid data, which is the opposite of what happened.
if numel(unique(converted)) >= numel(unique(coldata(:)))
    return;
end

error('labels:numericlabelcollision',... %Error code and associated error
    ['Label column ''%s'' has numeric values that convert to the SAME text ' ...
    'label, so two or more distinct levels would be merged into one. ' ...
    'MATLAB renders numbers with about five significant digits by default, ' ...
    'so values differing only beyond that point become indistinguishable ' ...
    'labels. How to fix: supply the label column as text yourself (a cellstr ' ...
    'or string column), rounding or naming the levels the way you want them ' ...
    'reported, rather than relying on numeric-to-text conversion.'],colname);

end

function local_reject_missing(coldata,colname)
%Error if a label column carries missing values.
%
%char and cell are skipped: ismissing() is not meaningful for a char matrix
%(it compares element-wise and never flags the blank padding that char
%matrices carry), and a cell that is not cellstr has no text form at all, so
%it is caught by the conversion below with a clearer message.
if ischar(coldata) || iscell(coldata)
    return;
end

missmask = ismissing(coldata);

if any(missmask(:))
    badrows = find(missmask(:));
    error('labels:missingvalues',... %Error code and associated error
        strcat('Error: The ''',colname,''' column contains missing values\n',...
        'Missing values were found in ',[' ' num2str(numel(badrows))],...
        ' of ',[' ' num2str(numel(missmask))],' rows',...
        ' (first at row',[' ' num2str(badrows(1))],')\n\n',...
        'Label columns (id, group, event, time) identify levels of a facet,',...
        ' so a missing\nlabel cannot be filled in safely: every missing entry',...
        ' would collapse into a\nsingle level, pooling observations that',...
        ' belong to different participants,\ngroups, events or occasions.\n\n',...
        'How to fix: supply a label for every row, or drop the affected rows',...
        ' before\nloading them.\n'));
end

end

function out = local_to_cellstr(coldata,colname)
%Convert one label column to cellstr, preserving sort order and matching the
%text form psyrat_loadfile already produces for numeric columns.

if isdatetime(coldata)

    %ISO 8601, so that sorting the TEXT gives the same order as sorting the
    %dates. unique() sorts, and MATLAB's default datetime text form is
    %day-first ('20-Jan-2024'), which sorts by day-of-month: 20-Jan-2024 and
    %15-Jun-2024 come back in the OPPOSITE order from the dates themselves
    %(measured). Occasion levels are built by unique(), so the default form
    %would reorder the occasions of a test-retest design.
    %
    %BE ACCURATE ABOUT THE STAKES. This protects LABELS, not coefficients.
    %PsyRAT treats occasion as a fully exchangeable random facet:
    %psyrat_rel_trt takes occasions only as a scalar count nocc, REL.time is
    %consumed only through numel(), and the Stan occasion effect is i.i.d., so
    %permuting occasion levels does not move a variance component. What a
    %wrong sort corrupts is which session a reader sees as "Occasion 1" in the
    %output tables and figures. Worth getting right; not a numerical bug.
    out = cellstr(coldata(:),local_datetime_format(coldata));

elseif ischar(coldata)

    %A char matrix is one label per ROW, so it must NOT be linearized first:
    %cellstr(['ab';'cd'](:)) returns {'a';'c';'b';'d'}, four labels instead of
    %two. cellstr() on the matrix itself does the right thing and strips the
    %blank padding.
    out = cellstr(coldata);

else

    %string, categorical, duration, and cell arrays holding string scalars.
    try
        out = cellstr(coldata(:));
    catch ME
        error('labels:unsupportedtype',... %Error code and associated error
            strcat('Error: The ''',colname,''' column is of type ''',...
            class(coldata),''', which cannot be converted to text\n\n',...
            'Label columns (id, group, event, time) must hold text or',...
            ' numbers.\n',...
            'How to fix: convert the column before loading, for example with',...
            '\nstring() or cellstr().\n\n',...
            'The underlying error was:',[' ' ME.message],'\n'));
    end

end

out = out(:);

end

function fmt = local_datetime_format(coldata)
%Pick the shortest ISO 8601 form that does not merge distinct instants.
%Date-only occasions read far better as '2024-01-20' than as
%'2024-01-20 00:00:00', but dropping the time would collapse two sessions
%recorded on the same day into one occasion, so the time is kept whenever any
%value carries one.
tod = timeofday(coldata(~ismissing(coldata)));

if isempty(tod) || all(tod == 0)
    fmt = 'uuuu-MM-dd';
else
    fmt = 'uuuu-MM-dd HH:mm:ss';
end

end
