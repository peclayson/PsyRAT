function report = psyrat_preflight_validate(varargin)
%Validate PsyRAT inputs before running expensive estimation workflows.
%
% report = psyrat_preflight_validate('datatable',tbl,'sserrvar',1,'diffest',2)
% report = psyrat_preflight_validate('savename','out.psyrat','savepath','/tmp', ...
%     'requireSaveName',true,'requireSavePath',true)
%
% Returns a struct with fields:
%   ok       - true when no blocking validation errors were found
%   errors   - struct array with id, what, howto, and message
%   warnings - reserved for non-fatal diagnostics
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

if mod(length(varargin),2)
    error('varargin:incomplete', ...
        ['Input arguments are incomplete. ' ...
        'Provide name/value pairs, for example ' ...
        '''psyrat_preflight_validate(''datatable'',tbl,''sserrvar'',1)''.']);
end

opts = localDefaults();

for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};

    if ~(ischar(key) || (isstring(key) && isscalar(key)))
        error('varargin:keytype', ...
            'Parameter names must be character vectors or scalar strings.');
    end

    switch lower(char(key))
        case 'datatable'
            opts.datatable = val;
        case 'idheader'
            opts.idheader = char(string(val));
        case 'measheader'
            opts.measheader = char(string(val));
        case 'eventheader'
            opts.eventheader = char(string(val));
        case 'groupheader'
            opts.groupheader = char(string(val));
        case 'sserrvar'
            opts.sserrvar = val;
        case 'diffest'
            opts.diffest = val;
        case 'diffwpcov'
            opts.diffwpcov = val;
        case 'diffrescor'
            opts.diffrescor = val;
        case 'dynrel'
            opts.dynrel = val;
        case 'splits'
            opts.splits = val;
        case 'family'
            % Observation likelihood family: 'gaussian' (default) or 'gamma'
            % (scaled chi-square, log-linked mean, strictly positive support).
            opts.family = lower(char(string(val)));
        case 'savename'
            opts.savename = char(string(val));
        case 'savepath'
            opts.savepath = char(string(val));
        case 'requiresavename'
            opts.requireSaveName = logical(val);
        case 'requiresavepath'
            opts.requireSavePath = logical(val);
        case 'checkdiffpairs'
            % Legacy option accepted for compatibility. Standard
            % difference-score preflight no longer requires paired counts.
        otherwise
            error('varargin:unknown', ...
                'Unknown validation option ''%s''. Check help for valid options.', ...
                char(string(key)));
    end
end

report = localEmptyReport();

% Save-file validations.
if opts.requireSaveName && isempty(opts.savename)
    report = localAddError(report, 'psyrat_data:savename', ...
        'Output filename was not provided.', ...
        'Set psyrat_data.proc.savename to a non-empty filename such as ''study.psyrat''.');
end

if ~isempty(opts.savename) && any(isspace(opts.savename))
    report = localAddError(report, 'psyrat_data:savename', ...
        sprintf('Output filename contains whitespace: ''%s''.', opts.savename), ...
        'Rename the output file to remove spaces (CmdStan does not support whitespace in filenames).');
end

if opts.requireSavePath && isempty(opts.savepath)
    report = localAddError(report, 'psyrat_data:savepath', ...
        'Output save path was not provided.', ...
        'Set psyrat_data.proc.savepath to an existing directory path.');
end

if ~isempty(opts.savepath)
    if ~exist(opts.savepath, 'dir')
        report = localAddError(report, 'psyrat_data:savepath', ...
            sprintf('Output save path does not exist: ''%s''.', opts.savepath), ...
            'Create the directory first or choose an existing folder for psyrat_data.proc.savepath.');
    elseif any(isspace(opts.savepath))
        report = localAddError(report, 'psyrat_data:savepath', ...
            sprintf('Output save path contains whitespace: ''%s''.', opts.savepath), ...
            'Use a save path with no spaces (CmdStan does not support whitespace in paths).');
    end
end

% Table/data validations.
%
%istable is tested BEFORE isempty, deliberately. A table carrying column
%headers but no data rows has size [0 nvar], so prod(size) == 0 and isempty()
%is TRUE for it -- indistinguishable, under the old "if ~isempty(...)" gate,
%from the opts.datatable = [] default that means "no table was supplied". That
%collapsed the two cases and skipped EVERY table validation below for an empty
%dataset. Measured on a 0x3 table before this guard existed: report.ok = 1 with
%zero errors and zero warnings for the plain case, for diffest = 2, for
%diffest = 3, and even when the column headers were entirely wrong -- the
%missing-header check never ran either, so nothing stopped an empty dataset
%from going on to estimation. psyrat_loadfile has no row-count guard of its
%own, so a headers-only file reaches here intact.
%
%One error and no column-level checks is the intended shape: with no rows the
%FILE is the problem, and advice about which column to assign cannot be acted
%on until that is fixed.
if istable(opts.datatable)
    if height(opts.datatable) == 0
        report = localAddError(report, 'preflight:emptyTable', ...
            'The dataset is a table with column headers but no data rows.', ...
            'Check that the file was exported with its single-trial rows, then reload it.');
    else
        report = localValidateTable(report, opts);
    end
elseif ~isempty(opts.datatable)
    report = localAddError(report, 'varargin:datatabletype', ...
        'The provided datatable input is not a MATLAB table.', ...
        'Convert your dataset to a table before running estimation.');
end

report.ok = isempty(report.errors);

end

function opts = localDefaults()
opts = struct();
opts.datatable = [];
opts.idheader = 'id';
opts.measheader = 'meas';
opts.eventheader = 'event';
opts.groupheader = 'group';
opts.sserrvar = [];
opts.diffest = [];
opts.diffwpcov = [];
opts.diffrescor = [];
opts.dynrel = [];
opts.splits = [];
opts.family = 'gaussian';
opts.savename = '';
opts.savepath = '';
opts.requireSaveName = false;
opts.requireSavePath = false;
end

function report = localValidateTable(report, opts)
tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasID, ~] = localFindColumn(colnames, opts.idheader);
[hasMeas, meascol] = localFindColumn(colnames, opts.measheader);

if ~hasID || ~hasMeas
    missing = {};
    if ~hasID
        missing{end+1} = opts.idheader; %#ok<AGROW>
    end
    if ~hasMeas
        missing{end+1} = opts.measheader; %#ok<AGROW>
    end

    report = localAddError(report, 'varargin:colheaders', ...
        sprintf('Dataset is missing required column header(s): %s.', strjoin(missing, ', ')), ...
        sprintf('Ensure your table includes both ''%s'' and ''%s'' columns, then rerun.', ...
        opts.idheader, opts.measheader));
end

if hasMeas && ~isnumeric(tbl.(meascol))
    report = localAddError(report, 'meas:notnumeric', ...
        sprintf('Measurement column ''%s'' is not numeric.', meascol), ...
        'Select a numeric measurement column (ERP amplitudes/latencies) and reload the dataset.');
end

% Family-appropriateness: the Gamma / scaled-chi-square family has strictly
% positive support (mean-linked, right-skewed data such as single-trial
% time-frequency power). Non-positive scores are outside its support and would
% make the log link / likelihood undefined, so block them before estimation.
% The default Gaussian family is unaffected (ERP amplitudes are routinely
% negative), so this guard is gamma-only.
if hasMeas && isnumeric(tbl.(meascol)) && strcmpi(opts.family, 'gamma')
    meas = tbl.(meascol);
    if any(meas(:) <= 0)
        report = localAddError(report, 'meas:notpositive', ...
            sprintf(['Measurement column ''%s'' contains non-positive values, ', ...
            'which the Gamma / scaled-chi-square family does not support.'], meascol), ...
            ['Use the Gamma family only with strictly positive scores (for ', ...
            'example single-trial time-frequency power). Shift, transform, or ', ...
            'switch to the Gaussian family for data that include zero or ', ...
            'negative values.']);
    end
end

% Statistics and Machine Learning Toolbox for the gamma concurrent difference
% designs (analyses 7, 8 and 10 with residual covariance on). Their
% Gaussian-copula extractors call normcdf and gaminv only AFTER sampling, so a
% MATLAB without the toolbox paid for a full CmdStan fit and then lost it to an
% undefined-function error (audit 2026-09-05). Refuse here, before anything is
% compiled, on both entry routes. The license feature string is the one
% psyrat_estimation_engines_available already tests; exist() covers a licensed
% but uninstalled toolbox.
if strcmpi(opts.family, 'gamma') && localDiffCovRequested(opts) && ...
        ~(license('test', 'Statistics_Toolbox') == 1 && exist('gaminv', 'file') > 0)
    report = localAddError(report, 'preflight:statisticsToolboxRequired', ...
        ['The gamma family with residual covariance (a concurrent difference design) ', ...
        'reads its results through a Gaussian copula that needs the Statistics and ', ...
        'Machine Learning Toolbox (normcdf, gaminv), which this MATLAB does not have.'], ...
        ['Install the Statistics and Machine Learning Toolbox, set residual covariance ', ...
        'to No (a non-concurrent difference), or use the Gaussian family.']);
end

% B27 (owner ruling 2026-08-17: reject, do not convert). A LOGICAL label
% column bypasses psyrat_normalize_label_columns -- its own comment
% deliberately declines to convert logical, because string(true) is "true"
% where the legacy numeric idiom spells '1', and inventing a third label
% spelling for a type no shipped path converts was outside that fix's
% scope. Unconverted, the downstream behavior is unpredictable and
% case-dependent: some group chunkers key labels through char() and fail on
% the control characters char(0)/char(1), one route samples to completion
% before dying at its label append, the gamma string()-keyed cases complete
% while INVENTING 'false'/'true' labels the user never typed, a logical id
% can silently no-op the test-retest participant-wide exclusion, and a
% logical time column can even run correctly by exchangeability. Rejecting
% the TYPE here replaces that spectrum with one clear, actionable message
% on every entry route, including the direct-table one.
%
% Resolution is by EXACT name over the normalizer's own column set
% (ismember, mirroring psyrat_normalize_label_columns and the engine's
% case-sensitive datatable.group dot access), deliberately NOT this file's
% localFindColumn: its first case-insensitive match would let a
% case-variant decoy (a non-logical 'Group' alongside a logical 'group')
% hide the column the engine actually reads, and would reject tables whose
% exact-named column is healthy.
labelcols = {'id','group','event','time'};
for lc = 1:numel(labelcols)
    thiscol = labelcols{lc};
    if ismember(thiscol, colnames) && islogical(tbl.(thiscol))
        report = localAddError(report, 'labels:logicaltype', ...
            sprintf(['Column ''%s'' is a logical (true/false) array, which ', ...
            'cannot be used as a label column.'], thiscol), ...
            sprintf(['Convert ''%s'' to numeric codes, a categorical, or text ', ...
            'labels (for example double(t.%s) or categorical(t.%s)), then rerun.'], ...
            thiscol, thiscol, thiscol));
    end
end

hasMissing = false;
try
    hasMissing = any(any(ismissing(tbl)));
catch
    % If ismissing is unavailable for a specific class, keep best-effort behavior.
    hasMissing = false;
end
if hasMissing
    report = localAddError(report, 'varargin:missingcells', ...
        'Input data include missing cells, which CmdStan cannot process.', ...
        'Remove or impute missing rows so every row has id, measurement, and facet values before estimation.');
end

% Data-integrity guards (audit 2026-09-05). Each is a condition the loader and
% the engine accepted without comment and that changed the estimand or the
% exported numbers silently: a non-finite score (ismissing does not flag Inf);
% a label carrying a CSV delimiter, quote, or line break (the hand-rolled CSV
% writers do not quote data cells, so the exported row misaligns); a
% participant listed under two groups (the engine numbers persons within
% group, so one person becomes two); and an occasion column with one level
% (the column's presence routes the run to the two-facet test-retest design,
% whose occasion components cannot be identified from one occasion). Two more
% are warnings because they can be legitimate: labels differing only in case
% or surrounding spaces, and exact duplicate rows.
report = localValidateDataIntegrity(report, tbl, colnames, hasMeas, meascol);

[hasTime, ~] = localFindColumn(colnames, 'time');

% Subject-level dynamic/conditional reliability (Rast & Clayson) of a trial +
% occasion two-facet two-event difference score (dynrel=2, diffest=2, sserrvar=2,
% analysis 16) IS supported with a time/occasion facet.
isDynrelSubjDiffTrt = isequal(opts.dynrel, 2) && isequal(opts.diffest, 2);

% Unconditional per-subject test-retest reliability (trt_sserrvar, analysis 25)
% is ALSO supported with a time/occasion facet: subject-level error variance in
% the crossed test-retest design when non-difference (diffest=1) and non-dynamic
% (dynrel~=2). This is the direct test-retest analog of the single-occasion
% subject-level model (analysis 6). Difference designs (diffest=2/3) and the
% dynamic non-difference subject-level combination stay blocked below / by the
% engine.
isTrtSserr = ~isequal(opts.diffest, 2) && ~isequal(opts.diffest, 3) && ...
    ~isequal(opts.dynrel, 2);

% Subject-level dynamic/conditional reliability (Rast & Clayson) of a NON-difference
% trial + occasion two-facet design (dynrel=2, diffest~=2/3, sserrvar=2, analysis 27)
% IS supported with a time/occasion facet: each participant's conditional reliability
% phi_s(z) is reported from the per-subject residual the location-scale model already
% estimates (the dynamic analog of trt_sserrvar, analysis 25).
isDynrelSubjNondiff = isequal(opts.dynrel, 2) && ~isequal(opts.diffest, 2) && ...
    ~isequal(opts.diffest, 3);

% Person-specific dynamic NONCONCURRENT difference-of-differences (dynrel=2,
% diffest=3, sserrvar=2; analysis 28 single-occasion, 29 with a time/occasion
% facet) IS supported: the engine routes it to the four-margin designs (both
% families since 2026-08-20) and enforces its own parameterization guards.
isDynrelSubjDoD = isequal(opts.dynrel, 2) && isequal(opts.diffest, 3);

if hasTime && ~isempty(opts.sserrvar) && isequal(opts.sserrvar, 2) && ...
        ~isDynrelSubjDiffTrt && ~isTrtSserr && ~isDynrelSubjNondiff && ...
        ~isDynrelSubjDoD
    report = localAddError(report, 'preflight:timeWithSubjectReliabilityUnsupported', ...
        ['Single-subject error variance estimation is not supported when a time/occasion facet is present ' ...
        '(test-retest mode).'], ...
        'Set sserrvar/ssubjrel to 1 (off) or remove the time facet for this run.');
end

% Two-event difference-score reliability (diffest=2) WITH a time/occasion facet
% is the supported two-facet (test-retest) difference pipeline (analysis 10), so
% it is permitted. Subject-level difference scores (sserrvar=2) with time are
% supported only for the dynamic/conditional reliability model (dynrel=2,
% analysis 16). Still block the genuinely unsupported combinations with time:
% difference-of-differences (diffest=3) and non-dynamic subject-level difference
% scores (diffest=2 with sserrvar=2 and dynrel ~= 2).
blockDiffWithTime = false;
if ~isempty(opts.diffest)
    if isequal(opts.diffest, 3)
        %the person-specific dynamic DoD (dynrel=2 with sserrvar=2, analysis
        %29) supports a time/occasion facet; every other DoD-with-time
        %combination stays blocked
        if ~(isequal(opts.dynrel, 2) && isequal(opts.sserrvar, 2))
            blockDiffWithTime = true;
        end
    elseif isequal(opts.diffest, 2) && ~isempty(opts.sserrvar) && ...
            isequal(opts.sserrvar, 2) && ~isequal(opts.dynrel, 2)
        blockDiffWithTime = true;
    end
end
if hasTime && blockDiffWithTime
    report = localAddError(report, 'preflight:timeWithDiffUnsupported', ...
        ['This difference-score configuration is not supported when a ', ...
        'time/occasion facet is present. Two-facet (test-retest) two-event ', ...
        'difference scores are supported; difference-of-differences (diffest=3) ', ...
        'and subject-level difference scores (sserrvar=2) with time are not.'], ...
        'Disable difference-score estimation (diffest=1), run on a single-occasion dataset, or use two-event difference scores without subject-level error variance.');
end

if ~isempty(opts.diffest) && isequal(opts.diffest, 3) && ...
        ~isempty(opts.sserrvar) && isequal(opts.sserrvar, 2) && ...
        ~isequal(opts.dynrel, 2)
    %dynrel=2 is exempt: the person-specific dynamic DoD (analyses 28/29) is
    %exactly the subject-level combination, and the engine enforces its own
    %gamma location-scale requirements
    report = localAddError(report, 'preflight:dodWithSubjectReliabilityUnsupported', ...
        ['Difference-of-differences estimation with subject-level error variance ' ...
        '(diffest=3 with sserrvar=2) is not supported in the current estimation workflow ' ...
        'without dynamic/conditional reliability (dynrel=2).'], ...
        'Use two-event difference-score mode (diffest=2) for subject-level reliability, disable subject-level error variance, or enable dynrel=2 for the person-specific dynamic difference-of-differences (supported in both the Gaussian and Gamma location-scale families).');
end

if ~isempty(opts.diffest) && isequal(opts.diffest, 2)
    report = localValidateDiffMode(report, opts);
elseif ~isempty(opts.diffest) && isequal(opts.diffest, 3)
    report = localValidateDoDiffMode(report, opts);
end

% Reliability-with-splits validations (splits=2 parallel, 3 nonparallel).
if ~isempty(opts.splits) && (isequal(opts.splits, 2) || isequal(opts.splits, 3))
    report = localValidateSplits(report, opts);
end

% Non-fatal warnings for permitted-but-shaky configurations. These never
% set report.ok to false; they flag designs where variance components may be
% weakly identified so the user can decide whether to proceed (API-008). Each
% helper guards on the columns it needs and is a no-op when they are absent.
report = localWarnFewParticipants(report, opts);
report = localWarnLowMedianTrials(report, opts);
report = localWarnIncompleteEventTimeCells(report, opts);
end

function report = localValidateDiffMode(report, opts)
tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasEvent, eventcol] = localFindColumn(colnames, opts.eventheader);
if ~hasEvent
    report = localAddError(report, 'preflight:diffNeedsEventColumn', ...
        ['Difference-score reliability requires an event/condition column, ' ...
        sprintf('but ''%s'' was not found.', opts.eventheader)], ...
        'Specify an event column in the processing setup before enabling diffest=2.');
    return;
end

events = unique(localToStringColumn(tbl.(eventcol)), 'stable');
nevents = numel(events);

if nevents < 2
    report = localAddError(report, 'preflight:diffNeedsTwoEvents', ...
        sprintf('Difference-score reliability found %d event type(s); exactly 2 are required.', nevents), ...
        'Select two event conditions (X and Y) for difference-score estimation.');
    return;
elseif nevents > 2
    report = localAddError(report, 'preflight:diffNeedsTwoEvents', ...
        sprintf('Difference-score reliability found %d event types; exactly 2 are required.', nevents), ...
        'Limit selected events to two conditions before running diffest=2.');
    return;
end

% Standard G-theory difference-score estimation uses long-form observations
% and permits unequal event counts. Paired-count checks belong only to
% residual-covariance paths that explicitly reshape events into paired rows.
if localDiffCovRequested(opts)
    [hasID, idcol] = localFindColumn(colnames, opts.idheader);
    [hasGroup, groupcol] = localFindColumn(colnames, opts.groupheader);
    if hasID
        [hasPairs, pairWhat] = localHasDiffPairs(tbl, idcol, eventcol, hasGroup, groupcol, events);
        if ~hasPairs
            report = localAddError(report, 'preflight:diffCovNoPairs', pairWhat, ...
                ['Residual covariance estimation for subject-level difference scores needs at least one ', ...
                'within-participant paired observation for the selected events in each modeled group.']);
        end
    end
end
end

function report = localValidateDoDiffMode(report, opts)
tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasEvent, eventcol] = localFindColumn(colnames, opts.eventheader);
if ~hasEvent
    report = localAddError(report, 'preflight:dodNeedsEventColumn', ...
        ['Difference-of-differences reliability requires an event column, ' ...
        sprintf('but ''%s'' was not found.', opts.eventheader)], ...
        'Specify an event column in the processing setup before enabling diffest=3.');
    return;
end

events = unique(localToStringColumn(tbl.(eventcol)), 'stable');
nevents = numel(events);
if nevents ~= 4
    report = localAddError(report, 'preflight:dodNeedsFourEvents', ...
        sprintf('Difference-of-differences reliability found %d event type(s); exactly 4 are required.', nevents), ...
        'Limit selected events to exactly 4 conditions before running diffest=3.');
end

%PER-PERSON coverage, for the person-specific dynamic DoD only (analyses
%28/29): each participant's own per-cell trial counts are that design's
%D-study divisors, so a participant with ZERO trials in any cell makes their
%conditional coefficient undefined. Without this check the fit completes and
%every later summary/viewer attempt dies in the composite kernel's n > 0
%guard with no participant named - a completed multi-hour run rendered
%unviewable (pre-merge review, 2026-08-13). Fail loudly BEFORE fitting, the
%same doctrine as the analogous per-person check in the difference path.
%Scoped to dynrel=2 + sserrvar=2 so the group-level static DoD (analysis 9),
%whose kernel uses group-level counts, is byte-unchanged.
if nevents == 4 && isequal(opts.dynrel, 2) && isequal(opts.sserrvar, 2)
    [hasId, idcol] = localFindColumn(colnames, 'id');
    if hasId
        ids = localToStringColumn(tbl.(idcol));
        evs = localToStringColumn(tbl.(eventcol));
        uid = unique(ids, 'stable');
        bad = {};
        for k = 1:numel(uid)
            have = unique(evs(strcmp(ids, uid{k})));
            if numel(have) < 4
                bad{end+1} = char(uid{k}); %#ok<AGROW>
            end
        end
        if ~isempty(bad)
            shown = strjoin(bad(1:min(5, numel(bad))), ', ');
            if numel(bad) > 5
                shown = sprintf('%s, ... (%d total)', shown, numel(bad));
            end
            report = localAddError(report, 'preflight:dodPersonMissingCell', ...
                sprintf(['Person-specific dynamic difference-of-differences ' ...
                'requires every participant to have trials in all four ' ...
                'events; missing for participant(s): %s.'], shown), ...
                ['Exclude participants without all four conditions, or ' ...
                'choose events every participant completed.']);
        end
    end
end
end

function report = localValidateSplits(report, opts)
% Validate reliability-with-splits inputs. The loader renames the user's
% items-per-split (n_i) column to 'weight' (the number of items/trials
% averaged into each split mean). Splits are a Stage-1 design (non-difference,
% group-level, no dynrel); the engine and the GUI both enforce that
% combination, so only the n_i column itself is validated here.
tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasWeight, weightcol] = localFindColumn(colnames, 'weight');
if ~hasWeight
    % Required for both parallel and nonparallel splits: the parallel path
    % checks n_i equality, the nonparallel path uses 1/n_i as a precision
    % weight. The loader also requires it, but flag it here so the GUI/CLI
    % preflight gives the same up-front diagnostic.
    report = localAddError(report, 'preflight:splitsNeedsWeight', ...
        'Reliability with data splits requires an items-per-split (n_i) column, but none was found.', ...
        'Select the split-weight (n_i) column in the processing setup before running a splits analysis.');
    return;
end

% The loader has already validated that weight is a positive whole number;
% drop any non-finite entries defensively before measuring variation.
w = double(tbl.(weightcol));
w = w(:);
w = w(isfinite(w));
if isempty(w)
    return;
end

allEqual = (max(w) - min(w)) == 0;

if isequal(opts.splits, 3)
    % Nonparallel splits identify the per-item residual (sig_err) and the
    % person-by-split interaction (sig_posxid/sig_splitxid) ONLY through
    % differences in the 1/n_i precision weights. Equal n_i leaves them
    % confounded, so hard-block that degenerate case. When n_i varies but
    % only a little, the components may still be weakly identified
    % (empirically n_i in 26-30, CV ~ 0.05, collapsed them while
    % n_i in [8,100], CV ~ 0.7, recovered) -- warn but allow the run.
    if allEqual
        report = localAddError(report, 'preflight:splitsNoVariation', ...
            ['Nonparallel splits were requested, but every split has the same ' ...
            'items-per-split (n_i) count.'], ...
            ['Provide unequal n_i across splits (the nonparallel observed-design ' ...
            'model is identified only by n_i variation), or declare the splits parallel.']);
    else
        cv = std(w) / mean(w);
        if cv < localSplitsLowCV()
            report = localAddWarning(report, 'preflight:splitsLowVariation', ...
                sprintf(['Items-per-split (n_i) variation is low (coefficient of ' ...
                'variation %.2f). The per-item residual and the person-by-split ' ...
                'interaction may be weakly identified.'], cv), ...
                ['The run will still proceed; interpret the nonparallel-splits ' ...
                'estimates with caution, or use splits with more variable n_i.']);
        end
    end
elseif isequal(opts.splits, 2)
    % Parallel splits assume equal n_i: the split mean is the score and
    % n' = number of splits, reusing the standard one-facet/test-retest
    % machinery. Unequal n_i suggests the nonparallel model fits better;
    % warn but do not block (the user's declaration is honored).
    if ~allEqual
        report = localAddWarning(report, 'preflight:splitsUnequalWeight', ...
            ['Parallel splits were declared, but the items-per-split (n_i) ' ...
            'counts are unequal across splits.'], ...
            ['Parallel splits assume equal n_i; consider declaring the splits ' ...
            'nonparallel so the unequal weights are modeled.']);
    end
end
end

function thr = localSplitsLowCV()
% Heuristic threshold for flagging low items-per-split (n_i) variation under
% nonparallel splits (splits=3). Below this coefficient of variation, the
% per-item residual and the person-by-split interaction are weakly
% identified. Anchored to the empirical evidence: n_i in 26-30 (CV ~ 0.05)
% collapsed the components, while n_i in [8,100] (CV ~ 0.7) recovered
% cleanly; 0.25 sits between them, flagging the collapse region while
% staying quiet for well-spread designs. Heuristic -- tune if needed.
thr = 0.25;
end

function report = localWarnFewParticipants(report, opts)
% Warn when a modeled group has few participants. Variance-component
% estimates (especially between-person and interaction terms) are unstable
% and prior-dominated at small N. Estimation fits a separate model per group,
% so the count is evaluated per group when a group column is present.
minParticipants = 20;  % below this, flag the design as underpowered

tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasID, idcol] = localFindColumn(colnames, opts.idheader);
if ~hasID
    return;
end
idvals = localToStringColumn(tbl.(idcol));

[hasGroup, groupcol] = localFindColumn(colnames, opts.groupheader);
if hasGroup
    groupvals = localToStringColumn(tbl.(groupcol));
    groups = unique(groupvals, 'stable');
else
    groupvals = strings(size(idvals));
    groups = "";
end

for g = 1:numel(groups)
    if hasGroup
        groupmask = groupvals == groups(g);
    else
        groupmask = true(size(idvals));
    end
    npart = numel(unique(idvals(groupmask)));
    if npart < minParticipants
        if hasGroup
            what = sprintf(['Group ''%s'' has only %d participant(s), ' ...
                'fewer than the recommended minimum of %d.'], ...
                char(groups(g)), npart, minParticipants);
        else
            what = sprintf(['The dataset has only %d participant(s), ' ...
                'fewer than the recommended minimum of %d.'], ...
                npart, minParticipants);
        end
        report = localAddWarning(report, 'preflight:fewParticipants', what, ...
            ['Variance-component estimates can be unstable and ' ...
            'prior-dominated with few participants. The run will still ' ...
            'proceed; interpret the reliability estimates with added ' ...
            'caution or collect more participants.']);
    end
end
end

function report = localWarnLowMedianTrials(report, opts)
% Warn when the median number of trials per design cell is low. Thin cells
% leave within-person measurement error poorly estimated. A cell is the
% crossing of the id facet with whichever facet columns are present
% (event, time, group).

% Splits runs are exempt: one row is a split mean (a chunk of trials), not a
% single trial, so the per-cell row count is the splits-per-person count, not
% the trials-per-cell quantity this check is about. The real per-split item
% count (n_i) and its required variation are handled by localValidateSplits.
if ~isempty(opts.splits) && (isequal(opts.splits, 2) || isequal(opts.splits, 3))
    return;
end

minMedianTrials = 10;  % below this median, flag the cells as thin

tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

[hasID, idcol] = localFindColumn(colnames, opts.idheader);
[hasMeas, meascol] = localFindColumn(colnames, opts.measheader);
if ~hasID || ~hasMeas
    return;
end

% Build the grouping (cell) key from id plus any present facet columns.
% 'time' is detected by its literal header to match the rest of preflight.
groupVars = {idcol};
candidateFacets = {opts.eventheader, 'time', opts.groupheader};
for f = 1:numel(candidateFacets)
    [hasFacet, facetcol] = localFindColumn(colnames, candidateFacets{f});
    if hasFacet
        groupVars{end+1} = facetcol; %#ok<AGROW>
    end
end

% Count observations (trials) per cell. The varfun(@length, ...) tally
% mirrors the per-id trial count used in psyrat_relsummary.
trltable = varfun(@length, tbl, 'GroupingVariables', groupVars, ...
    'InputVariables', meascol);
medianTrials = median(trltable.GroupCount);

if medianTrials < minMedianTrials
    facetNames = groupVars(2:end);
    if isempty(facetNames)
        cellDesc = 'each participant';
    else
        cellDesc = sprintf('id x %s', strjoin(facetNames, ' x '));
    end
    what = sprintf(['Median trials per cell is %.1f, fewer than the ' ...
        'recommended minimum of %d (cell = %s).'], ...
        medianTrials, minMedianTrials, cellDesc);
    report = localAddWarning(report, 'preflight:lowMedianTrials', what, ...
        ['Within-person measurement error is poorly estimated when cells ' ...
        'are thin. The run will still proceed; consider whether the trial ' ...
        'counts are adequate for the scored component.']);
end
end

function report = localWarnIncompleteEventTimeCells(report, opts)
% Warn when the observed event-by-occasion grid has empty (id, event, time)
% cells. The loader's coverage guards (meas:mismatchedevents,
% meas:mismatchedoccasions) enforce UNIFORMITY of per-participant counts --
% distinct-event counts and event-by-occasion cell counts respectively --
% not COMPLETENESS of the grid, so a dataset in which every participant
% misses the same NUMBER of cells (e.g. routine artifact rejection emptying
% one condition in one participant-occasion) loads cleanly and estimates
% (B24's reachability fixture, 2026-08-16). Estimation is valid with empty
% cells: Bayesian variance-component estimation does not require a complete
% design, and quantities for empty cells borrow strength from the model.
% The user should still be TOLD, because an empty cell has no observed-score
% summaries and its estimates are model-based. Warning only; never flips
% report.ok (owner ruling 2026-08-16: warn, do not block).

tbl = opts.datatable;
colnames = tbl.Properties.VariableNames;

% The grid is defined by event x occasion, so both facet columns must be
% present ('time' by its literal header, matching the rest of preflight).
[hasID, idcol] = localFindColumn(colnames, opts.idheader);
[hasEvent, eventcol] = localFindColumn(colnames, opts.eventheader);
[hasTime, timecol] = localFindColumn(colnames, 'time');
if ~hasID || ~hasEvent || ~hasTime
    return;
end

ids = localToStringColumn(tbl.(idcol));
events = localToStringColumn(tbl.(eventcol));
times = localToStringColumn(tbl.(timecol));

uid = unique(ids);
uev = unique(events);
utm = unique(times);

% Present cells as unique ROWS of the three label columns -- membership is
% tested row-wise (ismember over a table), never through a joined string
% key: a delimiter-based key can COLLIDE when a label itself contains the
% delimiter, which both fabricates missing cells in a cardinality count and
% hides genuinely missing ones in a key lookup (caught by the 2026-08-16
% refuter pass; the first cut of this check had exactly that defect).
present = unique(table(ids, events, times, ...
    'VariableNames', {'id','event','time'}));
nexpected = numel(uid) * numel(uev) * numel(utm);

% Enumerate the empty cells: the single source for the count, the affected
% participants, and the examples, so the message cannot disagree with
% itself. Preflight tables are small; the triple loop is fine.
missingdesc = {};
affected = false(numel(uid), 1);
for i = 1:numel(uid)
    for e = 1:numel(uev)
        for t = 1:numel(utm)
            combo = table(uid(i), uev(e), utm(t), ...
                'VariableNames', {'id','event','time'});
            if ~ismember(combo, present)
                affected(i) = true;
                missingdesc{end+1} = sprintf('%s: %s @ %s', ...
                    char(uid(i)), char(uev(e)), char(utm(t))); %#ok<AGROW>
            end
        end
    end
end
nmissing = numel(missingdesc);
if nmissing <= 0
    return;
end

nshow = min(3, numel(missingdesc));
examples = strjoin(missingdesc(1:nshow), '; ');
if numel(missingdesc) > nshow
    examples = sprintf('%s; ...', examples);
end

% "observations", not "trials": splits tables (one row per split mean)
% reach this warning too, and present/absent semantics hold there as well.
what = sprintf(['The observed event x occasion grid is incomplete: ' ...
    '%d of %d (participant x event x occasion) cells contain no ' ...
    'observations, affecting %d of %d participant(s) (e.g. %s).'], ...
    nmissing, nexpected, sum(affected), numel(uid), examples);
report = localAddWarning(report, 'preflight:incompleteEventTimeCells', what, ...
    ['Estimation proceeds: variance-component estimation permits ' ...
    'incomplete designs, and quantities for empty cells are model-based ' ...
    'rather than observed. Verify the empty cells reflect intended ' ...
    'exclusions (e.g. artifact rejection) rather than a data-preparation ' ...
    'error.']);
end

function report = localValidateDataIntegrity(report, tbl, colnames, hasMeas, meascol)
%Data-integrity guards; see the call site for the rationale. Every check is
%guarded on the columns it needs and is a no-op when they are absent, and the
%label columns are resolved by EXACT canonical name (the loader has already
%normalized them), mirroring the B27 logical-type check above.
labelcols = {'id','group','event','time'};

% 1. Non-finite scores (NaN is already reported as a missing cell above).
if hasMeas && isnumeric(tbl.(meascol))
    meas = tbl.(meascol);
    isinf_row = ~isfinite(meas(:)) & ~isnan(meas(:));
    if any(isinf_row)
        report = localAddError(report, 'meas:notfinite', ...
            sprintf(['Measurement column ''%s'' contains %d non-finite (Inf or -Inf) ', ...
            'value(s); the first is in row %d.'], meascol, nnz(isinf_row), ...
            find(isinf_row, 1, 'first')), ...
            'Remove or correct the non-finite scores before estimation; every score must be a finite number.');
    end
end

% 2. Reserved characters in labels (error) and case/whitespace variants (warning).
for lc = 1:numel(labelcols)
    thiscol = labelcols{lc};
    if ~ismember(thiscol, colnames) || islogical(tbl.(thiscol))
        continue;
    end
    labels = localToStringColumn(tbl.(thiscol));
    labels = labels(~ismissing(labels));
    if isempty(labels)
        continue;
    end

    bad = contains(labels, ',') | contains(labels, '"') | ...
        contains(labels, newline) | contains(labels, char(13));
    if any(bad)
        report = localAddError(report, 'labels:reservedcharacters', ...
            sprintf(['Column ''%s'' has %d label(s) containing a comma, a double quote, ', ...
            'or a line break (for example ''%s'').'], thiscol, nnz(bad), ...
            char(labels(find(bad, 1, 'first')))), ...
            ['Rename those labels without commas, quotes, or line breaks: the exported ', ...
            'CSV tables write labels unquoted, so such a label would shift every later ', ...
            'column of its row.']);
    end

    u = unique(labels);
    folded = lower(strtrim(u));
    [~, ~, ic] = unique(folded);
    counts = accumarray(ic(:), 1);
    k = find(counts > 1, 1, 'first');
    if ~isempty(k)
        pair = u(ic == k);
        report = localAddWarning(report, 'preflight:labelCaseVariants', ...
            sprintf(['Column ''%s'' has labels that differ only in letter case or ', ...
            'surrounding spaces (for example ''%s'' and ''%s''); PsyRAT treats them ', ...
            'as distinct levels.'], thiscol, char(pair(1)), char(pair(2))), ...
            ['If they name the same level, make the spelling identical before ', ...
            'estimation; otherwise ignore this warning.']);
    end
end

% 3. A participant listed under more than one group (WARNING). Filed as a refusal by
%    the 2026-09-05 audit and downgraded by owner ruling on 2026-09-07 when the merged
%    unit suite showed it rejecting the toolbox's own dev fixture (testdata.csv labels a
%    within-person task as group): ids numbered within each group are a legitimate
%    convention, and the core keys persons by (group, id), so the within-group estimates
%    are unaffected either way. The warning names the count and an example id and says
%    when the condition belongs in the event column instead; the GUI shows it modally
%    (G8) and the scripted route prints it, and the run proceeds.
if ismember('id', colnames) && ismember('group', colnames) && ...
        ~islogical(tbl.id) && ~islogical(tbl.group)
    ids = localToStringColumn(tbl.id);
    groups = localToStringColumn(tbl.group);
    [uid, ~, ic] = unique(ids);
    ngroups = zeros(numel(uid), 1);
    for k = 1:numel(uid)
        ngroups(k) = numel(unique(groups(ic == k)));
    end
    multi = uid(ngroups > 1);
    if ~isempty(multi)
        report = localAddWarning(report, 'preflight:groupWithinParticipant', ...
            sprintf(['%d participant id(s) appear under more than one group (for example ', ...
            '''%s''). PsyRAT treats group as a between-person factor and counts such a ', ...
            'participant once per group, as separate people.'], ...
            numel(multi), char(multi(1))), ...
            ['If these rows belong to the same person, use the group column only for ', ...
            'between-person grouping and put the within-person condition in the event ', ...
            'column. If participant ids are numbered within each group, no action is needed.']);
    end
end

% 4. An occasion column with a single level (error).
if ismember('time', colnames) && ~islogical(tbl.time)
    times = localToStringColumn(tbl.time);
    times = times(~ismissing(times));
    if ~isempty(times) && numel(unique(times)) < 2
        report = localAddError(report, 'preflight:singleOccasionLevel', ...
            sprintf(['The occasion column ''time'' has only one level (''%s''). With an ', ...
            'occasion column present PsyRAT fits the two-facet test-retest design, whose ', ...
            'occasion variance components cannot be identified from a single occasion.'], ...
            char(times(1))), ...
            ['Leave the Occasion ID column on none (omit timecol in psyrat_run) to fit ', ...
            'the single-occasion design, or supply data from at least two occasions.']);
    end
end

% 5. Exact duplicate rows (warning: identical scores can recur legitimately).
ndup = 0;
try
    ndup = height(tbl) - height(unique(tbl));
catch
    % a column type unique() cannot order: skip the check rather than fail preflight
    ndup = 0;
end
if ndup > 0
    report = localAddWarning(report, 'preflight:duplicateRows', ...
        sprintf('%d row(s) are exact duplicates of another row (same id, labels, and score).', ndup), ...
        ['If the duplicates are copies rather than genuinely repeated scores, remove ', ...
        'them: duplicated trials inflate every affected participant''s trial count and ', ...
        'the reliability projected at the observed counts.']);
end
end

function [tf, actualName] = localFindColumn(colnames, requested)
match = strcmpi(colnames, requested);
tf = any(match);
actualName = '';
if tf
    actualName = colnames{find(match, 1, 'first')};
end
end

function tf = localDiffCovRequested(opts)
tf = false;
if isfield(opts, 'diffwpcov') && ~isempty(opts.diffwpcov) && isequal(opts.diffwpcov, 2)
    tf = true;
end
if isfield(opts, 'diffrescor') && ~isempty(opts.diffrescor) && isequal(opts.diffrescor, 2)
    tf = true;
end
end

function [hasPairs, what] = localHasDiffPairs(tbl, idcol, eventcol, hasGroup, groupcol, events)
hasPairs = true;
what = '';

eventvals = localToStringColumn(tbl.(eventcol));
idvals = localToStringColumn(tbl.(idcol));

if hasGroup
    groupvals = localToStringColumn(tbl.(groupcol));
    groups = unique(groupvals, 'stable');
else
    groupvals = strings(size(idvals));
    groups = "";
end

for g = 1:numel(groups)
    if hasGroup
        groupmask = groupvals == groups(g);
    else
        groupmask = true(size(idvals));
    end

    ids = unique(idvals(groupmask), 'stable');
    npairs = 0;
    for i = 1:numel(ids)
        idmask = groupmask & (idvals == ids(i));
        n1 = sum(idmask & (eventvals == events(1)));
        n2 = sum(idmask & (eventvals == events(2)));
        npairs = npairs + min(n1, n2);
    end

    if npairs < 1
        hasPairs = false;
        if hasGroup
            what = sprintf(['No pairable observations were found for residual covariance ', ...
                'estimation in group=''%s''.'], char(groups(g)));
        else
            what = 'No pairable observations were found for residual covariance estimation.';
        end
        return;
    end
end
end

function out = localToStringColumn(v)
if isstring(v)
    out = v(:);
elseif iscellstr(v)
    out = string(v(:));
elseif ischar(v)
    out = string(cellstr(v));
elseif isnumeric(v) || islogical(v) || iscategorical(v)
    out = string(v(:));
else
    out = string(v(:));
end
end

function report = localEmptyReport()
report = struct();
report.ok = true;
report.errors = struct('id', {}, 'what', {}, 'howto', {}, 'message', {});
report.warnings = struct('id', {}, 'what', {}, 'howto', {}, 'message', {});
end

function report = localAddError(report, id, what, howto)
err = struct();
err.id = id;
err.what = what;
err.howto = howto;
err.message = sprintf('%s\nHow to fix: %s', what, howto);
report.errors(end+1) = err;
end

function report = localAddWarning(report, id, what, howto)
warn = struct();
warn.id = id;
warn.what = what;
warn.howto = howto;
warn.message = sprintf('%s\nRecommendation: %s', what, howto);
report.warnings(end+1) = warn;
end
