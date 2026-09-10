function report = psyrat_report(psyrat_data, varargin)
%Assemble a headless, render-free report of a reliability analysis.
%
%report = psyrat_report(psyrat_data)
%report = psyrat_report(psyrat_data, 'outdir', dir, 'format', '.xlsx')
%
%psyrat_report is the scripted/batch counterpart of the GUI results viewer. It
%does NOT recompute: the variance-component estimation (psyrat_computevarcomp)
%and the D-study summary (psyrat_relsummary) are the compute layer, and this
%function only assembles their results into return tables and, optionally,
%writes them to disk -- with no figures and no save dialogs. A typical headless
%workflow is therefore two clean steps:
%
%   [psyrat_data, ~] = psyrat_relsummary('psyrat_data', psyrat_data, ...);
%   report = psyrat_report(psyrat_data, 'outdir', '/path/to/out');
%
%Each exported table carries a provenance header (seed, engine, priors,
%convergence; see psyrat_provenance_lines) so a saved table can be interpreted
%on its own.
%
%Input
% psyrat_data - toolbox data structure AFTER psyrat_relsummary has populated
%   psyrat_data.relsummary. psyrat_data.rel.analysis selects which tables apply.
%
%Name/Value options
% 'outdir' - directory to write the tables into (default '': return only, no
%   files written).
% 'format' - '.xlsx' (default) or '.csv'; the export file format.
% 'CI'     - credible-interval width in (0,1). Governs every interval THIS
%   function computes from the posterior draws: the data-splits and
%   dynamic-reliability coefficients, the one-facet / test-retest
%   reliability-vs-n curves (.curve, and .diffcurve for trt_diff), and (G52)
%   both difference-of-differences tables (.dod_observed/.dod_varcomp). It
%   does NOT reach the .cutoff/.overall/.variance/.sscoeffs tables, whose
%   intervals were fixed at relsummary time (psyrat_data.relsummary.ciperc).
%   DEFAULT (G41): the width the summary was built at
%   (psyrat_data.relsummary.ciperc) when one is recorded, so a plain call
%   writes one width into every table of an export; .95 otherwise. The shared
%   provenance header records that stored width ("Credible interval: NN%");
%   the summary-free families (DoD, splits, dynrel) instead get a
%   "Credible interval" header line recording the width this report applied
%   (G52). An explicit 'CI' still wins, and when it differs from a stored
%   width the export header gains a note saying which intervals carry which
%   width (on the sserr family, whose report computes no CI-governed table,
%   the note says the option did not apply). An explicit EMPTY ('CI',[]) is
%   treated as unset and takes the default -- before G41 an empty propagated
%   into the quantile calls and hard-errored downstream. An explicit NaN
%   errors in the coefficient kernels on the summary-free families rather
%   than being disclosed as a mismatch (there is no stored width to
%   mismatch against).
% 'marginal' - dynamic-reliability only: 0 (default) typical-person surface
%   (delta_p = 0), 1 population-average (lognormal). Forced to 0 (with a warning
%   from psyrat_dynrel_summary) for the difference variants.
% 'reltype'  - dynamic-reliability two-facet variants only: 1/2/3 coefficient of
%   equivalence / stability / both-random (default 3). Inert otherwise.
% 'nocc'     - dynamic-reliability two-facet variants only: occasion n' ([]
%   default -> 1; or 'observed'). Inert otherwise.
% 'ngrid'    - dynamic-reliability only: grid points per dimension axis for the
%   reliability surface (default 25).
% 'ntrials'  - one-facet / test-retest families only: the maximum trial count for
%   the reliability-vs-n curve (default 50, matching the GUI plot default). The
%   curve is tabulated over n = 1..ntrials.
%
%Output
% report - struct with fields:
%   .analysis    - the resolved figure-layer analysis type ('sing'/'sing_diff'/
%                  'trt'/'trt_diff'/'sserr'/'dod'/'splits'/'dynrel').
%   .rel_analysis- the underlying psyrat_data.rel.analysis value.
%   .tables      - struct of returned MATLAB tables (e.g. .cutoff, .overall, and
%                  .curve [reliability vs n trials] for the one-facet/test-retest
%                  families; .diffcurve [difference-score reliability vs n trials]
%                  additionally for the two-facet difference design 'trt_diff';
%                  .variance/.sscoeffs for the subject-level (sserr)
%                  family; .coefficients/.varcomp for the data-splits family;
%                  .surface/.varcomp[/.ssrel] for the dynamic-reliability family).
%   .provenance  - the provenance header lines (cell of char).
%   .gbias_note  - (data-splits only) the generalizability lower-bound caveat,
%                  also written into the header of the exported splits tables.
%   .estimand_note - (dynamic-reliability only) the surface estimand (typical
%                  person vs population average; plus reltype/occasion n' for the
%                  two-facet variants), also written into each exported header.
%   .curve_note  - (one-facet / test-retest only) the reliability-vs-n curve
%                  coefficient label, also written into each exported header. That
%                  curve is PER-CONDITION on every analysis, the difference designs
%                  included, and the note says so.
%   .diffcurve_note - ('trt_diff' only) the difference-score curve's coefficient
%                  label plus its estimand caveat: the curve is drawn at equal
%                  trial counts in both conditions, while the overall table's
%                  difference score uses the observed per-condition counts, so the
%                  two coincide only when those counts are equal. Also written
%                  into each exported header.
%   .sscoeff_note- (subject-level sserr only) the per-subject coefficient label
%                  (dependability vs generalizability), also written into each
%                  exported header.
%   .cutoff_note - (one-facet / test-retest families only) explains the -1
%                  sentinel cells the .cutoff table can carry: a cutoff that
%                  could not be calculated (no trial count in the projected
%                  range reached the threshold) or one that extrapolates beyond
%                  every participant's observed count. Empty for a clean run;
%                  also written into each exported header. The GUI raises the
%                  same two conditions as dialogs (G1).
%   .strata      - (dynamic-reliability only) the raw psyrat_dynrel_summary
%                  per-stratum struct (the unflattened G/D/ICC surface) for
%                  programmatic use.
%   .files       - cellstr of written file paths ({} when 'outdir' is unset).
%
%Coverage note: this supports the families whose results are render-free: the
%one-facet / test-retest table builders ('gui',0) (sing, sing_diff, trt,
%trt_diff; cutoff + overall tables plus the reliability-vs-n curve tabulated from
%the calc layer -- psyrat_rel_sing / psyrat_rel_trt -- the headless analog of the
%GUI trials plot), the subject-level error-variance family (ic_sserrvar,
%trt_sserrvar, ic_diff_sserrvar; the group-level SD/SEM/ICC table via the
%render-free psyrat_variancet plus the per-subject coefficients flattened into
%one long table), difference-of-differences (ic_dodiff; per-group observed +
%variance-component tables), the nonparallel data-splits family (ic_splits,
%trt_splits; per-stratum coefficients + variance components via
%psyrat_splits_summary), and the dynamic-reliability family (ic_dynrel and its
%difference / two-facet / subject-level variants; flattened G(z)/D(z)/ICC
%surface + variance components [+ per-participant table for subject-level
%variants] via psyrat_dynrel_summary).

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

% ---- options ---------------------------------------------------------------
outdir = local_opt(varargin,'outdir','');
fmt    = local_opt(varargin,'format','.xlsx');
%G41: default the interval width to the one the relsummary was built at, so a
%plain psyrat_report call writes ONE width into every table of an export (the
%TestCiCoverageConsistency one-width rule). An explicit 'CI' still wins; when
%it differs from the stored width the header gains a disclosure note below.
%.95 remains the fallback for structs that carry no built summary (the
%splits/dynrel routes can arrive summary-free).
ciperc = local_opt(varargin,'CI',[]);
ciperc_stored = [];
if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'ciperc') && ...
        isnumeric(psyrat_data.relsummary.ciperc) && ...
        isscalar(psyrat_data.relsummary.ciperc) && ~isnan(psyrat_data.relsummary.ciperc)
    ciperc_stored = psyrat_data.relsummary.ciperc;
end
if isempty(ciperc)
    if ~isempty(ciperc_stored)
        ciperc = ciperc_stored;
    else
        ciperc = .95;
    end
end
% Dynamic-reliability passthrough options (inert for the other families).
marginal = local_opt(varargin,'marginal',0);
reltype  = local_opt(varargin,'reltype',3);
nocc     = local_opt(varargin,'nocc',[]);
ngrid    = local_opt(varargin,'ngrid',25);
% One-facet / test-retest reliability-vs-n curve upper bound (GUI plot default).
ntrials  = local_opt(varargin,'ntrials',50);
if ~isempty(fmt) && fmt(1) ~= '.'
    fmt = ['.' fmt];
end

% ---- set up the report container -------------------------------------------
report = struct();
report.analysis = '';
report.rel_analysis = local_rel_analysis(psyrat_data);
report.tables = struct();
report.provenance = psyrat_provenance_lines(psyrat_data);
report.files = {};

%G41 disclosure flag: an explicit 'CI' differing from the width the summary
%was built at still puts two widths in one export (the stored-width tables
%cannot be recomputed here), so the header must say which intervals carry
%which width. The note TEXT is composed after the table dispatch below,
%because what it should honestly say depends on whether this family's
%report computes any CI-governed table at all (the sserr route computes
%none, so the original wording asserted intervals that did not exist --
%adversarial wave, 2026-08-29). NaN-safe on purpose: the gate is spelled
%~(<= tol) rather than (> tol) so an explicit 'CI',NaN counts as a
%mismatch and gets disclosed instead of silently defeating the comparison.
report.ci_note = '';
ci_mismatch = ~isempty(ciperc_stored) && ...
    ~(abs(ciperc - ciperc_stored) <= 1e-12);

% Carry forward any admissibility note the summary built. psyrat_relsummary
% attaches it, so a caller that ran the summary first (which this function
% requires anyway) gets it for free. Empty for a clean run.
if isfield(psyrat_data,'relsummary') && ...
        isfield(psyrat_data.relsummary,'admissibility_note')
    report.admissibility_note = psyrat_data.relsummary.admissibility_note;
end

% ---- dispatch by analysis family (render-free; 'gui',0) --------------------
switch psyrat_analysis_family(report.rel_analysis)
    case 'dod'
        % Difference-of-differences: the DoD table builders read rel directly
        % (no relsummary) and are render-free with 'gui',0.
        report.analysis = 'dod';
        [report.tables.dod_observed, report.tables.dod_varcomp, dod_adm_note] = ...
            local_dod_tables(psyrat_data, ciperc);
        % The DoD kernel's admissibility record has no relsummary to travel
        % through, so the carry-forward above cannot pick it up and a headless
        % DoD export had no way to learn a variance came out negative (RC-43).
        % Appended rather than assigned: a run that also produced a summary note
        % keeps it, and both reach the same header block below.
        if ~isempty(dod_adm_note)
            if isfield(report,'admissibility_note') && ...
                    ~isempty(report.admissibility_note)
                report.admissibility_note = sprintf('%s\n%s', ...
                    report.admissibility_note, dod_adm_note);
            else
                report.admissibility_note = dod_adm_note;
            end
        end
    case 'splits'
        % Nonparallel data-splits: psyrat_splits_summary reads REL.out directly
        % (no relsummary, like DoD) and is render-free. Assemble the same two
        % tables the splits viewer shows -- reliability coefficients and variance
        % components -- and surface the generalizability lower-bound caveat.
        report.analysis = 'splits';
        [report.tables.coefficients, report.tables.varcomp, report.gbias_note] = ...
            local_splits_tables(psyrat_data, ciperc);
    case 'dynrel'
        % Dynamic reliability (Rast & Clayson): psyrat_dynrel_summary reads
        % REL.out directly (no relsummary, like splits/DoD) and is render-free.
        % The G(z)/D(z) surface the GUI only plots is flattened into a table so
        % the headless report carries the reliability values; the variance
        % components and (subject-level variants only) the per-participant table
        % mirror the dynrel viewer's exports. psyrat_analysis_family is the
        % shared classifier (same accepted set as psyrat_dynrel_summary).
        report.analysis = 'dynrel';
        [report.tables.surface, report.tables.varcomp, ssrelTbl, ...
            report.estimand_note, report.strata] = ...
            local_dynrel_tables(psyrat_data, ciperc, marginal, reltype, nocc, ngrid);
        if ~isempty(ssrelTbl)
            report.tables.ssrel = ssrelTbl;
        end
    otherwise
        % psyrat_analysis_family returned 'sing', 'trt', or '' (unrecognized).
        % A future dynrel/splits/dod variant not yet added to the classifier
        % also lands here (its family resolves to ''); relsummary is never
        % populated for those families, so guard with a clear "not yet covered"
        % error rather than letting it hit the relsummary check below and emit
        % the misleading "relsummary is not populated" message.
        famstr = lower(string(report.rel_analysis));
        if contains(famstr,"dynrel") || contains(famstr,"splits") || ...
                contains(famstr,"dodiff")
            error('psyrat_report:unsupported', ...
                ['Headless reporting is not yet available for analysis ''%s'' ', ...
                '(a recognized family variant not in the report dispatch list). ', ...
                'How to fix: use the GUI results viewer for this analysis, or add ', ...
                'this variant to the matching case in psyrat_report.'], ...
                report.rel_analysis);
        end
        % One-facet / test-retest families. These tables are derived from the
        % D-study summary, so relsummary must already be populated.
        if ~isstruct(psyrat_data) || ~isfield(psyrat_data,'relsummary') || ...
                isempty(psyrat_data.relsummary)
            error('psyrat_report:norelsummary', ...
                ['psyrat_data.relsummary is not populated. How to fix: run ', ...
                'psyrat_relsummary on the estimated data first (the GUI does ', ...
                'this when you open the results viewer), then call psyrat_report.']);
        end
        report.analysis = local_analysis_type(psyrat_data);
        switch report.analysis
            case {'sing','sing_diff','trt_diff'}
                report.tables.cutoff  = psyrat_depcutofft('psyrat_data',psyrat_data,'gui',0);
                report.tables.overall = psyrat_depoverallt('psyrat_data',psyrat_data,'gui',0);
                [report.tables.curve, report.curve_note] = ...
                    local_curve_table(psyrat_data, ciperc, ntrials);
                % RC-31. The two-facet difference design also exports a genuine
                % difference-score curve, the headless analog of
                % psyrat_trt_diffvtrialsplot. It is a SEPARATE table rather than
                % extra rows in the one above because it is stratified Group x n,
                % not Group x Event x n: a difference has no per-condition
                % stratification, so folding it in would need a pseudo-event and
                % would break the Event column for any consumer grouping on it.
                if strcmp(report.analysis,'trt_diff')
                    [report.tables.diffcurve, report.diffcurve_note] = ...
                        local_diffcurve_table(psyrat_data, ciperc, ntrials);
                end
            case 'trt'
                report.tables.cutoff  = psyrat_trt_relcutofft('psyrat_data',psyrat_data,'gui',0);
                report.tables.overall = psyrat_trt_reloverallt('psyrat_data',psyrat_data,'gui',0);
                [report.tables.curve, report.curve_note] = ...
                    local_curve_table(psyrat_data, ciperc, ntrials);
            case 'sserr'
                % Subject-level error-variance family. Two parts, matching the GUI
                % variance viewer: the group-level SD/SEM/ICC table (now render-free
                % via the C3 gui gate on psyrat_variancet) and the per-subject
                % coefficients flattened into one long table (a port of the
                % psyrat_savesscoeffs CSV loop). No reliability-vs-n curve here --
                % the GUI shows a per-subject plot for sserr, not a D-study
                % projection, and the per-subject table already carries those values.
                report.tables.variance = ...
                    psyrat_variancet('psyrat_data',psyrat_data,'gui',0);
                [report.tables.sscoeffs, report.sscoeff_note] = ...
                    local_sserr_table(psyrat_data);
            otherwise
                error('psyrat_report:unsupported', ...
                    ['Headless reporting is not yet available for analysis ''%s''. ', ...
                    'How to fix: use the GUI results viewer for this analysis, or ', ...
                    'wait for the later headless-report increment that covers it.'], ...
                    report.rel_analysis);
        end
        %G1 (2026-09-05): explain the D-study's -1 sentinel cells for a headless
        %reader. The GUI raises "Cutoff not calculable" and "Extrapolation
        %beyond data" dialogs from relerr (psyrat_relfigures,
        %psyrat_criterionfigures); this route has no dialog and relerr is not
        %stored on psyrat_data, so the per-cell markers are read instead. Empty
        %for a clean run, so ordinary exports are unchanged.
        if any(strcmp(report.analysis, {'sing','sing_diff','trt','trt_diff'}))
            report.cutoff_note = local_cutoff_note(psyrat_data);
        end
end

%G41/G52 disclosure text, composed now that report.tables says what this
%family actually computed. When CI-governed tables exist, name both widths:
%.curve/.diffcurve sit beside intervals fixed at summary time, whereas the
%two DoD tables (CI-governed since G52) are the WHOLE export, so their
%wording differs. The splits/dynrel coefficient tables never reach here
%with a stored width on any shipped path, so they cannot mismatch. When no
%CI-governed table exists (the sserr route: ciperc reaches nothing there),
%say the option did not apply -- silence would hide that the explicit 'CI'
%was ignored, and the original wording asserted report-computed intervals
%that did not exist.
if ci_mismatch
    if isfield(report.tables, 'dod_observed') || ...
            isfield(report.tables, 'dod_varcomp')
        %only reachable on a hand-built struct (no shipped path attaches a
        %relsummary to a DoD rel), but if one arrives, tell the truth: G52
        %threads ciperc into BOTH DoD tables, so nothing in this export
        %keeps the stored width
        report.ci_note = sprintf(['Note: every interval in this report ' ...
            'uses %g%%, the ''CI'' option, not the %g%% width recorded ' ...
            'above.'], 100 * ciperc, 100 * ciperc_stored);
    elseif isfield(report.tables, 'curve') || isfield(report.tables, 'diffcurve')
        report.ci_note = sprintf(['Note: the reliability-vs-n curve ' ...
            'intervals in this report use %g%%, the ''CI'' option; ' ...
            'intervals fixed at summary time use the %g%% width recorded ' ...
            'above.'], 100 * ciperc, 100 * ciperc_stored);
    else
        report.ci_note = sprintf(['Note: the ''CI'' option (%g%%) does ' ...
            'not affect this family''s tables; all intervals use the ' ...
            '%g%% width recorded above.'], 100 * ciperc, 100 * ciperc_stored);
    end
end

%G52 disclosure: families with no built summary (DoD always; the splits and
%dynrel routes arrive summary-free) take their interval width from ciperc,
%but psyrat_provenance_lines records a width only from relsummary.ciperc --
%so their exports could not say which width their intervals use, breaking
%the one-width self-description doctrine (TestCiCoverageConsistency) the
%moment the width became user-settable. Record the width this report
%actually applied whenever it computed CI-governed tables and no stored
%width exists for the provenance line to record.
if isempty(ciperc_stored) && ...
        (isfield(report.tables, 'dod_observed') || ...
        isfield(report.tables, 'dod_varcomp') || ...
        strcmp(report.analysis, 'splits') || ...
        strcmp(report.analysis, 'dynrel'))
    report.ci_width_note = sprintf('Credible interval: %g%%', 100 * ciperc);
end

% ---- write the tables, if requested ----------------------------------------
if ~isempty(outdir)
    if exist(outdir,'dir') ~= 7
        mkdir(outdir);
    end
    [base, base_is_default] = local_base(psyrat_data);
    % Analysis-specific caveats (when present) are appended to every exported
    % table header so a standalone CSV/XLSX records them: the data-splits
    % generalizability lower-bound caveat and the dynamic-reliability surface
    % estimand. Concatenated (not grown in a loop) to keep checkcode clean.
    extralines = {};
    %the resolved-width line leads: it plays the role the provenance block's
    %"Credible interval" line plays for summary-carrying families (G52)
    if isfield(report,'ci_width_note') && ~isempty(report.ci_width_note)
        extralines = [extralines; {report.ci_width_note}];
    end
    if isfield(report,'gbias_note') && ~isempty(report.gbias_note)
        extralines = [extralines; {report.gbias_note}];
    end
    if isfield(report,'estimand_note') && ~isempty(report.estimand_note)
        extralines = [extralines; {report.estimand_note}];
    end
    if isfield(report,'curve_note') && ~isempty(report.curve_note)
        extralines = [extralines; {report.curve_note}];
    end
    if isfield(report,'diffcurve_note') && ~isempty(report.diffcurve_note)
        extralines = [extralines; {report.diffcurve_note}];
    end
    if isfield(report,'ci_note') && ~isempty(report.ci_note)
        extralines = [extralines; {report.ci_note}];
    end
    if isfield(report,'sscoeff_note') && ~isempty(report.sscoeff_note)
        extralines = [extralines; {report.sscoeff_note}];
    end
    % G1: the -1 sentinel explanation travels in the same header block, for the
    % same reason as the admissibility note below: a headless export has no
    % dialog, and a -1 in a cutoff cell needs its meaning next to it.
    if isfield(report,'cutoff_note') && ~isempty(report.cutoff_note)
        extralines = [extralines; {report.cutoff_note}];
    end
    % The admissibility note goes in the SAME header block, which is what makes
    % it reach a headless user at all: a console warning does not survive into
    % an exported file, and there is no dialog on this path. Empty for a clean
    % run, so ordinary exports are unchanged.
    if isfield(report,'admissibility_note') && ~isempty(report.admissibility_note)
        extralines = [extralines; {report.admissibility_note}];
    end
    fn = fieldnames(report.tables);
    % Resolve every target path first. When base fell back to the generic
    % 'psyrat_report' stem (proc.savename unset, the headless/batch path),
    % refuse to overwrite pre-existing files: a second run into the same outdir
    % would otherwise silently clobber the first with no warning. A run with an
    % explicit savename keeps overwriting on purpose (deliberate re-runs).
    fpaths = cell(numel(fn),1);
    for k = 1:numel(fn)
        fpaths{k} = fullfile(outdir, sprintf('%s_%s%s', base, fn{k}, fmt));
    end
    if base_is_default
        % isfile (not exist(...,'file')==2) so the check is strictly the literal
        % target path: exist resolves bare names against the MATLAB path, which
        % could spuriously flag a same-named file elsewhere when outdir is
        % relative and block a legitimate first write.
        clash = cellfun(@isfile, fpaths);
        if any(clash)
            hit = fpaths(clash);
            error('psyrat_report:overwrite', ...
                ['Refusing to overwrite an existing report file (%s). The tables ', ...
                'are named generically (''%s_*'') because psyrat_data.proc.savename ', ...
                'is not set, so a second run into the same outdir would silently ', ...
                'clobber the first. How to fix: set proc.savename, or write each ', ...
                'run to its own (empty) outdir.'], hit{1}, base);
        end
    end
    for k = 1:numel(fn)
        head = local_table_header(psyrat_data, fn{k}, report.provenance, extralines);
        psyrat_write_table(fpaths{k}, head, report.tables.(fn{k}));
        report.files{end+1} = fpaths{k};
    end
end
end

%-------------------------------------------------------------------------
% Local helpers
%-------------------------------------------------------------------------
function v = local_opt(args,name,default)
%Return the value following the first occurrence of option NAME in a name/value
%varargin. Only odd (name) positions are scanned, so a VALUE that happens to
%equal an option keyword (e.g. an outdir literally named 'CI') cannot be
%mistaken for the flag and shift the value lookup by one.
v = default;
for i = 1:2:numel(args)-1
    a = args{i};
    if (ischar(a) || (isstring(a) && isscalar(a))) && strcmpi(name,a)
        v = args{i+1};
        return;
    end
end
end

function a = local_rel_analysis(psyrat_data)
if isfield(psyrat_data,'rel') && isstruct(psyrat_data.rel) && ...
        isfield(psyrat_data.rel,'analysis') && ~isempty(psyrat_data.rel.analysis)
    a = psyrat_data.rel.analysis;
else
    a = 'ic';
end
end

function analysis = local_analysis_type(psyrat_data)
%Map psyrat_data.rel.analysis to the figure-layer analysis string, mirroring
%the viewer routing (psyrat_startview_sing). Only the render-free families are
%resolved here; everything else falls through to '' and triggers the
%unsupported error in the caller.
a = local_rel_analysis(psyrat_data);
switch a
    case {'ic',''}
        analysis = 'sing';
    case 'ic_diff'
        analysis = 'sing_diff';
    case 'trt'
        analysis = 'trt';
    case 'trt_diff'
        analysis = 'trt_diff';
    case {'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}
        analysis = 'sserr'; % subject-level error-variance family (C3)
    otherwise
        analysis = a; % unsupported in this increment; caller errors with it
end
end

function note = local_cutoff_note(psyrat_data)
%Explain the -1 sentinel cells a .cutoff table can carry (G1, 2026-09-05).
%
%psyrat_relsummary records two D-study failures per group x event cell and
%raises neither as an error. When no trial count in the projected range
%reaches the threshold it stores trlcutoff = -1 and writes -1 into the
%coefficient slots ("Cutoff not calculable"); when the cutoff it found exceeds
%every participant's observed count it keeps that cutoff and overwrites the
%coefficient slots with -1 ("Extrapolation beyond data"). The one-facet
%families keep the coefficient under event(eloc).rel, the test-retest families
%under event(eloc).relcutoff. The GUI turns both conditions into dialogs from
%relerr, psyrat_relsummary's second output; a headless caller sees no dialog,
%and relerr is not stored on psyrat_data, so this function reads the per-cell
%markers instead, which is the contract psyrat_relsummary's own header states
%("No cutoff was found is recorded by relerr.trlcutoff and by the per-cell
%event(eloc).trlcutoff = -1"). A coefficient of exactly -1 cannot come out of
%the ratio formulas, so the value is unambiguous as a marker. The cells are
%left as the GUI export prints them (both routes write the same -1); this note
%only explains them. The difference-score row of the one-facet and two-facet
%difference designs carries its own search result in diffscore.trlcutoff.ntrials,
%-1 when no count up to the search bound reached the threshold
%(psyrat_diffrel_trlcutoff). Returns '' for a clean run. The wording mirrors the
%two GUI dialogs, with the coefficient named as the summary named it.
note = '';
if ~isfield(psyrat_data,'relsummary') || ~isstruct(psyrat_data.relsummary) || ...
        ~isfield(psyrat_data.relsummary,'group')
    return;
end
relsummary = psyrat_data.relsummary;
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits') && ...
        ~isempty(psyrat_data.rel.splits)
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);
relname = 'dependability';
if isfield(relsummary,'gcoeff_name') && strcmpi(relsummary.gcoeff_name,'gen')
    relname = 'generalizability';
end
notcalc = {};
extrap = {};
for gloc = 1:numel(relsummary.group)
    grp = relsummary.group(gloc);
    if isfield(grp,'event') && isstruct(grp.event)
        for eloc = 1:numel(grp.event)
            ev = grp.event(eloc);
            label = local_cutoff_cell_label(grp, ev);
            hascut = isfield(ev,'trlcutoff') && isnumeric(ev.trlcutoff) && ...
                isscalar(ev.trlcutoff);
            coef = [];
            if isfield(ev,'relcutoff') && isstruct(ev.relcutoff) && isfield(ev.relcutoff,'m')
                coef = ev.relcutoff.m;
            elseif isfield(ev,'rel') && isstruct(ev.rel) && isfield(ev.rel,'m')
                coef = ev.rel.m;
            end
            hascoef = isnumeric(coef) && isscalar(coef);
            if hascut && ev.trlcutoff == -1
                notcalc{end+1} = label; %#ok<AGROW>
            elseif hascut && hascoef && coef == -1 && ev.trlcutoff > 0
                extrap{end+1} = sprintf('%s (cutoff %d)', label, ev.trlcutoff); %#ok<AGROW>
            end
        end
    end
    %the difference-score row (one-facet and two-facet difference designs)
    if isfield(grp,'diffscore') && isstruct(grp.diffscore) && ...
            isfield(grp.diffscore,'trlcutoff') && isstruct(grp.diffscore.trlcutoff) && ...
            isfield(grp.diffscore.trlcutoff,'ntrials')
        dn = grp.diffscore.trlcutoff.ntrials;
        if isnumeric(dn) && isscalar(dn) && dn == -1
            notcalc{end+1} = local_cutoff_cell_label(grp, struct('name','diff score')); %#ok<AGROW>
        end
    end
end
parts = {};
if ~isempty(notcalc)
    parts{end+1} = sprintf(['Note: %s cutoffs for adequate %s could not be ' ...
        'calculated for %s (data are too variable or there are not enough %s); ' ...
        'the affected %s and coefficient cells print -1.'], ...
        L.Unit, relname, strjoin(notcalc, ', '), L.units, lower(L.Cutoff));
end
if ~isempty(extrap)
    parts{end+1} = sprintf(['Note: not enough %s are present in the current data ' ...
        'for %s; the cutoff is an extrapolation beyond the data and the ' ...
        'coefficient cells print -1.'], L.units, strjoin(extrap, ', '));
end
if ~isempty(parts)
    note = strjoin(parts, newline);
end
end

function label = local_cutoff_cell_label(grp, ev)
%Name a relsummary cell the way the cutoff tables label their rows: the group
%name, the event name, or 'group - event'; 'the measurement' when the run has
%neither (psyrat_relsummary names that lone event 'measure').
gname = '';
ename = '';
if isfield(grp,'name') && ischar(grp.name)
    gname = strtrim(grp.name);
end
if isfield(ev,'name') && ischar(ev.name) && ~strcmpi(ev.name,'measure')
    ename = strtrim(ev.name);
end
if isempty(gname) && isempty(ename)
    label = 'the measurement';
elseif isempty(ename)
    label = gname;
elseif isempty(gname)
    label = ename;
else
    label = [gname ' - ' ename];
end
end

function [base, isdefault] = local_base(psyrat_data)
%Filename stem for written tables: the saved .psyrat base when available, else
%the generic 'psyrat_report' default. isdefault flags the fallback so the writer
%can refuse to clobber existing generically-named files (see the write block).
base = 'psyrat_report';
isdefault = true;
if isfield(psyrat_data,'proc') && isstruct(psyrat_data.proc) && ...
        isfield(psyrat_data.proc,'savename') && ~isempty(psyrat_data.proc.savename)
    [~,base] = fileparts(psyrat_data.proc.savename);
    isdefault = false;
end
end

function [obsTbl, vcTbl, adm_note] = local_dod_tables(psyrat_data, ciperc)
%Build the difference-of-differences observed and variance-component tables for
%every group, render-free ('gui',0). The DoD builders require the contrast as
%numeric event indices, so the stored rel.dod_map (event NAMES) is converted to
%indices into rel.events. A 'Group' column is prepended so the per-group rows
%are identifiable in the combined table.
%G52: ciperc is threaded into both builders. They carry their own .95
%defaults, so before this an explicit 'CI' on psyrat_report's DoD route was
%silently ignored -- and with no relsummary there is no stored width, so
%neither the provenance line nor the mismatch note could disclose the drop.
%The honored width is now recorded in the export header by the caller's
%ci_width_note (the summary-free analog of the provenance block's
%"Credible interval" line), keeping the export self-describing.
rel = psyrat_data.rel;
if ~isfield(rel,'dod_map') || isempty(rel.dod_map)
    error('psyrat_report:dodmap', ...
        ['DoD reporting requires rel.dod_map (the four contrast events). ', ...
        'How to fix: re-run the difference-of-differences analysis, which ', ...
        'records the contrast map.']);
end
events = cellstr(string(rel.events(:)));
mapnames = cellstr(string(rel.dod_map(:)));
[tf, map] = ismember(mapnames, events);
if numel(map) ~= 4 || ~all(tf)
    error('psyrat_report:dodmap', ...
        'rel.dod_map must name four events present in rel.events.');
end
map = map(:)';

if isfield(rel,'groups') && ~isempty(rel.groups)
    groupnames = cellstr(string(rel.groups(:)));
else
    groupnames = {'All'};
end

obsParts = cell(numel(groupnames),1);
vcParts  = cell(numel(groupnames),1);
admParts = cell(numel(groupnames),1);
for g = 1:numel(groupnames)
    ot = psyrat_dod_observedt('psyrat_data',psyrat_data,'groupidx',g,'map',map, ...
        'CI',ciperc,'gui',0);
    % psyrat_dod_varcompt defaults savefile=1 (its export callback), so pass
    % savefile=0 explicitly to get a pure render-free, dialog-free return.
    % Its second output is the admissibility record for this group (RC-43).
    % Only this one is taken, not psyrat_dod_observedt's: called without 'obs'
    % both builders resolve the same observed trial counts, so the two records
    % describe ONE evaluation and adding both would double the counts.
    [vt, admParts{g}] = psyrat_dod_varcompt('psyrat_data',psyrat_data, ...
        'groupidx',g,'map',map,'CI',ciperc, ...
        'gui',0,'savefile',0);
    obsParts{g} = addvars(ot, repmat(groupnames(g),height(ot),1), ...
        'Before', 1, 'NewVariableNames', 'Group');
    vcParts{g}  = addvars(vt, repmat(groupnames(g),height(vt),1), ...
        'Before', 1, 'NewVariableNames', 'Group');
end
obsTbl = vertcat(obsParts{:});
vcTbl  = vertcat(vcParts{:});

% One note for the whole report, with the counts summed across groups -- the
% aggregation psyrat_admissibility_collect performs for every other analysis
% family. Records from different groups are different evaluations, so these DO
% add. Empty for a clean run.
adm_note = psyrat_admissibility_note([admParts{:}]);
end

function [coefTbl, vcTbl, gbias_note] = local_splits_tables(psyrat_data, ciperc)
%Build the nonparallel data-splits coefficient and variance-component tables,
%render-free, at the observed design (no D-study projection -- split means lose
%the item partition). psyrat_splits_summary is the compute layer; this only
%flattens its per-stratum summary into tables, mirroring the splits viewer's
%coefDisplay/concatVarcomp assembly (psyrat_startview_splits) so the headless
%output matches the GUI export.
rel = psyrat_data.rel;
summ = psyrat_splits_summary(rel, 'CI', ciperc, 'obs', [], 'nocc', []);
gbias_note = summ.gbias_note;

% Unit label ('split' vs 'trial') for the n-per-person column, matching the
% viewer. Real splits results carry rel.splits (= 2/3); default to nonparallel
% here because this branch only runs for the data-splits analyses.
if isfield(rel,'splits') && ~isempty(rel.splits)
    L = psyrat_unitlabels(rel.splits);
else
    L = psyrat_unitlabels(3);
end

% ---- coefficients: one row per stratum x coefficient (shared builder) -------
% psyrat_splits_coef_rows is the single source of truth, also used by the splits
% viewer's coefDisplay. cell2table with makeValidName matches the viewer's
% saveCoef export exactly.
[data, cols] = psyrat_splits_coef_rows(summ, L);
coefTbl = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));

% ---- variance components: stack per-stratum tables (shared builder) ----------
% psyrat_splits_varcomp_concat is the single source of truth, also used by the
% splits viewer's concatVarcomp (stratum label + observed design columns).
vcTbl = psyrat_splits_varcomp_concat(summ);
end

function [surfaceTbl, vcTbl, ssrelTbl, estimand_note, strata] = ...
    local_dynrel_tables(psyrat_data, ciperc, marginal, reltype, nocc, ngrid)
%Build the dynamic-reliability tables, render-free, from psyrat_dynrel_summary
%(the compute layer, which reads REL.out directly -- no relsummary / D-study
%projection). Three outputs: the flattened G(z)/D(z)/ICC reliability surface (the
%headless analog of the GUI plot, which the viewer never tabulates), the stacked
%variance components (port of the dynrel viewer's concatVarcomp), and, for
%subject-level variants only, the per-participant table (port of concatSsrel; []
%otherwise). marginal/reltype/nocc/ngrid are passed straight through; they are
%inert for the variants that do not consume them (psyrat_dynrel_summary's
%contract). The varcomp/ssrel ports match what the viewer's save callbacks
%export (the raw concatenated tables), not the on-screen display subset.
rel = psyrat_data.rel;
summ = psyrat_dynrel_summary(rel,'CI',ciperc,'marginal',marginal,...
    'reltype',reltype,'nocc',nocc,'ngrid',ngrid);
strata = summ.strata;
ndim = summ.ndim;

% ---- reliability surface: one row per z grid point (port of the plot) -------
% G and D are structs with ll/pt/ul over the grid; ICCg/ICCd point estimates.
% 1D -> ngrid rows; 2D -> ngrid^2 rows with z1 outer, z2 inner, matching the
% surface arrays (rows index z1, columns index z2 -- psyrat_rel_dynrel). Two-facet
% variants (non-empty stratum nocc) also record the occasion n' as a column: the
% G(z)/D(z) values depend on it, so without it a saved surface table could not be
% interpreted from its own columns (the estimand note alone carries it otherwise).
% The gamma difference variant computes its surface at PER-EVENT n' (st.obs_ev),
% which for an unbalanced design (the reference ERP design is 49 error vs 341
% correct trials) is the only honest annotation: the scalar st.obs is their mean
% and names a count neither event has. Emit the two exact counts INSTEAD of the
% scalar there, so the saved table stays interpretable from its own columns --
% the same rationale as the occasion column above.
hasocc = ~isempty(strata) && ~isempty(strata(1).nocc);
hasev  = ~isempty(strata) && isfield(strata,'obs_ev') && ...
    numel(strata(1).obs_ev) == 2;
%the DoD variants carry FOUR per-cell counts; itemize all four for the same
%reason the two-event pair is itemized
hasev4 = ~isempty(strata) && isfield(strata,'obs_ev') && ...
    numel(strata(1).obs_ev) == 4;
cols = {'Stratum','z1'};
if ndim == 2; cols = [cols, {'z2'}]; end
if hasev4
    cols = [cols, {'n trials cell 1','n trials cell 2',...
        'n trials cell 3','n trials cell 4'}];
elseif hasev
    cols = [cols, {'n trials event 1','n trials event 2'}];
else
    cols = [cols, {'n trials'}];
end
if hasocc; cols = [cols, {'n occasions'}]; end
cols = [cols, {'Generalizability','G low','G high',...
    'Dependability','D low','D high','ICC (rel)','ICC (abs)'}];

% Preallocate the cell (the DoD helper's idiom) rather than growing it row by
% row; a 2D ngrid=25 surface is 625 rows per stratum.
nrows = 0;
for s = 1:numel(strata)
    if ndim == 2
        nrows = nrows + numel(strata(s).z1) * numel(strata(s).z2);
    else
        nrows = nrows + numel(strata(s).z1);
    end
end
data = cell(nrows, numel(cols));
r = 0;
for s = 1:numel(strata)
    st = strata(s);
    z1 = st.z1(:);
    if ndim == 2
        z2 = st.z2(:);
        for a = 1:numel(z1)
            for c = 1:numel(z2)
                r = r + 1;
                data(r,:) = local_surface_row(st, hasocc, hasev, hasev4, ndim, z1(a), z2(c), a, c);
            end
        end
    else
        for a = 1:numel(z1)
            r = r + 1;
            data(r,:) = local_surface_row(st, hasocc, hasev, hasev4, ndim, z1(a), [], a, []);
        end
    end
end
% cell2table with makeValidName matches the splits helper's export idiom.
surfaceTbl = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));

% ---- variance components: stack per-stratum tables (shared builder) ----------
% psyrat_dynrel_varcomp_concat is the single source of truth, also used by the
% dynrel viewer's concatVarcomp (stratum label + observed trial count).
vcTbl = psyrat_dynrel_varcomp_concat(summ);

% ---- per-participant reliability: subject-level variants only (shared) -------
% psyrat_dynrel_ssrel_concat is the single source of truth, also used by the
% dynrel viewer's concatSsrel. It returns [] when no stratum has an ssrel_table
% (non-subject-level variants), so the caller omits the field below.
ssrelTbl = psyrat_dynrel_ssrel_concat(summ);

% ---- estimand note for the export header (mirrors the viewer) ---------------
%the nonconcurrent DoD estimand clauses apply in EVERY family (gamma
%location-scale and, since the 2026-08-20 owner ruling, Gaussian), so the DoD
%tests run BEFORE the family branch: a Gaussian DoD result (summ.isgamma
%false) would otherwise fall through to the generic typical-person line and
%silently drop the "reference curve only" and per-participant clauses. The
%wording is identical across families - the per-cell log residual SD
%conditioning statement is true for both, because the scale submodel is
%log-linked in both (SCIENTIFIC_FORMULA_AUDIT.md section 26). The surfaces
%are the TYPICAL-PERSON estimand (all person effects at zero); presence-keyed
%on the dod_labels carrier. Mirrors psyrat_startview_dynrel.
isdodsum = isfield(summ,'dod_labels') && ~isempty(summ.dod_labels);
%B15 (2026-08-16): the DoD designs are identified by the ANALYSIS string,
%not by the labels carrier. A DoD result whose dod_labels is missing
%(damaged or pre-release stored file) must NOT fall through to a family
%branch's estimand sentence. Decline to assert an estimand instead,
%mirroring the string-keyed DEGRADED fallback in psyrat_provenance_lines.
%Mirrors psyrat_startview_dynrel.
isdodrun = any(strcmp(rel.analysis, ...
    {'ic_dodiff_dynrel_sserrvar','ic_dodiff_dynrel_sserrvar_trt'}));
if isdodsum
    estimand_note = ['Estimand: surface = typical person (all person '...
        'effects at zero); reference curve only - the per-participant '...
        'table is the primary output'];
    if ~isempty(ssrelTbl)
        estimand_note = [estimand_note ...
            '; per-participant table = conditional on each person''s own '...
            'per-cell log residual SDs, at their own dimension value(s)'];
    end
elseif isdodrun
    estimand_note = ['Estimand: cannot be determined from this stored '...
        'result - the dod_labels carrier is missing (damaged or '...
        'pre-release file); see the provenance lines for the '...
        'string-keyed disclosure'];
    if ~isempty(ssrelTbl)
        %same degraded arm as the surface clause: without the carrier the
        %per-participant sentence would otherwise claim a family branch's
        %conditional wording
        estimand_note = [estimand_note ...
            '; per-participant table = undetermined without the '...
            'dod_labels carrier'];
    end
elseif isfield(summ,'isgamma') && summ.isgamma
    %Neither gamma parameterization has a delta_p, so the Gaussian typical-person
    %label never applies: the surface is always marginal over the person scale
    %population, and any per-participant table reported alongside it is the
    %complementary CONDITIONAL read-out.
    %
    %WHAT they marginalize and condition on DIFFERS, so the wording is branched.
    %Under gammascale = 1 there is no log-residual submodel at all and the person
    %quantity is the log-nu dispersion. Under gammascale = 2 there IS one, and the
    %per-participant table conditions on each person's own log residual SD; their
    %own log-mean does not enter it, because the residual is decoupled from the
    %mean there (psyrat_ssrel_dynrel_gamma_ls). Mirrors psyrat_startview_dynrel.
    isls = isfield(summ,'isgamma_ls') && summ.isgamma_ls;
    if isls
        estimand_note = ['Estimand: surface = population average over person '...
            'residual-SD heterogeneity (gamma location-scale)'];
    else
        estimand_note = ['Estimand: surface = population average over person '...
            'dispersion (gamma estimand #2)'];
    end
    if ~isempty(ssrelTbl)
        if isls
            estimand_note = [estimand_note ...
                '; per-participant table = conditional on each person''s own '...
                'log residual SD'];
        else
            estimand_note = [estimand_note ...
                '; per-participant table = conditional on each person''s own '...
                'log-mean and log-nu'];
        end
    end
elseif isfield(summ,'marginal') && summ.marginal == 1
    estimand_note = 'Estimand: population average (lognormal)';
else
    estimand_note = 'Estimand: typical person (delta_p = 0)';
end
% Two-facet variants set strata(s).nocc; record the coefficient + occasion n'.
if ~isempty(strata) && ~isempty(strata(1).nocc)
    rtlabels = {'equivalence (occasion fixed)',...
        'stability (trial fixed)','equivalence and stability (both random)'};
    rt = reltype; if rt < 1 || rt > 3; rt = 3; end
    estimand_note = sprintf('%s; coefficient: %s; occasion n''=%d',...
        estimand_note, rtlabels{rt}, strata(1).nocc);
end
end

function row = local_surface_row(st, hasocc, hasev, hasev4, ndim, z1v, z2v, a, c)
%One flattened dynamic-reliability surface row: identity/design columns then the
%G/D/ICC values at grid index (a) [1D] or (a,c) [2D]. Factoring this out keeps
%the 1D and 2D paths from duplicating the cell assembly and the round-to-4dp.
%hasev (gamma difference variant) emits the two EXACT per-event n' the surface was
%computed at instead of the scalar st.obs, which is their mean; hasev4 (the DoD
%variants) emits all four per-cell counts for the same reason.
row = {char(string(st.label)), z1v};
if ndim == 2; row = [row, {z2v}]; end
if hasev4
    row = [row, {st.obs_ev(1), st.obs_ev(2), st.obs_ev(3), st.obs_ev(4)}];
elseif hasev
    row = [row, {st.obs_ev(1), st.obs_ev(2)}];
else
    row = [row, {st.obs}];
end
if hasocc; row = [row, {st.nocc}]; end
if ndim == 2
    row = [row, {round(st.G.pt(a,c),4), round(st.G.ll(a,c),4), round(st.G.ul(a,c),4), ...
        round(st.D.pt(a,c),4), round(st.D.ll(a,c),4), round(st.D.ul(a,c),4), ...
        round(st.ICCg.pt(a,c),4), round(st.ICCd.pt(a,c),4)}];
else
    row = [row, {round(st.G.pt(a),4), round(st.G.ll(a),4), round(st.G.ul(a),4), ...
        round(st.D.pt(a),4), round(st.D.ll(a),4), round(st.D.ul(a),4), ...
        round(st.ICCg.pt(a),4), round(st.ICCd.pt(a),4)}];
end
end

function [T, note] = local_sserr_table(psyrat_data)
%Flatten the subject-level (sserr) per-participant coefficients into ONE long
%table -- a render-free port of the psyrat_savesscoeffs CSV loop. Each stratum's
%relsummary.group(gloc).event(eloc).ssrel_table already holds the per-subject
%columns; this selects the exported subset (matching psyrat_savesscoeffs's column
%order), prepends the Group/Event stratum labels, and stacks them (addvars +
%vertcat, the local_dod_tables idiom). The dep_* columns hold the SELECTED
%coefficient (dependability or generalizability); the returned note (also written
%into the export header) records which, mirroring the viewer's export relabeling.
rs  = psyrat_data.relsummary;
rel = psyrat_data.rel;
[gnames, ngroups] = local_stratum_names(rel.groups);
[enames, nevents] = local_stratum_names(rel.events);

% The per-subject columns psyrat_savesscoeffs writes, in its order.
cols = {'id','ind2include','trls','dep_pt','dep_ll','dep_ul',...
    'icc_pt','icc_ll','icc_ul','sem_pt','sem_ll','sem_ul','ss_errvar'};

parts = cell(ngroups*nevents,1);
k = 0;
for gloc = 1:ngroups
    for eloc = 1:nevents
        t = rs.group(gloc).event(eloc).ssrel_table;
        t = t(:, cols);
        % prepend Event then Group so the final order is Group, Event, id, ...
        t = addvars(t, repmat(enames(eloc),height(t),1), ...
            'Before', 1, 'NewVariableNames', 'Event');
        t = addvars(t, repmat(gnames(gloc),height(t),1), ...
            'Before', 1, 'NewVariableNames', 'Group');
        k = k + 1;
        parts{k} = t;
    end
end
T = vertcat(parts{:});

if isfield(rs,'gcoeff_name') && strcmpi(rs.gcoeff_name,'gen')
    note = 'Subject-level coefficient (dep_* columns): Generalizability';
else
    note = 'Subject-level coefficient (dep_* columns): Dependability';
end
end

function [T, note] = local_curve_table(psyrat_data, ciperc, ntrials)
%Tabulate the reliability-as-a-function-of-trial-count curve -- the headless
%analog of the GUI trials plot, which the viewers only render. One-facet designs
%use psyrat_rel_sing; test-retest designs use psyrat_rel_trt: the SAME calc calls,
%components, and coefficient sourcing the plot functions use (psyrat_depvtrialsplot
%/ psyrat_trt_relvtrialsplot as driven by psyrat_relfigures), so the table matches
%the plot. That parity is an invariant to MAINTAIN, not an observation: until
%RC-10 the 'trt' branch and the plot matched each other only because both omitted
%nocc, so both disagreed with the summary tables in the same export. If you change
%a calc call here, change the matching plot function, and check both against the
%coefficient psyrat_relsummary computed.
%Flat long layout: one row per stratum (Group x Event) x n trials, with
%the selected coefficient's lower/point/upper limit over n = 1..ntrials.
rs  = psyrat_data.relsummary;
rel = psyrat_data.rel;
a   = rel.analysis;
[gnames, ngroups] = local_stratum_names(rel.groups);
[enames, nevents] = local_stratum_names(rel.events);

% Coefficient index (1 = dependability/absolute, 2 = generalizability/relative),
% sourced exactly as psyrat_relfigures feeds the plots. RC-31: that is gcoeff for
% EVERY design tabulated here, trt_diff included. This table is per-condition in
% all of them -- see local_curve_series, where each branch reads per-condition
% components -- and per-condition output follows gcoeff, the contract stated at
% psyrat_trt_coeflabel.m:94-99. trt_diff was special-cased to rs.diffgcoeff, which
% applied the difference-score selector to a per-condition quantity and put this
% export out of step with the per-condition summary rows in the same export. The
% genuine difference-score curve is a separate table (local_diffcurve_table); that
% one follows diffgcoeff.
gc = rs.gcoeff;
if gc == 2
    clab = 'Generalizability';
else
    clab = 'Dependability';
end

cols = {'Group','Event','n_trials', ...
    [clab '_LowerLimit'], [clab '_PointEstimate'], [clab '_UpperLimit']};
data = cell(ngroups*nevents*ntrials, numel(cols));
r = 0;
for gloc = 1:ngroups
    for eloc = 1:nevents
        [ll,pt,ul] = local_curve_series(a, rs, gc, gloc, eloc, ntrials, ciperc);
        for i = 1:ntrials
            r = r + 1;
            data(r,:) = {gnames{gloc}, enames{eloc}, i, ...
                round(ll(i),4), round(pt(i),4), round(ul(i),4)};
        end
    end
end
T = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));

%RC-31. This series is PER-CONDITION for every analysis tabulated here, the two
%difference designs included, so it must not be labeled a difference score. Both
%were: ic_diff reads sd_id, which psyrat_relsummary.m:4553 builds as
%sqrt(squeeze(idvar(:,eloc,eloc))) -- the DIAGONAL of the cross-condition
%covariance, i.e. one condition -- and trt_diff reads the per-condition SD draws
%the same way. A reader who took the old note at its word would have read a
%per-condition coefficient as the reliability of the contrast they care about.
%On a difference design the qualifier still earns its place, because a genuine
%difference-score curve is exported alongside this one and the two must be
%distinguishable in a standalone file.
if strcmp(a,'ic_diff') || strcmp(a,'trt_diff')
    note = ['Reliability-vs-n curve coefficient: ' clab ...
        ' (per condition, not the difference score)'];
else
    note = ['Reliability-vs-n curve coefficient: ' clab];
end
end

function [T, note] = local_diffcurve_table(psyrat_data, ciperc, ntrials)
%Tabulate the DIFFERENCE-SCORE reliability-as-a-function-of-trial-count curve for
%the two-facet difference design -- the headless analog of
%psyrat_trt_diffvtrialsplot, and the difference-score companion to
%local_curve_table's per-condition series. Figure/table parity is the same
%invariant to MAINTAIN here as there: same kernel, same components, same
%coefficient sourcing. If you change one, change the other.
%
%Flat long layout: one row per Group x n. There is no Event column, and that is
%the point -- a difference between two conditions is not stratified by condition.
%
%Coefficient: diffgcoeff, the separately selected difference-score error type,
%per the contract at psyrat_trt_coeflabel.m:94-99. The per-condition table
%follows gcoeff.
rs  = psyrat_data.relsummary;
rel = psyrat_data.rel;
[gnames, ngroups] = local_stratum_names(rel.groups);

if isfield(rs,'diffgcoeff') && ~isempty(rs.diffgcoeff)
    gc = rs.diffgcoeff;
else
    gc = 1;
end
if gc == 2
    clab = 'Generalizability';
else
    clab = 'Dependability';
end

cols = {'Group','n_trials', ...
    [clab '_LowerLimit'], [clab '_PointEstimate'], [clab '_UpperLimit']};
data = cell(ngroups*ntrials, numel(cols));
r = 0;
for gloc = 1:ngroups
    [ll,pt,ul] = local_diffcurve_series(rs, gc, gloc, ntrials, ciperc);
    for i = 1:ntrials
        r = r + 1;
        data(r,:) = {gnames{gloc}, i, ...
            round(ll(i),4), round(pt(i),4), round(ul(i),4)};
    end
end
T = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));

%The estimand caveat rides in the note because it cannot be recovered from the
%columns: n is the count in BOTH conditions, while the headline difference score
%in the overall table is computed at the observed (possibly unequal) per-
%condition counts, so the two need not coincide.
note = ['Difference-score reliability-vs-n curve coefficient: ' clab ...
    ' (n = trials per condition, equal in both; the overall table''s ' ...
    'difference score uses the observed per-condition counts)'];
end

function [ll,pt,ul] = local_diffcurve_series(rs, gc, gloc, ntrials, ciperc)
%Compute one group's difference-score curve over n = 1..ntrials, replicating the
%exact calc call psyrat_trt_diffvtrialsplot makes. Returns column vectors.
%
%psyrat_diffrel_trt takes a scalar pair obs = [nX nY] and has no range mode
%(unlike psyrat_rel_sing/psyrat_rel_trt, which the per-condition series above
%calls once for the whole range), so this loops and adapts the returned struct to
%(ll,pt,ul). obs = [n n] projects equal counts in both conditions -- the same
%convention psyrat_rel_diffdynrel_trt:178 uses, and the only one on which n
%remains a plain trial count.
if ~isfield(rs,'group') || ~isfield(rs.group,'diffcomp')
    error('psyrat_report:nodiffcomp', ...
        ['relsummary.group.diffcomp is not populated, so the difference-score ' ...
        'components are unavailable. How to fix: run psyrat_relsummary with ' ...
        '''analysis'',''trt_diff'' on these data before calling psyrat_report.']);
end

if isfield(rs,'nocc') && ~isempty(rs.nocc)
    nocc = rs.nocc;
else
    nocc = 1;
end

if gc == 2
    diff_est = 'gen';
else
    diff_est = 'dep';
end

dc = rs.group(gloc).diffcomp;
blockargs = {'bp',dc.bp,'bpi',dc.bpi,'bpo',dc.bpo, ...
    'bt',dc.bt,'bo',dc.bo,'boi',dc.boi, ...
    'er_var',dc.er_var,'er_cov',dc.er_cov};

ll = zeros(ntrials,1); pt = zeros(ntrials,1); ul = zeros(ntrials,1);

%One warning per CALL from the kernel, and this loops it ntrials times, so
%suppress the shared identifier, aggregate, and report once -- the idiom
%psyrat_admissibility_warn's own help prescribes for loop callers and
%psyrat_rel_diffdynrel_trt.m:152-202 already uses. At obs = [n n] every adjusted
%block reduces to raw/n, so admissibility here does not vary with n: this either
%never fires or fires at all ntrials points identically, which is exactly what
%aggregation is for. See the fuller note in psyrat_trt_diffvtrialsplot.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];
for n = 1:ntrials
    ds = psyrat_diffrel_trt(blockargs{:}, ...
        'obs',[n n],'nocc',[nocc nocc], ...
        'reltype',rs.reltype,'est',diff_est,'CI',ciperc);
    ll(n) = ds.ll; pt(n) = ds.pt; ul(n) = ds.ul;
    if isfield(ds,'admissibility')
        adm_all = [adm_all, ds.admissibility]; %#ok<AGROW>
    end
end
warning(adm_ws);
psyrat_admissibility_warn(adm_all);
end

function sigi = local_sigtrl(d, wp_resid)
%Trial main effect (sigma_i) draws for one one-facet stratum, returned as an SD
%vector shaped like the residual SD it will be combined with.
%
%The component is stored under a different field name per analysis, because each
%branch of psyrat_relsummary builds its own data struct: the one-facet model
%stores an SD directly (sig_trl); the single-session difference stores a per-
%event SD (sd_trl); the subject-level difference stores the trial covariance and
%the per-event value is its diagonal, a VARIANCE, hence the sqrt. Checked in that
%order because the difference structs carry trl_varcov as well, and sd_trl is the
%unambiguous per-event SD.
%
%Legacy result structs produced before the trial main effect was modeled carry
%none of these; those return zeros, which reduces the coefficient to the prior
%behavior (sigma_i = 0, dependability == generalizability) rather than erroring.
sigi = zeros(size(wp_resid));
if isfield(d,'sig_trl') && isfield(d.sig_trl,'raw') && ~isempty(d.sig_trl.raw)
    sigi = d.sig_trl.raw;
elseif isfield(d,'sd_trl') && isfield(d.sd_trl,'raw') && ~isempty(d.sd_trl.raw)
    sigi = cell2mat(d.sd_trl.raw);
elseif isfield(d,'trl_varcov') && isfield(d.trl_varcov,'raw') && ...
        ~isempty(d.trl_varcov.raw)
    sigi = sqrt(max(cell2mat(d.trl_varcov.raw),0));
end
sigi = sigi(:);
end

function [ll,pt,ul] = local_curve_series(a, rs, gc, gloc, eloc, ntrials, ciperc)
%Compute one stratum's reliability curve over n = 1..ntrials, replicating the
%exact calc call the matching plot function makes. Returns column vectors.
d = rs.data.g(gloc).e(eloc);
switch a
    case 'ic'
        % wp is the TOTAL within-person SD and sigma_i is passed separately so
        % psyrat_rel_sing can subtract it back out for the relative coefficient.
        % Both are required: passing 'i' against a residual-only wp leaves
        % dependability unchanged and corrupts generalizability (RC-02).
        sigi = local_sigtrl(d, d.sig_e.raw);
        [ll,pt,ul] = psyrat_rel_sing('gcoeff',gc,'metric','global',...
            'bp',d.sig_u.raw,'wp',sqrt(d.sig_e.raw.^2 + sigi.^2),'i',sigi,...
            'obs',[1 ntrials],'CI',ciperc);
    case 'ic_diff'
        % RC-31. PER-CONDITION between-/within-person SDs, not difference-score
        % ones as this comment used to claim: psyrat_relsummary.m:4553 builds
        % sd_id as sqrt(squeeze(idvar(:,eloc,eloc))), the diagonal of the
        % cross-condition covariance. Same components psyrat_depvtrialsplot uses,
        % and correctly driven by gcoeff; only the description was wrong.
        bp = cell2mat(d.sd_id.raw);
        wp = exp(cell2mat(d.b_sigma.raw));
        sigi = local_sigtrl(d, wp);
        [ll,pt,ul] = psyrat_rel_sing('gcoeff',gc,'metric','global',...
            'bp',bp,'wp',sqrt(wp.^2 + sigi.^2),'i',sigi,...
            'obs',[1 ntrials],'CI',ciperc);
    case 'trt'
        % RC-10. Plain test-retest at the SELECTED occasion estimand, matching
        % the trt_diff branch below and psyrat_trt_relvtrialsplot. This branch
        % previously passed no nocc, so it defaulted to n'_o = 1 while the
        % summary tables in the same export used the k the user chose. The
        % comment that stood here justified that by pointing at the plot
        % function -- which had the same omission -- and called n'_o = 1 "the
        % toolbox's deliberate default". The deliberate default is
        % noccmode = 1; when the user selects the multi-occasion composite,
        % every other consumer honors k and this curve now does too.
        if isfield(rs,'nocc') && ~isempty(rs.nocc)
            nocc = rs.nocc;
        else
            nocc = 1;
        end
        [ll,pt,ul] = psyrat_rel_trt('gcoeff',gc,'reltype',rs.reltype,...
            'bp',d.sig_id.raw,'bo',d.sig_occ.raw,'bt',d.sig_trl.raw,...
            'txp',d.sig_trlxid.raw,'oxp',d.sig_occxid.raw,'txo',d.sig_trlxocc.raw,...
            'err',d.sig_err.raw,'obs',[1 ntrials],'nocc',nocc,'CI',ciperc);
    case 'trt_diff'
        % RC-31. PER-CONDITION two-facet reliability at the selected occasion
        % estimand, per psyrat_depvtrialsplot's trt_diff branch. This calls
        % psyrat_rel_trt, the plain two-facet kernel, on the per-condition
        % components (sig_id and friends are the covariance diagonals stored as
        % SD draws) -- it is NOT a difference kernel, and the comment that stood
        % here called it "two-facet difference", which is what made feeding it
        % diffgcoeff look reasonable. The genuine difference-score curve is
        % local_diffcurve_series.
        if isfield(rs,'nocc') && ~isempty(rs.nocc)
            nocc = rs.nocc;
        else
            nocc = 1;
        end
        [ll,pt,ul] = psyrat_rel_trt('gcoeff',gc,'reltype',rs.reltype,...
            'bp',d.sig_id.raw,'bo',d.sig_occ.raw,'bt',d.sig_trl.raw,...
            'txp',d.sig_trlxid.raw,'oxp',d.sig_occxid.raw,'txo',d.sig_trlxocc.raw,...
            'err',d.sig_err.raw,'obs',[1 ntrials],'nocc',nocc,'CI',ciperc);
    otherwise
        error('psyrat_report:curve', ...
            'No reliability-vs-n curve is defined for analysis ''%s''.', a);
end
ll = ll(:); pt = pt(:); ul = ul(:);
end

function [names, n] = local_stratum_names(field)
%Stratum labels for a group/event field, mirroring the plot functions' handling:
%the literal 'none' sentinel means a single unlabeled stratum (''); otherwise the
%field is the list of stratum names. Returns a cellstr and its count.
if strcmpi(field,'none')
    names = {''};
    n = 1;
else
    names = cellstr(string(field(:)));
    n = numel(names);
end
end

function name = local_dataset_name(psyrat_data)
%Dataset identity for the export header. Prefer the saved data's recorded source
%filename (rel.filename, matching the GUI viewers); fall back to the loaded
%long-format table's source (raw.filename) so in-memory / scripted runs that
%never set rel.filename still record a dataset; then the saved .psyrat name. ''
%when none is available, so the caller omits the Dataset line.
name = '';
if isfield(psyrat_data,'rel') && isstruct(psyrat_data.rel) && ...
        isfield(psyrat_data.rel,'filename') && ~isempty(psyrat_data.rel.filename)
    name = char(string(psyrat_data.rel.filename));
    return;
end
if isfield(psyrat_data,'raw') && isstruct(psyrat_data.raw) && ...
        isfield(psyrat_data.raw,'filename') && ~isempty(psyrat_data.raw.filename) && ...
        (ischar(psyrat_data.raw.filename) || isstring(psyrat_data.raw.filename))
    name = char(string(psyrat_data.raw.filename));
    return;
end
if isfield(psyrat_data,'proc') && isstruct(psyrat_data.proc) && ...
        isfield(psyrat_data.proc,'savename') && ~isempty(psyrat_data.proc.savename)
    name = char(string(psyrat_data.proc.savename));
end
end

function v = local_toolbox_version(psyrat_data)
%Toolbox version for the header: the saved data's recorded version if present,
%else the running toolbox version, else 'unknown'.
if isfield(psyrat_data,'ver') && ~isempty(psyrat_data.ver)
    v = psyrat_data.ver;
    return;
end
try
    v = psyrat_defineversion();
catch
    v = 'unknown';
end
if isempty(v)
    v = 'unknown';
end
end

function head = local_table_header(psyrat_data, tabletag, provlines, extralines)
%Build the header block written above a table: toolbox + dataset identity, the
%table tag, the shared provenance lines, then any analysis-specific extra lines
%(e.g. the data-splits generalizability lower-bound caveat).
if nargin < 4 || isempty(extralines)
    extralines = {};
end
head = cell(0,1);
head{end+1,1} = sprintf('Table Generated on %s', char(datetime('now')));
% Always identify the toolbox version: the version the data were created with
% when recorded, otherwise the running toolbox version, otherwise 'unknown'.
head{end+1,1} = sprintf('PsyRAT Toolbox v%s', local_toolbox_version(psyrat_data));
dataset = local_dataset_name(psyrat_data);
if ~isempty(dataset)
    head{end+1,1} = sprintf('Dataset: %s', dataset);
end
head{end+1,1} = sprintf('Table: %s', tabletag);
head = [head; provlines(:); extralines(:)];
head{end+1,1} = '';
end
