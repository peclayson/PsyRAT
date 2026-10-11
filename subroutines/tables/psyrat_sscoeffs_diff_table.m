function [T, info] = psyrat_sscoeffs_diff_table(psyrat_data)
%Validate and assemble the per-participant difference-score rows of a
%subject-level difference analysis (ic_diff_sserrvar) for export (G107).
%
%[T, info] = psyrat_sscoeffs_diff_table(psyrat_data)
%
%Both subject-coefficient exports call this before anything is written: the
%Save Subject Coefficients callback (psyrat_variancet) right after its Cancel
%guard, and psyrat_report before it builds any table. Every group is checked
%first. A missing or malformed difference table raises an error that names the
%group and how to rebuild the summary; event rows are never exported in its
%place, and nothing is written.
%
%The rows are the stored summary tables (relsummary.group(g).diffscore
%.ssrel_table), copied without recalculation, group by group in the stored row
%order, each row carrying the run-level labels it needs to stand alone. A typed
%zero-row table is valid: its group contributes no rows and gets an explicit
%empty status. Rows whose participant is excluded, or whose values are NaN,
%are kept as stored.
%
%Input
% psyrat_data - toolbox data structure after psyrat_relsummary
%   ('sing_diff_sserr'); rel.analysis must be 'ic_diff_sserrvar'.
%
%Output
% T    - table, one row per stored difference row, in these variables:
%        Group ('' for an ungrouped run), ScoreType ('difference'), EventA and
%        EventB (the fitted events, in rel.events order), Direction
%        ('<EventA> minus <EventB>'), Coefficient (Dependability or
%        Generalizability), id, ind2include (the DIFFERENCE-score flag, which
%        decides retention), Cutoff, CutoffBound, trls1 and trls2 (the
%        participant's own trial count in EventA and in EventB), dep_pt, dep_ll
%        and dep_ul (the selected difference coefficient; the dep_* convention
%        of tables.sscoeffs), IntervalKind, IntervalWidth (a proportion), icc_*,
%        sem_*, and diff_rel_errvar (the stored ss_errvar of the difference
%        table).
% info - struct: groups (stratum labels), sectionLabels (one per group),
%        rows (per group, the row numbers of T), emptyStatus (per group, '' or
%        the empty-status line), notes (run-level lines every difference export
%        carries in its header), and eventA, eventB, direction, coefficient,
%        intervalKind, intervalWidth, cutoff and bound.

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

route = ['How to fix: summarize the saved draws again (open the result in View ', ...
    'Results, or call psyrat_relsummary with ''analysis'',''sing_diff_sserr''), ', ...
    'then export again. No refit is needed.'];

rel = psyrat_data.rel;
if ~isfield(rel,'analysis') || ~strcmp(rel.analysis,'ic_diff_sserrvar')
    error('psyrat:sscoeffsDiff:analysis', ...
        ['Difference-score subject rows exist only for the subject-level ', ...
        'difference analysis (ic_diff_sserrvar).']);
end

%strata, named the way the subject-coefficient exports name them
if isequal(strcmpi(rel.groups,'none'), true)
    gnames = {''};
else
    gnames = cellstr(string(rel.groups(:)))';
end
ngroups = numel(gnames);
enames = cellstr(string(rel.events(:)))';
if numel(enames) ~= 2
    error('psyrat:sscoeffsDiff:malformed', ...
        'This result does not record exactly two events, so it has no difference score. %s', route);
end

%---- validate everything before assembling anything ------------------------
if ~isfield(psyrat_data,'relsummary') || ~isstruct(psyrat_data.relsummary) || ...
        ~isfield(psyrat_data.relsummary,'group')
    error('psyrat:sscoeffsDiff:missing', ...
        'This result carries no summary, so no difference rows can be exported and nothing was saved. %s', route);
end
rs = psyrat_data.relsummary;
coef = local_coefficient(rs);
problem = '';
if isempty(coef)
    problem = 'no difference coefficient choice (diffgcoeff) is recorded';
elseif ~isfield(rs,'depcutoff') || ~isnumeric(rs.depcutoff) || ~isscalar(rs.depcutoff) || ...
        ~isfinite(rs.depcutoff)
    problem = 'no reliability cutoff is recorded';
elseif ~isfield(rs,'meascutoff') || ~(ischar(rs.meascutoff) || isstring(rs.meascutoff)) || ...
        strlength(string(rs.meascutoff)) == 0
    problem = 'no cutoff bound is recorded';
elseif ~isfield(rs,'ciperc') || ~isnumeric(rs.ciperc) || ~isscalar(rs.ciperc) || ...
        ~(rs.ciperc > 0 && rs.ciperc < 1)
    problem = 'no credible-interval width is recorded';
end
if ~isempty(problem)
    error('psyrat:sscoeffsDiff:malformed', ...
        'The stored summary is incomplete (%s), so nothing was saved. %s', problem, route);
end
for g = 1:ngroups
    where = local_where(gnames{g});
    if numel(rs.group) < g || ~isfield(rs.group(g),'diffscore') || ...
            ~isstruct(rs.group(g).diffscore) || ~isscalar(rs.group(g).diffscore) || ...
            ~isfield(rs.group(g).diffscore,'ssrel_table')
        error('psyrat:sscoeffsDiff:missing', ...
            ['No subject-level difference-score table was found for %s, so the ', ...
            'difference rows cannot be exported and nothing was saved. %s'], where, route);
    end
    problem = local_table_problem(rs.group(g).diffscore.ssrel_table);
    if ~isempty(problem)
        error('psyrat:sscoeffsDiff:malformed', ...
            ['The difference-score table for %s is malformed (%s), so nothing ', ...
            'was saved. %s'], where, problem, route);
    end
end

%---- assemble --------------------------------------------------------------
direction = sprintf('%s minus %s', enames{1}, enames{2});
kind = ['equal-tailed ' psyrat_interval_label(rel,'lower')];
bound = char(string(rs.meascutoff));
names = {'Group','ScoreType','EventA','EventB','Direction','Coefficient','id', ...
    'ind2include','Cutoff','CutoffBound','trls1','trls2','dep_pt','dep_ll','dep_ul', ...
    'IntervalKind','IntervalWidth','icc_pt','icc_ll','icc_ul','sem_pt','sem_ll', ...
    'sem_ul','diff_rel_errvar'};
parts = cell(1,ngroups);
rows = cell(1,ngroups);
labels = cell(1,ngroups);
emptyStatus = repmat({''},1,ngroups);
next = 0;
for g = 1:ngroups
    D = rs.group(g).diffscore.ssrel_table;
    n = height(D);
    parts{g} = table(repmat(gnames(g),n,1), repmat({'difference'},n,1), ...
        repmat(enames(1),n,1), repmat(enames(2),n,1), repmat({direction},n,1), ...
        repmat({coef},n,1), cellstr(string(D.id)), logical(D.ind2include), ...
        repmat(rs.depcutoff,n,1), repmat({bound},n,1), D.trls1, D.trls2, ...
        D.dep_pt, D.dep_ll, D.dep_ul, repmat({kind},n,1), repmat(rs.ciperc,n,1), ...
        D.icc_pt, D.icc_ll, D.icc_ul, D.sem_pt, D.sem_ll, D.sem_ul, D.ss_errvar, ...
        'VariableNames', names);
    rows{g} = next + (1:n);
    next = next + n;
    base = sprintf('Difference Score (%s)', direction);
    if isempty(gnames{g})
        labels{g} = base;
    else
        labels{g} = [gnames{g} ' - ' base];
    end
    if n == 0
        if isempty(gnames{g})
            emptyStatus{g} = ['The stored difference-score table has zero rows ', ...
                '(empty), so no difference rows are written.'];
        else
            emptyStatus{g} = sprintf(['Group %s: the stored difference-score table ', ...
                'has zero rows (empty), so no difference rows are written for this ', ...
                'group.'], gnames{g});
        end
    end
end
T = vertcat(parts{:});

%run-level notes: what the rows are, which flag decides retention, and the
%three quantities whose scale is easy to misread
if strcmp(coef,'Dependability')
    errword = 'absolute';
else
    errword = 'relative';
end
notes = { ...
    sprintf(['Difference-score rows: the %s coefficient of each participant''s %s ', ...
    'difference score (posterior mean, %g%% %s), at the participant''s own trial ', ...
    'count in each event. The coefficient does not depend on the direction of the ', ...
    'difference.'], lower(coef), direction, 100*rs.ciperc, kind), ...
    ['Difference-score flag: it compares the difference-score coefficient with the ', ...
    'reliability cutoff at the cutoff bound, and it is the flag that decides retention ', ...
    'and the Save IDs lists. Each event row''s flag compares that event''s own ', ...
    'coefficient with the same cutoff. Flags compare unrounded values. The Save ', ...
    'Subject Coefficients file prints three decimals (the cutoff four), and ', ...
    'psyrat_report returns the unrounded values in tables.sscoeffs_diff, so a value ', ...
    'printed at the cutoff can carry either flag.'], ...
    sprintf(['Difference-score SEM: the square root of the error variance of the ', ...
    'selected coefficient, absolute error for dependability and relative error for ', ...
    'generalizability (%s error here), at the participant''s own trial counts. A draw ', ...
    'with a negative error variance has no SEM; the posterior mean is then NaN and the ', ...
    'limits use the remaining draws.'], errword), ...
    ['Difference error variance: always the relative error variance of the difference ', ...
    'between the participant''s two event means at the participant''s own trial ', ...
    'counts, whichever coefficient is selected. It is not the event rows'' Error ', ...
    'Variance (ss_errvar), which is the participant''s single-trial relative residual ', ...
    'variance. Its covariance term is divided by the harmonic mean of the two counts, ', ...
    'following Rocha et al. (2026) Table 6. When the counts differ, that divisor is ', ...
    'the adopted convention, not an exact result.'], ...
    sprintf(['Difference-score ICC (n''=1): at one observation per event, with %s ', ...
    'error. It is a derived quantity, not an entry of Rocha et al. (2026) Table 6.'], ...
    errword)};

info = struct();
info.groups = gnames;
info.sectionLabels = labels;
info.rows = rows;
info.emptyStatus = emptyStatus;
info.notes = notes;
info.eventA = enames{1};
info.eventB = enames{2};
info.direction = direction;
info.coefficient = coef;
info.intervalKind = kind;
info.intervalWidth = rs.ciperc;
info.cutoff = rs.depcutoff;
info.bound = bound;
end

function coef = local_coefficient(rs)
%The difference coefficient the summary selected ('' when none is recorded).
coef = '';
if isfield(rs,'diffgcoeff_name') && (ischar(rs.diffgcoeff_name) || isstring(rs.diffgcoeff_name))
    switch lower(char(rs.diffgcoeff_name))
        case 'dep'
            coef = 'Dependability';
        case 'gen'
            coef = 'Generalizability';
    end
elseif isfield(rs,'diffgcoeff') && isnumeric(rs.diffgcoeff) && isscalar(rs.diffgcoeff)
    if rs.diffgcoeff == 1
        coef = 'Dependability';
    elseif rs.diffgcoeff == 2
        coef = 'Generalizability';
    end
end
end

function where = local_where(gname)
if isempty(gname)
    where = 'this run';
else
    where = sprintf('group "%s"', gname);
end
end

function problem = local_table_problem(D)
%'' for a usable stored difference table, else what is wrong with it. Types and
%shapes are checked column by column; a typed zero-row table passes.
problem = '';
if ~istable(D)
    problem = 'it is not a table';
    return;
end
required = {'id','trls1','trls2','dep_pt','dep_ll','dep_ul','icc_pt','icc_ll', ...
    'icc_ul','sem_pt','sem_ll','sem_ul','ss_errvar','ind2include'};
absent = required(~ismember(required, D.Properties.VariableNames));
if ~isempty(absent)
    problem = ['it has no ' strjoin(absent, ', ') ' column'];
    return;
end
v = D.id;
if ~((iscellstr(v) || (isstring(v) && ~any(ismissing(v)))) && size(v,2) == 1)
    problem = 'id is not a single column of text';
    return;
end
numeric = required(2:13);
for k = 1:numel(numeric)
    v = D.(numeric{k});
    if ~(isnumeric(v) && isreal(v) && ismatrix(v) && size(v,2) == 1)
        problem = sprintf('%s is not a single numeric column', numeric{k});
        return;
    end
end
for k = 1:2
    v = D.(numeric{k});
    if any(~isfinite(v) | v < 0 | v ~= round(v))
        problem = sprintf('%s is not a trial count', numeric{k});
        return;
    end
end
v = D.ind2include;
if ~((islogical(v) || (isnumeric(v) && all(v == 0 | v == 1))) && size(v,2) == 1)
    problem = 'ind2include is not a single column of 0/1 flags';
end
end
