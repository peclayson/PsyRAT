function manifestPath = G107_write_reopen_artifacts(outdir)
%Write analysis-8 subject exports for a reopen OUTSIDE MATLAB (G107's
%independent reopen of the CSV and XLSX files, and G109's round trip of
%special IDs through a non-MATLAB CSV reader).
%
% manifestPath = G107_write_reopen_artifacts(outdir)
%
%Runs, with no fit, on six draw-backed fixtures (single group at
%dependability and at generalizability, the grouped fixture, the special IDs,
%and two typed zero-row difference tables injected at the summary seam: the
%single group, whose difference CSV and XLSX therefore hold no rows, and the
%first of two groups): the real Save Subject Coefficients callback (its save
%dialog answered by the G107 uiputfile double) and the real psyrat_report in
%CSV and in XLSX. Each difference section in the manifest says whether it is
%empty, so the reader can require the explicit empty status and no rows.
%It then writes manifest.json with what an independent reader must find:
%every event row from the stored summary, every difference row from the
%stored summary, and beside it the posterior mean and limits recomputed from
%the raw draws (G107ExportFixtures.diffReference). It asserts nothing: a
%callback or report error is recorded in the manifest, and the Python reader
%(csv and openpyxl, run after MATLAB exits) is the check. Before the fix it
%therefore documents the omission; after the fix it is the outside-MATLAB
%acceptance for (g).
%
%OUTDIR must not exist, so an earlier artifact set is never overwritten.
%Development-side only; nothing in the toolbox calls this.
%
%Input
% outdir - a path that does not exist yet.
%
%Output
% manifestPath - full path of the manifest.json written into OUTDIR.

if isfolder(outdir) || isfile(outdir)
    error('G107:artifactsExist', ...
        ['Refusing to write into an existing path (%s). How to fix: pass a ' ...
        'new folder name.'], outdir);
end
mkdir(outdir);

oldVisible = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
shadow = G107ExportFixtures.installShadows();
cleaner = onCleanup(@() localRestore(shadow, oldVisible)); %#ok<NASGU>

%{name, draws, coefficient (1 dependability, 2 generalizability), bound
%(1 lower, 2 point, 3 upper), interval width, groups whose stored difference
%table is replaced by its typed zero-row version after the summary}
scen = { ...
    'single_dependability',    G107ExportFixtures.baseRel(),       1, 2, 0.80, []; ...
    'single_generalizability', G107ExportFixtures.baseRel(),       2, 1, 0.80, []; ...
    'grouped',                 G107ExportFixtures.groupedRel(),    1, 2, 0.80, []; ...
    'special_ids',             G107ExportFixtures.specialIdRel(),  1, 2, 0.80, []; ...
    'typed_empty',             G107ExportFixtures.baseRel(),       1, 2, 0.80, 1; ...
    'typed_empty_grouped',     G107ExportFixtures.groupedRel(),    1, 2, 0.80, 1};

[guiKeys, guiNames] = G107ExportFixtures.guiDiffColumns();
manifest = struct();
manifest.generated_by = 'tests/helpers/G107_write_reopen_artifacts.m';
manifest.gui_diff_keys = guiKeys;
manifest.gui_diff_columns = guiNames;
manifest.headless_diff_columns = G107ExportFixtures.headlessDiffColumns();
manifest.headless_event_columns = G107ExportFixtures.headlessEventColumns();
manifest.scenarios = cell(1, size(scen, 1));

for s = 1:size(scen, 1)
    name = scen{s, 1};
    gc = scen{s, 3};
    meas = scen{s, 4};
    ci = scen{s, 5};
    cut = G107ExportFixtures.mixedCutoffValue(scen{s, 2}, gc, ci, meas);
    pd = G107ExportFixtures.summarize(scen{s, 2}, gc, ci, cut, meas);
    for g = scen{s, 6}
        D = pd.relsummary.group(g).diffscore.ssrel_table;
        pd.relsummary.group(g).diffscore.ssrel_table = D([], :);
    end
    pd.proc.savename = [name '.psyrat'];

    e = struct();
    e.scenario = name;

    %GUI: the real Save Subject Coefficients callback
    e.gui_csv = fullfile(outdir, [name '_gui.csv']);
    before = findall(groot, 'Type', 'figure');
    psyrat_variancet('psyrat_data', pd, 'gui', 1);
    after = findall(groot, 'Type', 'figure');
    figs = after(~ismember(after, before));
    button = localButton(figs, 'Save Subject');
    G107ExportFixtures.resetDialogs({e.gui_csv}, 'decline');
    out = G107ExportFixtures.fire(button);
    e.gui_error = out.errorText;
    after = findall(groot, 'Type', 'figure');
    delete(after(~ismember(after, before)));

    %headless: the real psyrat_report, both formats, each to its own folder
    for fmt = {'.csv', '.xlsx'}
        tag = fmt{1}(2:end);
        d = fullfile(outdir, [name '_report_' tag]);
        e.(['report_' tag '_dir']) = d;
        e.(['report_' tag '_expected']) = struct( ...
            'variance', fullfile(d, [name '_variance' fmt{1}]), ...
            'sscoeffs', fullfile(d, [name '_sscoeffs' fmt{1}]), ...
            'sscoeffs_diff', fullfile(d, [name '_sscoeffs_diff' fmt{1}]));
        try
            r = psyrat_report(pd, 'outdir', d, 'format', fmt{1});
            e.(['report_' tag '_files']) = r.files;
            e.(['report_' tag '_error']) = '';
        catch ME
            e.(['report_' tag '_files']) = {};
            e.(['report_' tag '_error']) = ME.message;
        end
    end

    e = localExpected(e, pd, gc, ci, cut, meas);
    manifest.scenarios{s} = e;
end

manifestPath = fullfile(outdir, 'manifest.json');
fid = fopen(manifestPath, 'w');
fwrite(fid, jsonencode(manifest, 'PrettyPrint', true));
fclose(fid);
fprintf('G107 reopen artifacts written: %s\n', manifestPath);
end


function e = localExpected(e, pd, gc, ci, cut, meas)
%What an independent reader must find in the files. NaN becomes JSON null.
rs = pd.relsummary;
e.coefficient_label = G107ExportFixtures.coefLabel(gc);
e.event_coefficient_label = G107ExportFixtures.eventCoefLabel(pd);
e.ci = ci;
e.cutoff = cut;
e.bound = rs.meascutoff;
e.event_header = G107ExportFixtures.guiEventHeader(e.event_coefficient_label);
e.provenance = psyrat_provenance_lines(pd)';

g = G107ExportFixtures.groupLabels(pd);
ev = G107ExportFixtures.eventNames(pd);
labels = G107ExportFixtures.eventSectionLabels(pd);
sections = {};
k = 0;
for gi = 1:numel(rs.group)
    for ei = 1:numel(rs.group(gi).event)
        k = k + 1;
        T = rs.group(gi).event(ei).ssrel_table;
        rows = cell(1, height(T));
        for j = 1:height(T)
            rows{j} = struct('id', char(string(T.id(j))), ...
                'include', double(T.ind2include(j)), 'trls', T.trls(j), ...
                'values', [T.dep_pt(j), T.dep_ll(j), T.dep_ul(j), ...
                T.icc_pt(j), T.icc_ll(j), T.icc_ul(j), ...
                T.sem_pt(j), T.sem_ll(j), T.sem_ul(j), T.ss_errvar(j)]);
        end
        sections{end+1} = struct('label', labels{k}, 'group', g{gi}, ...
            'event', ev{ei}, 'rows', {rows}); %#ok<AGROW>
    end
end
e.event_sections = sections;

dlabels = G107ExportFixtures.diffSectionLabels(pd);
dsections = cell(1, numel(rs.group));
for gi = 1:numel(rs.group)
    D = rs.group(gi).diffscore.ssrel_table;
    R = G107ExportFixtures.diffReference(pd, gi, gc, ci, cut, meas);
    rows = cell(1, height(D));
    for j = 1:height(D)
        id = char(string(D.id(j)));
        kk = find(strcmp(R.id, id), 1);
        rows{j} = struct('id', id, 'include', double(D.ind2include(j)), ...
            'trls1', D.trls1(j), 'trls2', D.trls2(j), ...
            'pt', D.dep_pt(j), 'll', D.dep_ll(j), 'ul', D.dep_ul(j), ...
            'icc_pt', D.icc_pt(j), 'icc_ll', D.icc_ll(j), 'icc_ul', D.icc_ul(j), ...
            'sem_pt', D.sem_pt(j), 'sem_ll', D.sem_ll(j), 'sem_ul', D.sem_ul(j), ...
            'errvar', D.ss_errvar(j), ...
            'ref_pt', R.pt(kk), 'ref_ll', R.ll(kk), 'ref_ul', R.ul(kk), ...
            'ref_trls1', R.trls1(kk), 'ref_trls2', R.trls2(kk));
    end
    dsections{gi} = struct('label', dlabels{gi}, 'group', g{gi}, ...
        'eventA', ev{1}, 'eventB', ev{2}, ...
        'direction', G107ExportFixtures.directionText(ev{1}, ev{2}), ...
        'empty', height(D) == 0, 'rows', {rows});
end
e.diff_sections = dsections;
end


function b = localButton(figs, labelfrag)
b = [];
btns = findall(figs, 'Style', 'pushbutton');
for h = btns(:)'
    s = get(h, 'String');
    if iscell(s)
        s = strjoin(s, ' ');
    end
    if contains(string(s), labelfrag)
        b = h;
        return;
    end
end
error('G107:noButton', 'No "%s" button on the window.', labelfrag);
end


function localRestore(shadow, oldVisible)
G107ExportFixtures.removeShadows(shadow);
set(groot, 'DefaultFigureVisible', oldVisible);
end
