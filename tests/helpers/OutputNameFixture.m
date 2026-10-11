classdef OutputNameFixture
    %OUTPUTNAMEFIXTURE Test-only harness for TestOutputNamePreflight (row 31,
    %LA-2). It measures what the filesystem actually accepts, drives the real
    %Specify Inputs -> Analyze callback and psyrat_run with only the
    %estimation core shadowed, and records everything a drive left behind.
    %
    %NO FIT. psyrat_computevarcomp (the engine core) and psyrat_checkconv are
    %replaced by stubs, and psyrat_detect_cmdstan_version by a stub that never
    %runs stanc. install() and every drive prove each resolution with which()
    %BEFORE anything can reach them. The REAL psyrat_computevarcompwarp runs,
    %so its preflight, its save and its sidecar block are the code under test.
    %
    %NOTHING OUTSIDE THE TEST'S OWN FOLDERS. Every entry this harness writes
    %or removes lives in a folder the test created with tempname + mkdir, and
    %the error-report writer's tempdir fallback is redirected into one of
    %them. The ground truth for "can this name be written here" is the
    %pipeline's own writer forms (save, fopen 'w') tried on the exact name in
    %such a folder on the same volume. That never calls toolbox code, so it is
    %independent of whatever check the production code uses.

    properties (Constant)
        %search bound for the component-name limit, in characters of the stem
        Kmax = 1100;
        %dialogs a drive can leave behind; removed before and after each drive
        DialogNames = {'Specify Inputs', 'Processing Complete', ...
            'All Measurements Failed', 'File Error', 'Whitespace detected', ...
            'Preflight warnings'};
        %functions that must resolve to the test's shadows before any drive
        Shadowed = {'uigetfile', 'uiputfile', 'tempdir', 'uiwait', ...
            'psyrat_startview', 'psyrat_detect_cmdstan_version', ...
            'psyrat_computevarcomp', 'psyrat_checkconv'};
        %functions that must be the checkout's own copies
        Real = {'psyrat_computevarcompwarp', 'psyrat_startproc', ...
            'psyrat_preflight_validate', 'psyrat_run', 'psyrat_write_runconfig', ...
            'psyrat_fopen_unused'};
    end

    methods (Static)

        % ---- measured filesystem facts (task-owned folder only) ------------

        function lim = measureLimits(folder)
            %Largest accepted component, in characters, for an output written
            %by save (<k>.psyrat) and for a sidecar written by fopen 'w'
            %(<k>_runconfig.json), found by bisection in FOLDER, a fresh
            %folder the test created. Each trial writes one file and deletes it.
            lim = struct();
            [lim.save, lim.saveFound, lim.saveNote] = ...
                OutputNameFixture.bisect(folder, '.psyrat', 'save');
            [lim.sidecar, lim.sidecarFound, lim.sidecarNote] = ...
                OutputNameFixture.bisect(folder, '_runconfig.json', 'fopen');
        end

        function [L, found, note] = bisect(folder, suffix, writer)
            kmax = OutputNameFixture.Kmax;
            ok = @(k) OutputNameFixture.writable(folder, [repmat('a', 1, k) suffix], writer);
            found = false;
            note = '';
            if ~ok(1)
                error('OutputNameFixture:environment', ...
                    'Even a one-character %s name cannot be written in %s.', writer, folder);
            end
            if ok(kmax)
                L = kmax + numel(suffix);
                note = sprintf('no limit up to %d characters', L);
                return;
            end
            lo = 1;
            hi = kmax;
            while hi - lo > 1
                mid = floor((lo + hi) / 2);
                if ok(mid); lo = mid; else; hi = mid; end
            end
            %a boundary is quoted only when its neighbours agree with it
            if ~(ok(lo) && ~ok(hi) && ok(max(lo - 1, 1)) && ~ok(min(hi + 1, kmax)))
                error('OutputNameFixture:environment', ...
                    'The %s writer is not monotone in name length near %d characters in %s.', ...
                    writer, lo + numel(suffix), folder);
            end
            L = lo + numel(suffix);
            found = true;
        end

        function [tf, msg] = writable(folder, name, writer)
            %Write NAME in FOLDER with one of the pipeline's writer forms, then
            %remove it. FOLDER is task-owned and holds nothing else.
            p = fullfile(folder, name);
            tf = false;
            msg = '';
            switch writer
                case 'save'    %the warp's form: save(fullfile(savepath,savename),'psyrat_data')
                    psyrat_data = struct('probe', 1); %#ok<NASGU>
                    try
                        save(p, 'psyrat_data');
                        tf = isfile(p);
                    catch saveME
                        msg = saveME.message;
                    end
                case 'fopen'   %the sidecar writer's form: fopen(sidecarfile,'w')
                    [fid, msg] = fopen(p, 'w');
                    if fid >= 3
                        fclose(fid);
                        tf = isfile(p);
                    end
            end
            if isfile(p); delete(p); end
        end

        function [ok, detail] = oracle(folder, savename)
            %Whether a run saving SAVENAME could write every file it derives
            %from it: the output (save; like save, a name without an extension
            %gets .mat), the <stem>_runconfig.json sidecar and the
            %<stem>_error.txt report name (fopen). Asked of the filesystem in
            %FOLDER, a task-owned folder on the save folder's volume.
            [~, stem, ext] = fileparts(savename);
            outname = savename;
            if isempty(ext); outname = [savename '.mat']; end
            [o, mo] = OutputNameFixture.writable(folder, outname, 'save');
            [s, ms] = OutputNameFixture.writable(folder, [stem '_runconfig.json'], 'fopen');
            [r, mr] = OutputNameFixture.writable(folder, [stem '_error.txt'], 'fopen');
            ok = o && s && r;
            detail = struct('output', o, 'sidecar', s, 'report', r, ...
                'messages', {{mo, ms, mr}});
        end

        function [k, found] = largestValid(folder, makeName, lo, hi)
            %Largest K in LO..HI for which the oracle finds every file of the
            %request MAKENAME(K) writable in FOLDER, by bisection. No counting
            %rule (bytes, UTF-16 units or code points) is assumed: only the
            %writers answer. FOUND is false when even LO is not writable
            %(K = LO - 1) or HI still is (K = HI). As in bisect, a boundary is
            %returned only when its neighbours agree with it.
            ok = @(n) OutputNameFixture.oracle(folder, makeName(n));
            found = false;
            if ~ok(lo); k = lo - 1; return; end
            if ok(hi); k = hi; return; end
            a = lo;
            b = hi;
            while b - a > 1
                mid = floor((a + b) / 2);
                if ok(mid); a = mid; else; b = mid; end
            end
            if ~(ok(max(a - 1, lo)) && ~ok(min(b + 1, hi)))
                error('OutputNameFixture:environment', ...
                    'The writers are not monotone in request length near %d in %s.', a, folder);
            end
            k = a;
            found = true;
        end

        function m = member(label, requested, valid)
            %one request: the output name the run uses, its stem, and whether
            %the filesystem can hold every name derived from it
            [~, stem] = fileparts(requested);
            m = struct('label', label, 'requested', requested, 'stem', stem, ...
                'valid', logical(valid));
        end

        function h = header(n, last)
            %a valid column header of exactly n characters ending in LAST
            h = ['m' repmat('a', 1, n - 2) last];
        end

        % ---- shadows ---------------------------------------------------------

        function ctx = install(testCase, root)
            %Folders under ROOT (task-owned and canonical), the shadows, and the
            %proof that they, and not the real functions, will answer.
            ctx = struct();
            ctx.root = root;
            ctx.shadow = fullfile(root, 'shadow');
            ctx.input = fullfile(root, 'input');
            ctx.save = fullfile(root, 'save');
            ctx.tmpfallback = fullfile(root, 'tmpfallback');
            for f = {ctx.shadow, ctx.input, ctx.save, ctx.tmpfallback}
                mkdir(f{1});
            end
            OutputNameFixture.writeShadows(ctx.shadow);
            global PSYRAT_ONP_TEMPDIR %#ok<GVMIS>
            PSYRAT_ONP_TEMPDIR = [ctx.tmpfallback filesep];
            OutputNameFixture.watch({});
            OutputNameFixture.resetCounters();
            addpath(ctx.shadow, '-begin');
            rehash;
            testCase.addTeardown(@() OutputNameFixture.uninstall(ctx));
            OutputNameFixture.assertResolution(testCase, ctx);
        end

        function assertResolution(testCase, ctx)
            for n = OutputNameFixture.Shadowed
                testCase.assertEqual(which(n{1}), fullfile(ctx.shadow, [n{1} '.m']), ...
                    sprintf('the %s shadow is not first on the path; no drive may run', n{1}));
            end
            root = testCase.projectRoot();
            for n = OutputNameFixture.Real
                testCase.assertTrue(startsWith(which(n{1}), root), ...
                    sprintf('%s must be the checkout''s own copy (got %s)', n{1}, which(n{1})));
            end
        end

        function uninstall(ctx)
            names = {'PSYRAT_ONP_UIGETFILE', 'PSYRAT_ONP_UIPUTFILE', 'PSYRAT_ONP_TEMPDIR', ...
                'PSYRAT_ONP_ENGINE', 'PSYRAT_ONP_DETECT', 'PSYRAT_ONP_UIWAIT', ...
                'PSYRAT_ONP_STARTVIEW', 'PSYRAT_ONP_WATCH'};
            clear('global', names{:});
            %the base-class teardown may already have restored the path
            if any(strcmp(strsplit(path, pathsep), ctx.shadow))
                rmpath(ctx.shadow);
            end
            rehash;
            OutputNameFixture.closeDialogs();
        end

        function writeShadows(folder)
            OutputNameFixture.writeLines(fullfile(folder, 'uigetfile.m'), {
                'function [filename,pathname,filterindex] = uigetfile(varargin)'
                '% TEST SHADOW (TestOutputNamePreflight): answers the open dialog from a global.'
                'global PSYRAT_ONP_UIGETFILE'
                '[p,n,e] = fileparts(PSYRAT_ONP_UIGETFILE);'
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'uiputfile.m'), {
                'function [filename,pathname,filterindex] = uiputfile(varargin)'
                '% TEST SHADOW (TestOutputNamePreflight): answers the save dialog from a global; 0 is Cancel.'
                'global PSYRAT_ONP_UIPUTFILE'
                'if isempty(PSYRAT_ONP_UIPUTFILE) || isequal(PSYRAT_ONP_UIPUTFILE, 0)'
                '    filename = 0; pathname = 0; filterindex = 0; return;'
                'end'
                '[p,n,e] = fileparts(PSYRAT_ONP_UIPUTFILE);'
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'tempdir.m'), {
                'function d = tempdir'
                '% TEST SHADOW (TestOutputNamePreflight): the report writer''s last fallback, inside the test''s folder.'
                'global PSYRAT_ONP_TEMPDIR'
                'if ~isempty(PSYRAT_ONP_TEMPDIR)'
                '    d = PSYRAT_ONP_TEMPDIR; return;'
                'end'
                'here = fileparts(mfilename(''fullpath''));'
                'rmpath(here);'
                'restore = onCleanup(@() addpath(here, ''-begin'')); %#ok<NASGU>'
                'd = tempdir;'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'uiwait.m'), {
                'function uiwait(h, varargin)'
                '% TEST SHADOW (TestOutputNamePreflight): records the dialog the user acknowledges, then dismisses it.'
                'global PSYRAT_ONP_UIWAIT'
                'if nargin >= 1 && ~isempty(h) && all(ishghandle(h))'
                '    PSYRAT_ONP_UIWAIT{end + 1} = get(h(1), ''Name'');'
                '    delete(h);'
                'end'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'psyrat_startview.m'), {
                'function psyrat_startview(varargin)'
                '% TEST SHADOW (TestOutputNamePreflight): a successful single run would open the viewer; count it.'
                'global PSYRAT_ONP_STARTVIEW'
                'PSYRAT_ONP_STARTVIEW = PSYRAT_ONP_STARTVIEW + 1;'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'psyrat_detect_cmdstan_version.m'), {
                'function [verstr, home, ranok] = psyrat_detect_cmdstan_version(varargin)'
                '% TEST SHADOW (TestOutputNamePreflight): never runs stanc; reports CmdStan absent and counts the call.'
                'global PSYRAT_ONP_DETECT'
                'PSYRAT_ONP_DETECT = PSYRAT_ONP_DETECT + 1;'
                'verstr = ''''; home = ''''; ranok = false;'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'psyrat_computevarcomp.m'), {
                'function REL = psyrat_computevarcomp(varargin) %#ok<INUSD>'
                '% TEST SHADOW (TestOutputNamePreflight): no model is built, compiled or sampled, and'
                '% CmdStan is never called. Records the output the calling warp is estimating for,'
                '% and the state of every watched pre-existing entry at that moment, then returns'
                '% the minimal converged result the REAL warp carries on to its own save.'
                'global PSYRAT_ONP_ENGINE PSYRAT_ONP_WATCH'
                'try'
                '    nm = evalin(''caller'', ''psyrat_data.proc.savename'');'
                'catch'
                '    nm = ''<the caller holds no psyrat_data.proc.savename>'';'
                'end'
                'PSYRAT_ONP_ENGINE{end + 1} = char(nm);'
                'for w = 1:numel(PSYRAT_ONP_WATCH)'
                '    PSYRAT_ONP_WATCH(w).atengine{end + 1} = OutputNameFixture.entryState(PSYRAT_ONP_WATCH(w).path);'
                'end'
                'REL = struct(''analysis'', ''stub'', ''out'', struct(''conv'', struct(''converged'', 1)));'
                'end'});
            OutputNameFixture.writeLines(fullfile(folder, 'psyrat_checkconv.m'), {
                'function REL = psyrat_checkconv(REL, varargin) %#ok<INUSD>'
                '% TEST SHADOW (TestOutputNamePreflight): marks the stub result converged so the real warp saves.'
                'REL.out.conv.converged = 1;'
                'end'});
        end

        function writeLines(file, lines)
            fid = fopen(file, 'w');
            assert(fid >= 3, 'OutputNameFixture:shadow', 'cannot write %s', file);
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

        function resetCounters()
            global PSYRAT_ONP_ENGINE PSYRAT_ONP_DETECT PSYRAT_ONP_UIWAIT PSYRAT_ONP_STARTVIEW %#ok<GVMIS>
            PSYRAT_ONP_ENGINE = {};
            PSYRAT_ONP_DETECT = 0;
            PSYRAT_ONP_UIWAIT = {};
            PSYRAT_ONP_STARTVIEW = 0;
        end

        function watch(paths)
            %entries whose state the engine stub records at estimation entry
            global PSYRAT_ONP_WATCH %#ok<GVMIS>
            PSYRAT_ONP_WATCH = struct('path', {}, 'before', {}, 'atengine', {});
            for i = 1:numel(paths)
                PSYRAT_ONP_WATCH(end + 1) = struct('path', paths{i}, ...
                    'before', OutputNameFixture.entryState(paths{i}), 'atengine', {{}});
            end
        end

        function W = watched()
            global PSYRAT_ONP_WATCH %#ok<GVMIS>
            W = PSYRAT_ONP_WATCH;
        end

        function s = entryState(p)
            %existence, kind, bytes (files) and child names (folders) of one path
            s = struct('exists', false, 'isdir', false, ...
                'bytes', zeros(1, 0, 'uint8'), 'children', {{}});
            if isfolder(p)
                s.exists = true;
                s.isdir = true;
                d = dir(p);
                s.children = sort(setdiff({d.name}, {'.', '..'}));
            elseif isfile(p)
                s.exists = true;
                fid = fopen(p, 'r');
                s.bytes = reshape(fread(fid, Inf, '*uint8'), 1, []);
                fclose(fid);
            end
        end

        % ---- fixtures and drives --------------------------------------------

        function csv = writeFixture(ctx, headers, nancol)
            %The balanced single-trial file of the Analyze drive tests: 6 ids x
            %2 groups x 2 events x 2 occasions x 3 trials, deterministic
            %scores. NANCOL (0 for none; a vector for several) writes NA in
            %those measurement columns at data row 3, so only those members
            %fail the preflight.
            csv = fullfile(ctx.input, 'drive.csv');
            fid = fopen(csv, 'w');
            fprintf(fid, 'id,%s,group,event,time\n', strjoin(headers, ','));
            row = 0;
            for g = 1:2
                for s = 1:6
                    for e = 1:2
                        for t = 1:2
                            for r = 1:3
                                row = row + 1;
                                vals = cell(1, numel(headers));
                                for m = 1:numel(headers)
                                    vals{m} = sprintf('%0.3f', mod((0.29 * m + 0.37) * row, 5) - 2.5);
                                    if any(m == nancol) && row == 3; vals{m} = 'NA'; end
                                end
                                fprintf(fid, 'g%ds%d,%s,%s,e%d,t%d\n', g, s, strjoin(vals, ','), ...
                                    char('A' + g - 1), e, t);
                            end
                        end
                    end
                end
            end
            fclose(fid);
        end

        function prefs = guiPrefs(testCase, csv, headers, nancol)
            %Preferences as the Analyze drive tests build them, group and
            %event selected. The weld proves every header survived the import
            %unchanged: a header readtable altered would build another name.
            raw = psyrat_readtable('file', csv);
            collist = raw.Properties.VariableNames;
            for m = 1:numel(headers)
                testCase.assertTrue(ismember(headers{m}, collist), sprintf( ...
                    ['fixture: readtable must keep the %d-character header unchanged ' ...
                    '(namelengthmax is %d on this release)'], numel(headers{m}), namelengthmax));
                col = raw.(headers{m});
                testCase.assertTrue(isnumeric(col), 'fixture: every measurement column must import as numeric');
                if any(m == nancol)
                    testCase.assertTrue(isnan(col(3)), 'fixture: the NA token must import as NaN');
                end
            end
            allcols = [collist {'(none)'}];
            nonecol = numel(allcols);
            colof = @(n) find(strcmp(allcols, n), 1);
            prefs = psyrat_defaults;
            prefs.ver = psyrat_defineversion;
            prefs.guis.fsize = psyrat_guifsize;
            inp = struct();
            inp.id = colof('id');
            inp.meas = cellfun(@(n) colof(n), headers);
            inp.group = colof('group');
            inp.event = colof('event');
            inp.time = nonecol;
            inp.dim1 = nonecol;
            inp.dim2 = nonecol;
            inp.weight = nonecol;
            inp.whichgroups = '';
            inp.whichgroupscol = [];
            inp.whichevents = '';
            inp.whicheventscol = [];
            inp.whichtimes = '';
            inp.whichtimescol = [];
            prefs.proc.inp = inp;
        end

        function obs = driveGui(testCase, ctx, csv, prefs, typedname)
            %One Analyze through the real Specify Inputs screen; the save
            %dialog answers TYPEDNAME in the save folder.
            global PSYRAT_ONP_UIGETFILE PSYRAT_ONP_UIPUTFILE %#ok<GVMIS>
            OutputNameFixture.assertResolution(testCase, ctx);
            OutputNameFixture.resetCounters();
            PSYRAT_ONP_UIGETFILE = csv;
            PSYRAT_ONP_UIPUTFILE = fullfile(ctx.save, typedname);
            OutputNameFixture.closeDialogs();
            close all force;
            figure('Tag', 'psyrat_gui', 'Visible', 'off');   %the stand-in home screen
            psyrat_startproc('psyrat_prefs', prefs);
            fig = findall(groot, 'Type', 'figure', 'Name', 'Specify Inputs');
            testCase.assertNumElements(fig, 1, 'expected exactly one Specify Inputs screen');
            btn = findall(fig, 'Style', 'pushbutton', 'String', 'Analyze');
            testCase.assertNumElements(btn, 1, 'Analyze button not found');
            cb = btn.Callback;
            fcn = cb{1}; %#ok<NASGU>
            cbME = [];
            console = evalc('try, fcn(btn, [], cb{2:end}); catch cbME, end');
            obs = OutputNameFixture.collect(console, cbME);
        end

        function obs = driveScripted(testCase, ctx, csv, savename)
            %psyrat_run with SAVENAME in the save folder; its error, if any,
            %is what the caller would see.
            OutputNameFixture.assertResolution(testCase, ctx);
            OutputNameFixture.resetCounters();
            args = {'file', csv, 'idcol', 'id', 'meascol', 'meas1', 'groupcol', 'group', ...
                'eventcol', 'event', 'savepath', ctx.save, 'savename', savename}; %#ok<NASGU>
            runME = [];
            console = evalc('try, psyrat_run(args{:}); catch runME, end');
            obs = OutputNameFixture.collect(console, runME);
        end

        function obs = collect(console, escaped)
            global PSYRAT_ONP_ENGINE PSYRAT_ONP_DETECT PSYRAT_ONP_UIWAIT PSYRAT_ONP_STARTVIEW %#ok<GVMIS>
            obs = struct();
            obs.console = console;
            obs.escaped = escaped;
            obs.engine = reshape(PSYRAT_ONP_ENGINE, 1, []);
            obs.detect = PSYRAT_ONP_DETECT;
            obs.uiwait = PSYRAT_ONP_UIWAIT;
            obs.startview = PSYRAT_ONP_STARTVIEW;
            obs.dialogs = struct( ...
                'complete', OutputNameFixture.dialogText('Processing Complete'), ...
                'allfailed', OutputNameFixture.dialogText('All Measurements Failed'), ...
                'fileerror', OutputNameFixture.dialogText('File Error'));
        end

        function closeDialogs()
            for n = OutputNameFixture.DialogNames
                delete(findall(groot, 'Type', 'figure', 'Name', n{1}));
            end
            delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui'));
        end

        function txt = dialogText(name)
            %a dialog's text as one whitespace-normalized line ('' when absent);
            %message boxes wrap long lines, so the pieces are joined with spaces
            txt = '';
            fig = findall(groot, 'Type', 'figure', 'Name', name);
            if isempty(fig); return; end
            t = findall(fig, 'Type', 'text');
            parts = {};
            for i = 1:numel(t)
                parts = [parts; cellstr(char(t(i).String))]; %#ok<AGROW>
            end
            txt = regexprep(strjoin(parts', ' '), '\s+', ' ');
        end

        % ---- observation ----------------------------------------------------

        function S = snapshot(folder)
            %every entry under FOLDER, recursively: full path, folder flag, bytes
            d = dir(fullfile(folder, '**', '*'));
            d = d(~ismember({d.name}, {'.', '..'}));
            S = struct('path', {}, 'isdir', {}, 'bytes', {});
            for i = 1:numel(d)
                p = fullfile(d(i).folder, d(i).name);
                b = zeros(1, 0, 'uint8');
                if ~d(i).isdir
                    fid = fopen(p, 'r');
                    b = reshape(fread(fid, Inf, '*uint8'), 1, []);
                    fclose(fid);
                end
                S(end + 1) = struct('path', p, 'isdir', d(i).isdir, 'bytes', b); %#ok<AGROW>
            end
        end

        function R = newReports(before, after)
            %per-measurement error reports the drive added, wherever they went
            R = struct('path', {}, 'name', {}, 'requested', {}, 'identifier', {}, 'message', {});
            added = after(~ismember({after.path}, {before.path}));
            for i = 1:numel(added)
                a = added(i);
                if a.isdir || ~endsWith(a.path, '.txt'); continue; end
                txt = native2unicode(a.bytes, 'UTF-8');
                if ~startsWith(txt, 'PsyRAT measurement processing error'); continue; end
                [~, n, e] = fileparts(a.path);
                R(end + 1) = struct('path', a.path, 'name', [n e], ... %#ok<AGROW>
                    'requested', OutputNameFixture.token(txt, 'Requested output file: ([^\n]*)'), ...
                    'identifier', OutputNameFixture.token(txt, 'Error identifier: ([^\n]*)'), ...
                    'message', OutputNameFixture.token(txt, 'Error message:\n(.*?)(?:\n\nStack trace:|$)'));
            end
        end

        function s = token(txt, pattern)
            t = regexp(txt, pattern, 'tokens', 'once');
            if isempty(t); s = ''; else; s = strtrim(t{1}); end
        end

        function verifyOutcome(testCase, ctx, obs, members, before, after, route, mayChange)
            %The contract for one drive. Writable requests run and save under
            %their requested names; a request the folder cannot hold is refused
            %BEFORE estimation, visibly, with the requested name and a remedy;
            %nothing else appears, disappears or changes. MAYCHANGE lists
            %pre-existing files the run may legitimately rewrite (a consented
            %overwrite target and its existing sidecar).
            if nargin < 8; mayChange = {}; end
            rowsort = @(c) reshape(sort(c), 1, []);
            valid = members([members.valid]);
            invalid = members(~[members.valid]);

            %1. estimation reached exactly the writable requests, once each
            testCase.verifyEqual(rowsort(obs.engine), rowsort({valid.requested}), ...
                'estimation must start for every writable request and for no other');
            for m = invalid(:)'
                testCase.verifyFalse(any(strcmp(obs.engine, m.requested)), sprintf( ...
                    '%s: estimation started for an output the save folder cannot hold', m.label));
            end

            %2. writable requests saved under their requested names
            for m = valid(:)'
                testCase.verifyTrue(isfile(fullfile(ctx.save, m.requested)), sprintf( ...
                    '%s: the output must be saved under its requested name', m.label));
                testCase.verifyTrue(isfile(fullfile(ctx.save, [m.stem '_runconfig.json'])), sprintf( ...
                    '%s: its sidecar must be written', m.label));
            end

            %3. no surplus entry (probe leftovers, truncated or renamed
            %outputs), nothing removed, no pre-existing file rewritten
            reports = OutputNameFixture.newReports(before, after);
            allowed = {reports.path};
            for m = valid(:)'
                allowed = [allowed, {fullfile(ctx.save, m.requested), ...
                    fullfile(ctx.save, [m.stem '_runconfig.json'])}]; %#ok<AGROW>
            end
            added = setdiff({after.path}, {before.path});
            surplus = setdiff(added, allowed);
            testCase.verifyEmpty(surplus, sprintf( ...
                'the drive left entries nobody asked for: %s', strjoin(surplus, ' | ')));
            removed = setdiff({before.path}, {after.path});
            testCase.verifyEmpty(removed, sprintf( ...
                'pre-existing entries disappeared: %s', strjoin(removed, ' | ')));
            for b = before(:)'
                if b.isdir || any(strcmp(b.path, mayChange)); continue; end
                k = find(strcmp({after.path}, b.path), 1);
                if isempty(k); continue; end
                testCase.verifyEqual(after(k).bytes, b.bytes, sprintf( ...
                    'a pre-existing file changed: %s', b.path));
            end

            %4. the refusal is visible, names the request and says how to fix it
            switch route
                case 'gui'
                    if ~isempty(obs.escaped)
                        testCase.verifyFail(['an error escaped the Analyze callback: ' obs.escaped.message]);
                    end
                    for m = invalid(:)'
                        mine = reports(strcmp({reports.requested}, m.requested));
                        testCase.verifyNumElements(mine, 1, sprintf( ...
                            '%s: exactly one error report must name the refused output', m.label));
                        if isscalar(mine)
                            testCase.verifyTrue(startsWith(mine.identifier, 'psyrat_data:save'), sprintf( ...
                                '%s: the report must carry the output-name refusal (got %s)', ...
                                m.label, mine.identifier));
                            testCase.verifyTrue(contains(mine.message, 'How to fix:'), sprintf( ...
                                '%s: the refusal must say how to fix it', m.label));
                            testCase.verifyTrue(contains(mine.message, m.stem), sprintf( ...
                                '%s: the refusal must name the output it refuses', m.label));
                        end
                    end
                    for m = valid(:)'
                        testCase.verifyFalse(any(strcmp({reports.requested}, m.requested)), sprintf( ...
                            '%s: a writable request must not fail', m.label));
                    end
                case 'script'
                    if isempty(invalid)
                        if ~isempty(obs.escaped)
                            testCase.verifyFail(['a writable scripted request raised: ' obs.escaped.message]);
                        end
                    else
                        m = invalid(1);
                        %a verification, not an assertion: a missing refusal is
                        %the defect itself (the unfixed sidecar case raised
                        %nothing), so it must FAIL the test, not end it as
                        %Incomplete; the guard only skips the checks that need
                        %an error object to read
                        testCase.verifyNotEmpty(obs.escaped, sprintf( ...
                            '%s: a refused scripted request must raise to its caller', m.label));
                        if ~isempty(obs.escaped)
                            testCase.verifyTrue(startsWith(obs.escaped.identifier, 'psyrat_data:save'), sprintf( ...
                                '%s: the scripted refusal must be the output-name error (got %s: %s)', ...
                                m.label, obs.escaped.identifier, obs.escaped.message));
                            testCase.verifyTrue(contains(obs.escaped.message, 'How to fix:'), sprintf( ...
                                '%s: the refusal must say how to fix it', m.label));
                            testCase.verifyTrue(contains(obs.escaped.message, m.stem), sprintf( ...
                                '%s: the refusal must name the output it refuses', m.label));
                        end
                    end
            end
        end

        function verifyNothingElse(testCase, before, after, allowed)
            %Part 3 of verifyOutcome on its own, for drives whose members the
            %test checks itself (a data error, or several failing members): no
            %entry appears outside ALLOWED, none disappears, and no
            %pre-existing file changes.
            surplus = setdiff(setdiff({after.path}, {before.path}), allowed);
            testCase.verifyEmpty(surplus, sprintf( ...
                'the drive left entries nobody asked for: %s', strjoin(surplus, ' | ')));
            removed = setdiff({before.path}, {after.path});
            testCase.verifyEmpty(removed, sprintf( ...
                'pre-existing entries disappeared: %s', strjoin(removed, ' | ')));
            for b = before(:)'
                if b.isdir; continue; end
                k = find(strcmp({after.path}, b.path), 1);
                if isempty(k); continue; end
                testCase.verifyEqual(after(k).bytes, b.bytes, sprintf( ...
                    'a pre-existing file changed: %s', b.path));
            end
        end

    end
end
