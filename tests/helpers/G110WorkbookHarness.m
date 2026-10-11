classdef G110WorkbookHarness
    %G110WORKBOOKHARNESS Routes, dialog doubles, I/O seams and independent
    %readers for TestWorkbookReplacement (G110) and its local reopen probe.
    %
    %G110: every GUI XLSX Save callback wrote its header and table with
    %writecell/writetable straight into the file the save dialog returned. For
    %an existing workbook those writers update only the cells they cover, so
    %the old workbook's trailing rows, wider columns and extra sheets survive a
    %"Replace". The approved correction writes a fresh sibling workbook,
    %checks that the write completed, and only then replaces the consented
    %target; a writer failure or a failed final replacement leaves the original
    %bytes untouched and removes only the task's own scratch file.
    %
    %What this harness provides (no CmdStan, no fit, no desktop):
    % * ROUTES. The 16 button routes that reach the 13 census writer sites in
    %   10 files. Each is built from the shared synthetic fixtures
    %   (PsyRATTestDataFactory) and fired through the button's own Callback.
    % * DIALOG DOUBLES. uiputfile answers from a global (a path, or 0 for
    %   Cancel) and counts its calls. questdlg and uiwait never block; they are
    %   a hang guard for a headless run, and their calls are recorded.
    % * I/O SEAMS. writecell, writetable, movefile and copyfile doubles that
    %   log every call made by code outside MATLAB's own toolboxes (with the
    %   caller's file:line), optionally fail one of them, and otherwise forward
    %   to the real function by taking their own folder off the path for the
    %   call (the TestNativeEngine fitlme pattern). Calls made from inside a
    %   MATLAB toolbox are forwarded uncounted. The seams are on the path only
    %   between armSeams and disarmSeams. The two corruption modes model a
    %   SILENT writer defect: every writer call into the destination folder
    %   is forwarded and returns normally, and then, if that call's own
    %   payload holds a given required text cell (the needle), exactly that
    %   one cell of the file just written is overwritten with a marker, or
    %   (numeric mode) its first finite number becomes x + eps(x). The
    %   workbook stays valid; only an exact read-back comparison can tell. The
    %   shift mode moves every writer call into the destination folder the
    %   same distance down and right (each call's 'Range' origin), so values
    %   and their relative layout are intact but nothing starts at A1. The
    %   extra mode adds, after each writer call into the destination folder, a
    %   second sheet or one text cell outside the payload (or nothing: its
    %   control), so only a one-sheet or used-range check can tell. The
    %   damage modes of the move/copy seam change or remove the chosen file
    %   before reporting the failure, so the helper's restore is exercised;
    %   the wrongButReportsSuccess mode carries out the first move onto the
    %   chosen file, overwrites the file with other bytes and reports success,
    %   so only reading the chosen file back after the move can tell.
    % * READERS. staleReport unzips the workbook and scans its XML, which does
    %   not depend on readcell; content reads the first sheet with the
    %   timestamps masked; completeness checks a workbook against expectations
    %   taken from the callback's own table, the fixture's IDs, or the viewer's
    %   on-screen table -- never from the writer under test.
    %
    %Methods that can fail raise an error with a 'G110Harness:' identifier or
    %return a list of problems; the test class turns those into verifications,
    %so one broken route cannot hide the others.

    properties (Constant)
        %Where the planted "old workbook" puts its stale cells. Each lies
        %outside every census export (TestWorkbookReplacement asserts that on
        %a save to a new file before relying on it), so a save that really
        %replaces the workbook removes them and one that writes in place
        %leaves them.
        LongCell = 'A400'        % long-to-short: a row below every export
        LongRow = 400
        WideCell = 'CV2'         % wide-to-narrow: column 100, right of every export
        WideCol = 100
        NearCell = 'B1'          % wide-to-narrow next door: row 1 of every export is one cell wide
        OldSheet = 'G110_Old'    % a sheet no export writes
        TargetName = 'g110_target.xlsx'
        CsvTargetName = 'g110_target.csv'
        SeamKey = 'G110_io_seam' % getappdata(0, ...) key of the armed seam state
    end

    methods (Static)

        % ---------------------------------------------------------------
        % Environment
        % ---------------------------------------------------------------

        function env = install(varargin)
            %Create the folders, write the doubles and put the dialog doubles
            %first on the path. Name-value options (used by the local probe):
            %'Root' - folder that holds work/ and reference/ (default: a new
            %tempname); 'Keep' - true leaves Root in place at uninstall. The
            %doubles always live in their own temporary folders, never beside
            %the artifacts. Pair every install with uninstall (tests register
            %it with addTeardown, so it runs when a verification fails too).
            global G110_SAVE_ANSWER G110_SAVE_CALLS G110_QUESTDLG_CALLS G110_UIWAIT_CALLS %#ok<GVMIS>
            p = inputParser;
            p.addParameter('Root', '', @(x) ischar(x) || isstring(x));
            p.addParameter('Keep', false, @(x) islogical(x) || isnumeric(x));
            p.parse(varargin{:});

            env = struct();
            env.keep = logical(p.Results.Keep);
            if strlength(string(p.Results.Root)) == 0
                env.root = tempname;
            else
                env.root = char(p.Results.Root);
            end
            env.work = fullfile(env.root, 'work');        % the destination folder the seams watch
            env.ref = fullfile(env.root, 'reference');    % saves to brand-new files
            env.dialogs = tempname;
            env.seams = tempname;
            env.faults = tempname;   % fread/ferror/fwrite/delete doubles (armFaults)
            env.recorders = tempname; % errordlg call recorder (armErrordlgRecorder)
            %the whole warning state, restored at uninstall (planting turns off
            %the add-sheet warning, arming turns off the name-conflict warning)
            env.warnings = warning;

            folders = {env.root, env.work, env.ref, env.dialogs, env.seams, env.faults, ...
                env.recorders};
            for k = 1:numel(folders)
                if ~isfolder(folders{k})
                    mkdir(folders{k});
                end
            end
            G110WorkbookHarness.writeDialogDoubles(env.dialogs);
            G110WorkbookHarness.writeSeams(env.seams);
            G110WorkbookHarness.writeFaults(env.faults);
            G110WorkbookHarness.writeRecorders(env.recorders);

            G110_SAVE_ANSWER = 0;
            G110_SAVE_CALLS = 0;
            G110_QUESTDLG_CALLS = {};
            G110_UIWAIT_CALLS = 0;

            addpath(env.dialogs, '-begin');
            rehash;
            if ~startsWith(which('uiputfile'), env.dialogs)
                G110WorkbookHarness.uninstall(env);
                error('G110Harness:double', ...
                    'The uiputfile double is not the resolved uiputfile (which: %s).', ...
                    which('uiputfile'));
            end
        end

        function uninstall(env)
            %Undo install: seams and doubles off the path, their functions out
            %of memory, globals and appdata cleared, the warning state
            %restored, the temporary folders deleted (Root only when Keep is
            %false). Safe to call twice and after the base class has already
            %restored the path.
            G110WorkbookHarness.disarmSeams(env);
            if isfield(env, 'faults')
                G110WorkbookHarness.disarmFaults(env);
            end
            if isfield(env, 'recorders')
                G110WorkbookHarness.disarmErrordlgRecorder(env);
            end
            clear global G110_SAVE_ANSWER G110_SAVE_CALLS G110_QUESTDLG_CALLS G110_UIWAIT_CALLS
            G110WorkbookHarness.dropFromPath(env.dialogs);
            clear('uiputfile', 'questdlg', 'uiwait');
            rehash;
            warning(env.warnings);
            G110WorkbookHarness.removeFolder(env.dialogs);
            G110WorkbookHarness.removeFolder(env.seams);
            if isfield(env, 'faults')
                G110WorkbookHarness.removeFolder(env.faults);
            end
            if isfield(env, 'recorders')
                G110WorkbookHarness.removeFolder(env.recorders);
            end
            if ~env.keep
                G110WorkbookHarness.removeFolder(env.root);
            end
        end

        function armSeams(env, targetName, varargin)
            %Put the I/O seams first on the path and arm them for one save
            %into env.work whose chosen file is TARGETNAME. Options:
            %'Writer'  - 'none' (count only), 'after' (the FailAt-th writer
            %            call into env.work writes, then fails) or 'before'
            %            (it fails without writing);
            %'FailAt'  - which writer call into env.work fails (default 1);
            %'Replace' - 'none' or 'failFirst' (the first movefile/copyfile
            %            onto env.work/TARGETNAME fails the way the real
            %            function reports a failure; later ones proceed, so a
            %            design that restores a moved-aside original can).
            %'Writer' also takes 'corruptHeader' or 'corruptBody' with
            %'Needle' - the exact text of one required cell. Every writer
            %            call into env.work is forwarded and returns normally;
            %            when its payload holds the needle, that one cell of
            %            the file just written becomes G110_CORRUPT_HEADER or
            %            G110_CORRUPT_BODY (corruptWritten). 'corruptNumeric'
            %            needs no needle: the first finite number of a
            %            payload becomes the next larger double, x + eps(x),
            %            a change only an exact comparison can see. 'shift'
            %            forwards every writer call into env.work with its
            %            'Range' origin (A1 when it gives none) moved by
            %'Shift'   - [rows columns], default [1 0]: down and right. [0 0]
            %            rewrites each origin to itself (the shift control).
            %'Writer' also takes 'extra': every writer call into env.work is
            %            forwarded, and then the file just written gets
            %'Extra'   - 'sheet' (a sheet G110_Extra holding one cell), a cell
            %            reference such as 'C1' (G110_EXTRA in that cell of the
            %            first sheet), or 'none' (default; nothing, the
            %            control) (extraWritten).
            %'Replace' also takes 'damageFirst' (the first move/copy onto
            %            the target overwrites the target with other bytes,
            %            then reports failure), 'removeFirst' (it deletes the
            %            target, then reports failure) and 'damageAll' (as
            %            damageFirst, and every later move/copy onto the target
            %            fails too, so the original cannot be put back).
            %            'wrongButReportsSuccess' carries out the first
            %            move/copy onto the target, then overwrites the target
            %            with other bytes and reports success; later ones
            %            proceed.
            p = inputParser;
            p.addParameter('Writer', 'none', @(x) any(strcmp(x, ...
                {'none', 'after', 'before', 'corruptHeader', 'corruptBody', 'corruptNumeric', 'shift', ...
                'extra'})));
            p.addParameter('FailAt', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
            p.addParameter('Shift', [1 0], @(x) isnumeric(x) && numel(x) == 2 && ...
                all(x >= 0) && all(x == round(x)));
            p.addParameter('Replace', 'none', @(x) any(strcmp(x, ...
                {'none', 'failFirst', 'damageFirst', 'removeFirst', 'damageAll', 'wrongButReportsSuccess'})));
            p.addParameter('Needle', '', @(x) ischar(x) || isstring(x));
            p.addParameter('Extra', 'none', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            corrupting = any(strcmp(p.Results.Writer, {'corruptHeader', 'corruptBody'}));
            if corrupting && strlength(string(p.Results.Needle)) == 0
                error('G110Harness:needle', 'A corruption mode needs the text of the cell to corrupt.');
            end

            s = struct();
            s.armed = true;
            s.workDir = env.work;
            s.targetName = char(targetName);
            s.writerMode = p.Results.Writer;
            s.writerFailAt = p.Results.FailAt;
            s.writerCalls = 0;          % writer calls into env.work
            s.writerInjected = false;
            s.writerLog = cell(0, 3);   % {function, destination, caller file:line (name)}
            s.replaceMode = p.Results.Replace;
            s.replaceCalls = 0;         % movefile/copyfile calls onto env.work/targetName
            s.replaceInjected = false;
            s.moveLog = cell(0, 4);     % {function, source, destination, caller}
            s.corruptNeedle = '';
            s.corruptMarker = '';
            if corrupting
                s.corruptNeedle = char(p.Results.Needle);
                if strcmp(p.Results.Writer, 'corruptHeader')
                    s.corruptMarker = 'G110_CORRUPT_HEADER';
                else
                    s.corruptMarker = 'G110_CORRUPT_BODY';
                end
            end
            s.corruptApplied = 0;       % cells overwritten after a returned write
            s.corruptLog = cell(0, 4);  % {function, file, address, needle}
            s.shift = double(reshape(p.Results.Shift, 1, 2));   % [rows columns]
            s.shiftApplied = 0;         % writer calls whose origin was rewritten
            s.shiftLog = cell(0, 4);    % {function, file, origin, moved to}
            s.extra = char(p.Results.Extra);   % 'none', 'sheet' or a cell reference
            s.extraCalls = 0;           % writer calls the extra mode followed
            s.extraApplied = 0;         % of those, the ones after which something was added
            s.extraSheets = 0;          % sheets in the file after the last sheet was added
            s.extraLog = cell(0, 3);    % {function, file, what was added}
            setappdata(0, G110WorkbookHarness.SeamKey, s);

            %a folder that shadows the movefile/copyfile built-ins warns on
            %addpath; the warning state is restored at uninstall
            warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(env.seams, '-begin');
            rehash;
            names = {'writecell', 'writetable', 'movefile', 'copyfile'};
            for k = 1:numel(names)
                w = which(names{k});
                if ~startsWith(w, env.seams)
                    G110WorkbookHarness.disarmSeams(env);
                    error('G110Harness:seam', ...
                        'The %s seam is not the resolved %s (which: %s).', names{k}, names{k}, w);
                end
            end
        end

        function st = disarmSeams(env)
            %Take the seams off the path and return the state they recorded
            %([] when they were not armed).
            st = [];
            if isappdata(0, G110WorkbookHarness.SeamKey)
                st = getappdata(0, G110WorkbookHarness.SeamKey);
                rmappdata(0, G110WorkbookHarness.SeamKey);
            end
            G110WorkbookHarness.dropFromPath(env.seams);
            clear('writecell', 'writetable', 'movefile', 'copyfile', 'g110_seam_gate');
            rehash;
        end

        function armFaults(env, varargin)
            %Put the fread/ferror/fwrite/delete doubles first on the path.
            %They act ONLY on calls made directly from the file named by
            %'Caller' (default: the fresh-workbook helper,
            %psyrat_replace_workbook.m) on files in env.work; every other call,
            %MATLAB's own included, goes straight to the built-in. Options:
            %'CutRead'      - full path: the first fread of it returns only
            %                 half its bytes, and ferror on that file id then
            %                 reports an error (a read that failed part-way);
            %'FailWrite'    - name prefix: the first fwrite into a file whose
            %                 name starts with it writes half the data, then
            %                 raises;
            %'RefuseDelete' - name prefixes: delete of such a file is refused
            %                 the way MATLAB refuses it (a warning, no error,
            %                 the file stays).
            p = inputParser;
            p.addParameter('Caller', 'psyrat_replace_workbook.m', @(x) ischar(x) || isstring(x));
            p.addParameter('CutRead', '', @(x) ischar(x) || isstring(x));
            p.addParameter('FailWrite', '', @(x) ischar(x) || isstring(x));
            p.addParameter('RefuseDelete', {}, @(x) iscellstr(x) || ischar(x));
            p.parse(varargin{:});
            s = struct();
            s.armed = true;
            s.workDir = env.work;
            s.caller = char(p.Results.Caller);
            s.cutRead = char(p.Results.CutRead);
            s.readCut = false;          % the cut read happened
            s.errorFids = [];           % file ids whose read was cut
            s.ferrorCalls = 0;          % ferror calls by the caller
            s.failWrite = char(p.Results.FailWrite);
            s.writeFailed = false;      % the failing write happened
            s.refuseDelete = cellstr(p.Results.RefuseDelete);
            s.deleteRefused = {};       % full paths whose delete was refused
            setappdata(0, 'G110_fault_seam', s);
            warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(env.faults, '-begin');
            rehash;
            names = {'fread', 'ferror', 'fwrite', 'delete'};
            for k = 1:numel(names)
                w = which(names{k});
                if ~startsWith(w, env.faults)
                    G110WorkbookHarness.disarmFaults(env);
                    error('G110Harness:fault', ...
                        'The %s double is not the resolved %s (which: %s).', names{k}, names{k}, w);
                end
            end
            %Refresh the helper so that its next call finds these doubles. In the
            %G110 phase-9 probe (R2025b), after TestWorkbookReplacement had run
            %in the same session, the helper's own calls could reach the
            %built-ins with the doubles first on the path and named by which;
            %clearing only the helper restored delivery. Why is not established.
            %Only the helper is cleared, only when it is the aimed caller; it is
            %never on the call stack when armFaults runs.
            if strcmpi(s.caller, 'psyrat_replace_workbook.m')
                clear('psyrat_replace_workbook');
            end
        end

        function st = disarmFaults(env)
            %Take the fault doubles off the path and return what they recorded
            %([] when they were not armed).
            st = [];
            if isappdata(0, 'G110_fault_seam')
                st = getappdata(0, 'G110_fault_seam');
                rmappdata(0, 'G110_fault_seam');
            end
            G110WorkbookHarness.dropFromPath(env.faults);
            clear('fread', 'ferror', 'fwrite', 'delete', 'g110_fault_gate');
            rehash;
            %and the helper again if the doubles were aimed at it, so that its
            %later calls find the built-ins; a disarm with nothing armed (an
            %uninstall after the test's own disarm, or with no faults ever
            %armed) leaves the helper alone
            if ~isempty(st) && strcmpi(st.caller, 'psyrat_replace_workbook.m')
                clear('psyrat_replace_workbook');
            end
        end

        function armErrordlgRecorder(env)
            %Put the errordlg recorder first on the path: every errordlg call,
            %from any file, is recorded (caller, title, message) and then shown
            %by the real errordlg. Calls are what count: a modal errordlg
            %reuses an open dialog of the same title, so windows cannot be
            %counted.
            setappdata(0, 'G110_errordlg_calls', {});
            addpath(env.recorders, '-begin');
            rehash;
            w = which('errordlg');
            if ~startsWith(w, env.recorders)
                G110WorkbookHarness.disarmErrordlgRecorder(env);
                error('G110Harness:recorder', 'The errordlg recorder is not the resolved errordlg (which: %s).', w);
            end
        end

        function calls = disarmErrordlgRecorder(env)
            %Take the recorder off the path and return one struct per errordlg
            %call, in order: caller (file:line (name)), title, message.
            calls = struct('caller', {}, 'title', {}, 'message', {});
            if isappdata(0, 'G110_errordlg_calls')
                raw = getappdata(0, 'G110_errordlg_calls');
                rmappdata(0, 'G110_errordlg_calls');
                for k = 1:numel(raw)
                    args = raw{k}{2};
                    heading = '';
                    message = '';
                    if numel(args) >= 1
                        message = G110WorkbookHarness.joinText(args{1});
                    end
                    if numel(args) >= 2
                        heading = G110WorkbookHarness.textArg(args{2});
                    end
                    calls(end+1) = struct('caller', raw{k}{1}, 'title', heading, 'message', message); %#ok<AGROW>
                end
            end
            G110WorkbookHarness.dropFromPath(env.recorders);
            clear('errordlg');
            rehash;
        end

        function ids = openFileIds()
            %The file ids MATLAB has open. openedFiles (R2024a and later) is the
            %supported query; from R2025b fopen('all') raises
            %MATLAB:fopen:FopenAllStagedRemoval. Older supported releases fall
            %back to fopen('all').
            if exist('openedFiles') > 0 %#ok<EXIST>
                ids = openedFiles();
            else
                ids = fopen('all');
            end
            ids = reshape(double(ids), 1, []);
        end

        % ---------------------------------------------------------------
        % Routes (the census)
        % ---------------------------------------------------------------

        function keys = routeKeys()
            %The 16 button routes, in census order. Two routes reach the DoD
            %variance writer (the app's Export button and the variance
            %window's Save button) and each viewer writer has two buttons, so
            %16 routes cover the 13 writer sites in 10 files.
            keys = {'depoverallt_table', 'depoverallt_ids', ...
                'trt_reloverallt_table', 'trt_reloverallt_ids', ...
                'variancet_table', 'variancet_allsds', ...
                'depcutofft_table', 'trt_relcutofft_table', ...
                'dod_observed_table', 'dod_varcomp_export', 'dod_varcomp_window', ...
                'criterion_table', ...
                'dynrel_varcomp', 'dynrel_subject', ...
                'splits_coef', 'splits_varcomp'};
        end

        function spec = route(key)
            %One route: the writer site it reaches (file and callback), how its
            %window is opened, the Save button's label, and what a complete
            %export must contain. Fields:
            % file, site - the census writer site (relative path, callback name)
            % screen     - 'output' (a table builder run with gui = 1),
            %              'dodsetup', 'dynrel' or 'splits' (the app screens)
            % builder    - @(pd) call that opens an 'output' window
            % data       - @() the route's fixture
            % opener     - a button pressed first; the Save button is on the
            %              window it opens ('' when none)
            % save       - label fragment of the Save button
            % kind       - 'table' (the callback carries the table as its
            %              third Callback entry), 'ids', 'dodexport' or 'viewer'
            % first      - the export's first header line (cell A1)
            % provenance - whether the header carries the run provenance
            % show       - viewers: the button that shows the same table on
            %              screen (its row count is the expected body length)
            tables = fullfile('subroutines', 'tables');
            guis = fullfile('subroutines', 'guis');
            spec = struct('key', key, 'file', '', 'site', '', 'screen', 'output', ...
                'builder', [], 'data', [], 'opener', '', 'save', 'Save Table', ...
                'kind', 'table', 'first', 'Table Generated on', 'provenance', true, ...
                'show', '');
            switch key
                case 'depoverallt_table'
                    spec.file = fullfile(tables, 'psyrat_depoverallt.m');
                    spec.site = 'psyrat_saveoveralltable';
                    spec.builder = @(pd) psyrat_depoverallt('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.icData();
                    spec.first = 'Dependability Table Generated on';
                case 'depoverallt_ids'
                    spec.file = fullfile(tables, 'psyrat_depoverallt.m');
                    spec.site = 'psyrat_saveids';
                    spec.builder = @(pd) psyrat_depoverallt('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.icData();
                    spec.save = 'Save IDs';
                    spec.kind = 'ids';
                    spec.first = 'Data Generated on';
                    spec.provenance = false;
                case 'trt_reloverallt_table'
                    spec.file = fullfile(tables, 'psyrat_trt_reloverallt.m');
                    spec.site = 'psyrat_saveoveralltable';
                    spec.builder = @(pd) psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.trtData();
                    spec.first = 'Dependability Table Generated on';
                case 'trt_reloverallt_ids'
                    spec.file = fullfile(tables, 'psyrat_trt_reloverallt.m');
                    spec.site = 'psyrat_saveids';
                    spec.builder = @(pd) psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.trtData();
                    spec.save = 'Save IDs';
                    spec.kind = 'ids';
                    spec.first = 'Data Generated on';
                    spec.provenance = false;
                case 'variancet_table'
                    spec.file = fullfile(tables, 'psyrat_variancet.m');
                    spec.site = 'psyrat_savevartable';
                    spec.builder = @(pd) psyrat_variancet('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.icData();
                case 'variancet_allsds'
                    %the all-SDs Save lives on the window the test-retest
                    %variance table's "View All Standard Deviations" opens
                    spec.file = fullfile(tables, 'psyrat_variancet.m');
                    spec.site = 'psyrat_saveallvartable';
                    spec.builder = @(pd) psyrat_variancet('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.trtData();
                    spec.opener = 'View All';
                case 'depcutofft_table'
                    spec.file = fullfile(tables, 'psyrat_depcutofft.m');
                    spec.site = 'psyrat_saveinctrltable';
                    spec.builder = @(pd) psyrat_depcutofft('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.icData();
                case 'trt_relcutofft_table'
                    spec.file = fullfile(tables, 'psyrat_trt_relcutofft.m');
                    spec.site = 'psyrat_saveinctrltable';
                    spec.builder = @(pd) psyrat_trt_relcutofft('psyrat_data', pd, 'gui', 1);
                    spec.data = @() G110WorkbookHarness.trtData();
                case 'dod_observed_table'
                    %the app route: DoD setup screen, "Observed DoD Table",
                    %then that window's Save Table
                    spec.file = fullfile(tables, 'psyrat_dod_observedt.m');
                    spec.site = 'local_save_table';
                    spec.screen = 'dodsetup';
                    spec.data = @() G110WorkbookHarness.dodData();
                    spec.opener = 'Observed DoD';
                case 'dod_varcomp_export'
                    %the app route to the DoD variance writer: the setup
                    %screen's "Export DoD Var Components" calls
                    %psyrat_dod_varcompt with 'savefile',1,'gui',0
                    spec.file = fullfile(tables, 'psyrat_dod_varcompt.m');
                    spec.site = 'local_save_table';
                    spec.screen = 'dodsetup';
                    spec.data = @() G110WorkbookHarness.dodData();
                    spec.save = 'Export DoD';
                    spec.kind = 'dodexport';
                case 'dod_varcomp_window'
                    %the same writer through the variance window's Save Table
                    %(gui = 1; no app screen opens this window today)
                    spec.file = fullfile(tables, 'psyrat_dod_varcompt.m');
                    spec.site = 'local_save_table';
                    spec.builder = @(pd) psyrat_dod_varcompt('psyrat_data', pd, 'groupidx', 1, ...
                        'map', G110WorkbookHarness.dodMap(pd.rel), 'gui', 1);
                    spec.data = @() G110WorkbookHarness.dodData();
                case 'criterion_table'
                    spec.file = fullfile('subroutines', 'plotting', 'psyrat_criterionfigures.m');
                    spec.site = 'local_save_table';
                    spec.builder = @(pd) psyrat_criterionfigures('psyrat_prefs', ...
                        G110WorkbookHarness.criterionPrefs(), 'psyrat_data', pd, 'analysis', 'sing');
                    spec.data = @() G110WorkbookHarness.icData();
                case 'dynrel_varcomp'
                    spec.file = fullfile(guis, 'psyrat_startview_dynrel.m');
                    spec.site = 'exportDynrelTable';
                    spec.screen = 'dynrel';
                    spec.data = @() G110WorkbookHarness.withHeader( ...
                        PsyRATTestDataFactory.makeDynrelDiffSummaryData());
                    spec.save = 'Save variance table';
                    spec.kind = 'viewer';
                    spec.first = 'PsyRAT dynamic reliability: variance components';
                    spec.show = 'Show variance components';
                case 'dynrel_subject'
                    spec.file = fullfile(guis, 'psyrat_startview_dynrel.m');
                    spec.site = 'exportDynrelTable';
                    spec.screen = 'dynrel';
                    spec.data = @() G110WorkbookHarness.withHeader( ...
                        PsyRATTestDataFactory.makeDynrelSubjSummaryData());
                    spec.save = 'Save per-participant table';
                    spec.kind = 'viewer';
                    spec.first = 'PsyRAT dynamic reliability: per-participant reliability';
                    spec.show = 'Show per-participant table';
                case 'splits_coef'
                    spec.file = fullfile(guis, 'psyrat_startview_splits.m');
                    spec.site = 'exportSplitsTable';
                    spec.screen = 'splits';
                    spec.data = @() G110WorkbookHarness.withHeader( ...
                        PsyRATTestDataFactory.makeSplitsSummaryData());
                    spec.save = 'Save coefficients';
                    spec.kind = 'viewer';
                    spec.first = 'PsyRAT data-splits reliability: reliability coefficients';
                    spec.show = 'Show coefficients';
                case 'splits_varcomp'
                    spec.file = fullfile(guis, 'psyrat_startview_splits.m');
                    spec.site = 'exportSplitsTable';
                    spec.screen = 'splits';
                    spec.data = @() G110WorkbookHarness.withHeader( ...
                        PsyRATTestDataFactory.makeSplitsSummaryData());
                    spec.save = 'Save variance table';
                    spec.kind = 'viewer';
                    spec.first = 'PsyRAT data-splits reliability: variance components';
                    spec.show = 'Show variance components';
                otherwise
                    error('G110Harness:route', 'Unknown route "%s".', key);
            end
        end

        function fig = openRoute(spec, pd)
            %Open the window that carries the route's Save button. PD
            %overrides the route's fixture (the ID tests pass their own).
            if nargin < 2 || isempty(pd)
                pd = spec.data();
            end
            switch spec.screen
                case 'output'
                    delete(findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output'));
                    spec.builder(pd);
                    figs = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_output');
                    if isempty(figs)
                        error('G110Harness:open', '%s: the builder opened no output window.', spec.key);
                    end
                    fig = figs(1);
                case 'dodsetup'
                    fig = G110WorkbookHarness.openNamed(@() psyrat_startview_dodiff( ...
                        'psyrat_prefs', psyrat_defaults, 'psyrat_data', pd), ...
                        'DoD Reliability Setup');
                case 'dynrel'
                    fig = G110WorkbookHarness.openNamed(@() psyrat_startview_dynrel( ...
                        'psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd), ...
                        'Dynamic Reliability: View Results');
                case 'splits'
                    fig = G110WorkbookHarness.openNamed(@() psyrat_startview_splits( ...
                        'psyrat_prefs', psyrat_defaults(), 'psyrat_data', pd), ...
                        'Data-Splits Reliability: View Results');
                otherwise
                    error('G110Harness:open', 'Unknown screen "%s".', spec.screen);
            end
            if ~isempty(spec.opener)
                fig = G110WorkbookHarness.openChild(fig, spec.opener);
            end
        end

        function res = fireSave(fig, label, answer)
            %Press the Save button whose label contains LABEL, with the save
            %dialog double answering ANSWER (a full path, or 0 for Cancel).
            %Returns the error it raised (err, errText), the windows that
            %showed a message or opened during the press (dialogs) and their
            %text (dialogText; see dialogsSince, which also catches a reused
            %modal dialog), its console output (out), and how many times it
            %asked the save dialog (saveCalls).
            global G110_SAVE_ANSWER G110_SAVE_CALLS %#ok<GVMIS>
            G110_SAVE_ANSWER = answer;
            G110_SAVE_CALLS = 0;
            mark = G110WorkbookHarness.dialogMark();
            res = struct('err', [], 'errText', '', 'dialogs', {{}}, 'dialogText', '', ...
                'out', '', 'saveCalls', 0);
            try
                res.out = evalc('G110WorkbookHarness.pressButton(fig, label);');
            catch ME
                res.err = ME;
                res.errText = sprintf('%s: %s', ME.identifier, ME.message);
            end
            [res.dialogText, res.dialogs] = G110WorkbookHarness.dialogsSince(mark, false);
            res.saveCalls = G110_SAVE_CALLS;
        end

        function pressButton(fig, label)
            %Fire a push button's own Callback the way a click would: table
            %windows store {fcn, args...}, viewers store a function handle.
            b = G110WorkbookHarness.findButton(fig, label);
            cb = get(b, 'Callback');
            if iscell(cb)
                feval(cb{1}, b, [], cb{2:end});
            else
                feval(cb, b, []);
            end
        end

        function b = findButton(fig, label)
            %The one push button on FIG whose (joined) label contains LABEL.
            btns = findall(fig, 'Style', 'pushbutton');
            hit = gobjects(0);
            for k = 1:numel(btns)
                s = get(btns(k), 'String');
                if iscell(s)
                    s = strjoin(s, ' ');
                end
                if contains(char(s), label)
                    hit(end+1) = btns(k); %#ok<AGROW>
                end
            end
            if numel(hit) ~= 1
                error('G110Harness:button', 'Expected one "%s" button on "%s"; found %d.', ...
                    label, char(get(fig, 'Name')), numel(hit));
            end
            b = hit(1);
        end

        % ---------------------------------------------------------------
        % Planted files
        % ---------------------------------------------------------------

        function plantWorkbook(file)
            %An "old workbook" at FILE: A1 inside every export's header (so
            %any in-place write changes the file's bytes), plus the stale
            %cells and sheet listed in the Constant properties. Call only
            %while the seams are disarmed.
            state = warning('off', 'MATLAB:xlswrite:AddSheet');
            restore = onCleanup(@() warning(state)); %#ok<NASGU>
            writecell({'G110_PLANT_A1'}, file);
            writecell({'G110_STALE_NEAR'}, file, 'Range', G110WorkbookHarness.NearCell);
            writecell({'G110_STALE_WIDE'}, file, 'Range', G110WorkbookHarness.WideCell);
            writecell({'G110_STALE_LONG'}, file, 'Range', G110WorkbookHarness.LongCell);
            writecell({'G110_STALE_SHEET'}, file, 'Sheet', G110WorkbookHarness.OldSheet);
        end

        function plantCsv(file)
            %An "old CSV" longer and wider than any export.
            lines = cell(501, 1);
            lines{1} = 'G110_STALE_CSV first line, longer than any export header,,,,,,,,,,,,';
            for k = 1:500
                lines{k + 1} = sprintf('G110_STALE_CSV row %d,1,2,3,4,5,6,7,8,9,10,11,12', k);
            end
            G110WorkbookHarness.writeLines(file, lines);
        end

        function plantDecoys(folder, targetName)
            %Entries the save must leave exactly as they are. Two files: one
            %named like a scratch file beside the target and an unrelated
            %neighbor. And, named like the fresh-workbook helper's own
            %temporary files (psyrat_stage_<k>.xlsx, psyrat_backup_<k>.xlsx),
            %for k = 1..3: an empty (zero-byte) file, a FOLDER, and an empty
            %file whose name differs only in case. None is this task's file,
            %so none may be reused, written, moved or removed.
            [~, n, e] = fileparts(targetName);
            G110WorkbookHarness.writeLines(fullfile(folder, [n '.psyrat_partial_1' e]), ...
                {'G110 decoy: named like a scratch file, but this export did not create it.'});
            G110WorkbookHarness.writeLines(fullfile(folder, 'g110_neighbor.txt'), ...
                {'G110 neighbor: an unrelated file beside the chosen workbook.'});
            for prefix = {'psyrat_stage_', 'psyrat_backup_'}
                G110WorkbookHarness.emptyFile(fullfile(folder, [prefix{1} '1.xlsx']));
                mkdir(fullfile(folder, [prefix{1} '2.xlsx']));
                G110WorkbookHarness.emptyFile(fullfile(folder, upper([prefix{1} '3.xlsx'])));
            end
        end

        function emptyFile(file)
            %A zero-byte file.
            fid = fopen(file, 'w');
            if fid < 0
                error('G110Harness:write', 'Cannot write %s.', file);
            end
            fclose(fid);
        end

        function snap = snapshot(folder)
            %Every entry in FOLDER: files with their bytes, and folders (their
            %names end in a file separator, so a folder replaced by a file of
            %the same name shows as removed and added).
            L = dir(folder);
            L = L(~ismember({L.name}, {'.', '..'}));
            snap = struct('name', {}, 'bytes', {});
            for k = 1:numel(L)
                if L(k).isdir
                    snap(end+1).name = [L(k).name filesep]; %#ok<AGROW>
                    snap(end).bytes = [];
                else
                    snap(end+1).name = L(k).name; %#ok<AGROW>
                    snap(end).bytes = G110WorkbookHarness.fileBytes(fullfile(folder, L(k).name));
                end
            end
        end

        function p = compareSnapshots(before, after, targetName, targetMustMatch)
            %Problems between two snapshots: files added (left behind),
            %files removed, and files whose bytes changed. TARGETNAME is
            %compared only when TARGETMUSTMATCH is true.
            p = {};
            nb = sort({before.name});
            na = sort({after.name});
            added = setdiff(na, nb);
            removed = setdiff(nb, na);
            if ~isempty(added)
                p{end+1} = sprintf('entries left behind: %s', strjoin(added, ', '));
            end
            if ~isempty(removed)
                p{end+1} = sprintf('entries removed: %s', strjoin(removed, ', '));
            end
            for k = 1:numel(before)
                name = before(k).name;
                j = find(strcmp({after.name}, name), 1);
                if isempty(j) || (strcmp(name, targetName) && ~targetMustMatch)
                    continue;
                end
                if ~isequal(before(k).bytes, after(j).bytes)
                    p{end+1} = sprintf('%s changed (%d -> %d bytes)', name, ...
                        numel(before(k).bytes), numel(after(j).bytes)); %#ok<AGROW>
                end
            end
        end

        % ---------------------------------------------------------------
        % Seam records
        % ---------------------------------------------------------------

        function w = classifyWrites(st, workDir, targetName)
            %Sort the logged writer calls: into the chosen file itself
            %(inplace), into another file in its folder (sibling), or
            %anywhere else (elsewhere). Each entry names the call and the
            %caller's file:line.
            w = struct('inplace', {{}}, 'sibling', {{}}, 'elsewhere', {{}});
            for k = 1:size(st.writerLog, 1)
                [folder, n, e] = fileparts(st.writerLog{k, 2});
                what = sprintf('%s from %s into %s', st.writerLog{k, 1}, ...
                    st.writerLog{k, 3}, st.writerLog{k, 2});
                if ~G110WorkbookHarness.sameFolder(folder, workDir)
                    w.elsewhere{end+1} = what;
                elseif strcmp([n e], targetName)
                    w.inplace{end+1} = what;
                else
                    w.sibling{end+1} = what;
                end
            end
        end

        function lines = describeWrites(st)
            %The logged writer and move calls, one line each.
            lines = {};
            for k = 1:size(st.writerLog, 1)
                lines{end+1} = sprintf('%s from %s into %s', st.writerLog{k, 1}, ...
                    st.writerLog{k, 3}, st.writerLog{k, 2}); %#ok<AGROW>
            end
            for k = 1:size(st.moveLog, 1)
                lines{end+1} = sprintf('%s from %s: %s -> %s', st.moveLog{k, 1}, ...
                    st.moveLog{k, 4}, st.moveLog{k, 2}, st.moveLog{k, 3}); %#ok<AGROW>
            end
            if isfield(st, 'corruptLog')
                for k = 1:size(st.corruptLog, 1)
                    lines{end+1} = sprintf('silently corrupted %s of %s after %s (needle "%s")', ...
                        st.corruptLog{k, 3}, st.corruptLog{k, 2}, st.corruptLog{k, 1}, ...
                        st.corruptLog{k, 4}); %#ok<AGROW>
                end
            end
            if isfield(st, 'shiftLog')
                for k = 1:size(st.shiftLog, 1)
                    lines{end+1} = sprintf('%s into %s moved from %s to %s', st.shiftLog{k, 1}, ...
                        st.shiftLog{k, 2}, st.shiftLog{k, 3}, st.shiftLog{k, 4}); %#ok<AGROW>
                end
            end
            if isfield(st, 'extraLog')
                for k = 1:size(st.extraLog, 1)
                    lines{end+1} = sprintf('after %s into %s, added %s', st.extraLog{k, 1}, ...
                        st.extraLog{k, 2}, st.extraLog{k, 3}); %#ok<AGROW>
                end
            end
            if isempty(lines)
                lines = {'(none)'};
            end
        end

        % ---------------------------------------------------------------
        % Readers
        % ---------------------------------------------------------------

        function rep = staleReport(file)
            %Scan the workbook's own XML (unzipped, independent of readcell)
            %for the planted markers, and list its sheets two ways: from
            %xl/workbook.xml and from MATLAB's sheetnames.
            xml = G110WorkbookHarness.workbookXml(file);
            rep = struct();
            rep.long = contains(xml.all, 'G110_STALE_LONG');
            rep.wide = contains(xml.all, 'G110_STALE_WIDE');
            rep.near = contains(xml.all, 'G110_STALE_NEAR');
            rep.sheet = contains(xml.all, 'G110_STALE_SHEET');
            rep.plant = contains(xml.all, 'G110_PLANT_A1');
            rep.corrupt = contains(xml.all, 'G110_CORRUPT_');   % a silent-corruption marker
            tok = regexp(xml.workbook, '<(?:\w+:)?sheet\s[^>]*?name="([^"]*)"', 'tokens');
            rep.sheets = cell(1, numel(tok));
            for k = 1:numel(tok)
                rep.sheets{k} = tok{k}{1};
            end
            rep.sheetnames = reshape(cellstr(sheetnames(file)), 1, []);
        end

        function C = content(file)
            %The first sheet as readcell returns it, every cell rendered to
            %text: empty/missing -> '', numbers -> %.15g, date-time text
            %(the header's generation time) -> <TIME>.
            C = G110WorkbookHarness.canon(readcell(file));
        end

        function txt = csvText(file)
            %A CSV export's text with the generation time masked.
            txt = G110WorkbookHarness.maskTime(fileread(file));
        end

        function p = completeness(spec, fig, file)
            %Problems with FILE as a complete export of SPEC's route. The
            %expectations never come from the writer under test: the first
            %header line is the literal the callback writes, the body comes
            %from the callback's own table (table), the fixture's IDs (ids),
            %the DoD table builder with the screen's selections (dodexport),
            %or the row count of the same table shown on screen (viewer).
            p = {};
            C = G110WorkbookHarness.content(file);
            if isempty(C)
                p{end+1} = 'the workbook reads back empty';
                return;
            end
            if ~strcmp(C{1, 1}, spec.first)
                p{end+1} = sprintf('A1 reads "%s"; the export''s first header line is "%s"', ...
                    C{1, 1}, spec.first);
            end
            if spec.provenance && ~any(strcmp(C(:, 1), 'Seed: 12345'))
                p{end+1} = 'the provenance line "Seed: 12345" is missing from column A';
            end
            r = G110WorkbookHarness.columnNameRow(C);
            if isempty(r)
                p{end+1} = 'no column-name row (the first row with two or more filled cells)';
                return;
            end
            body = C(r + 1:end, :);
            switch spec.kind
                case 'ids'
                    cb = get(G110WorkbookHarness.findButton(fig, spec.save), 'Callback');
                    expected = G110WorkbookHarness.expectedIdBlock(cb{2});
                    if ~isequal(C(r, :), {'Good IDs', 'Bad IDs', 'Good Group', 'Bad Group'})
                        p{end+1} = sprintf('the ID header row reads "%s"', strjoin(C(r, :), '|'));
                    end
                    if ~isequal(body, expected)
                        p{end+1} = sprintf(['the ID block is %dx%d; the fixture''s groups and IDs ', ...
                            'give %dx%d (or a cell differs)'], size(body, 1), size(body, 2), ...
                            size(expected, 1), size(expected, 2));
                    end
                case {'table', 'dodexport'}
                    T = G110WorkbookHarness.expectedTable(spec, fig);
                    filled = nnz(~cellfun(@(x) isempty(strtrim(x)), C(r, :)));
                    if filled ~= width(T) || size(C, 2) ~= width(T)
                        p{end+1} = sprintf(['the column-name row has %d filled cells in a sheet %d ', ...
                            'columns wide; the table has %d columns'], filled, size(C, 2), width(T));
                    end
                    if size(body, 1) ~= height(T)
                        p{end+1} = sprintf('%d body rows; the table has %d', size(body, 1), height(T));
                    elseif ~isequal(body(:, 1), G110WorkbookHarness.canon(G110WorkbookHarness.column1(T)))
                        p{end+1} = 'the first body column differs from the table''s first column';
                    end
                case 'viewer'
                    n = G110WorkbookHarness.shownRows(fig, spec.show);
                    if size(body, 1) ~= n
                        p{end+1} = sprintf('%d body rows; the same table on screen has %d', ...
                            size(body, 1), n);
                    end
            end
        end

        function T = expectedTable(spec, fig)
            %The table a 'table' or 'dodexport' route writes, taken from the
            %button's own Callback (not from the writer).
            cb = get(G110WorkbookHarness.findButton(fig, spec.save), 'Callback');
            switch spec.kind
                case 'table'
                    %{fcn, psyrat_data, table, ...} on every table window
                    T = cb{3};
                case 'dodexport'
                    %{fcn, 'psyrat_data', pd, 'inputs', inputs, 'eventnames', names}:
                    %rebuild the table with savefile 0 from the screen's
                    %current group and ERP_1..ERP_4 selections, the way the
                    %callback's local_read_main_inputs reads them
                    args = cb(2:end);
                    pd = G110WorkbookHarness.namedArg(args, 'psyrat_data');
                    inputs = G110WorkbookHarness.namedArg(args, 'inputs');
                    gidx = 1;
                    if ~isempty(inputs.group)
                        gidx = inputs.group.Value;
                    end
                    map = [inputs.erp1.Value inputs.erp2.Value inputs.erp3.Value inputs.erp4.Value];
                    T = psyrat_dod_varcompt('psyrat_data', pd, 'groupidx', gidx, 'map', map, ...
                        'savefile', 0, 'gui', 0);
                otherwise
                    error('G110Harness:kind', 'No callback table for kind "%s".', spec.kind);
            end
        end

        function [n, firstText] = shownRows(fig, label)
            %Rows of the table a viewer shows on screen when LABEL is pressed
            %(the same table its Save button writes), and the first text cell
            %of its first row. The window is closed.
            before = findall(groot, 'Type', 'figure');
            G110WorkbookHarness.pressButton(fig, label);
            after = findall(groot, 'Type', 'figure');
            opened = after(~ismember(after, before));
            if numel(opened) ~= 1
                delete(opened);
                error('G110Harness:show', '"%s" opened %d windows; expected 1.', label, numel(opened));
            end
            t = findall(opened, 'Type', 'uitable');
            if numel(t) ~= 1
                delete(opened);
                error('G110Harness:show', '"%s" shows %d tables; expected 1.', label, numel(t));
            end
            data = get(t, 'Data');
            n = size(data, 1);
            firstText = '';
            if n > 0 && iscell(data)
                firstText = G110WorkbookHarness.firstTextOf(data(1, :));
            end
            delete(opened);
        end

        function needle = bodyNeedle(spec, fig)
            %The text of one required BODY cell of a route's export, taken
            %from the callback's own payload or the fixture, never from the
            %exporter: the first text value in the first row of the
            %callback's table ('table'), of the DoD table the screen's
            %selections build ('dodexport'), or of the table the viewer shows
            %on screen ('viewer'); for Save IDs, the first ID the fixture puts
            %in the first ID row. '' when there is none.
            switch spec.kind
                case {'table', 'dodexport'}
                    T = G110WorkbookHarness.expectedTable(spec, fig);
                    row = cell(1, width(T));
                    for c = 1:width(T)
                        row{c} = T{1, c};
                    end
                    needle = G110WorkbookHarness.firstTextOf(row);
                case 'ids'
                    cb = get(G110WorkbookHarness.findButton(fig, spec.save), 'Callback');
                    blk = G110WorkbookHarness.expectedIdBlock(cb{2});
                    needle = G110WorkbookHarness.firstTextOf(blk(1, :));
                case 'viewer'
                    [~, needle] = G110WorkbookHarness.shownRows(fig, spec.show);
                otherwise
                    needle = '';
            end
        end

        function r = columnNameRowOf(file)
            %Row of a workbook's column-name row (see columnNameRow).
            r = G110WorkbookHarness.columnNameRow(G110WorkbookHarness.content(file));
        end

        function damageTarget(dest, how)
            %Called by the movefile/copyfile seam in the damage modes, before
            %it reports the failure: a replacement that changed ('damage') or
            %removed ('remove') the chosen file part-way and then failed; and,
            %in the wrongButReportsSuccess mode, after the forwarded move,
            %which the seam then reports as a success ('damage'). Uses only
            %fopen/fwrite/delete, which no seam shadows.
            dest = G110WorkbookHarness.textArg(dest);
            if strcmp(how, 'remove')
                if isfile(dest)
                    delete(dest);
                end
            else
                fid = fopen(dest, 'w');
                if fid >= 0
                    fwrite(fid, uint8('G110 damaged by a failed replacement'), 'uint8');
                    fclose(fid);
                end
            end
        end

        function t = dialogText(figs)
            %Every text string shown in the windows FIGS (errordlg keeps its
            %message in an axes text object), joined by line breaks.
            parts = {};
            for k = 1:numel(figs)
                if ~isvalid(figs(k))
                    continue;
                end
                groups = {findall(figs(k), 'Type', 'text'), findall(figs(k), 'Style', 'text')};
                for g = 1:numel(groups)
                    objs = groups{g};
                    for m = 1:numel(objs)
                        parts{end+1} = G110WorkbookHarness.joinText(get(objs(m), 'String')); %#ok<AGROW>
                    end
                end
            end
            t = strjoin(parts, newline);
        end

        function mark = dialogMark()
            %What is open before a call: every figure, and every message of
            %the msgbox family (errordlg, warndlg, msgbox keep the message in
            %a text object tagged 'MessageBox'). Pass it to dialogsSince.
            mark = struct();
            mark.figs = findall(groot, 'Type', 'figure');
            mark.msgs = findall(groot, 'Type', 'text', 'Tag', 'MessageBox');
        end

        function [t, names] = dialogsSince(mark, closeOwned)
            %The messages shown since MARK, and the names of their windows.
            %A message counts when its text object did not exist at MARK.
            %That covers a window msgbox REUSED: a modal msgbox with the same
            %title takes over the open window and deletes and recreates its
            %contents (msgbox.m), so its message object is new and current,
            %while a window merely left open from before is not read again.
            %Windows opened since MARK are named, and read in full when they
            %hold no msgbox message. With CLOSEOWNED true the windows opened
            %since MARK are closed after they are read; windows that were
            %already open are left as they are.
            if nargin < 2
                closeOwned = false;
            end
            parts = {};
            names = {};
            msgs = findall(groot, 'Type', 'text', 'Tag', 'MessageBox');
            fresh = msgs(~ismember(msgs, mark.msgs));
            for k = 1:numel(fresh)
                parts{end+1} = G110WorkbookHarness.joinText(get(fresh(k), 'String')); %#ok<AGROW>
                owner = ancestor(fresh(k), 'figure');
                if ~isempty(owner)
                    names{end+1} = char(get(owner, 'Name')); %#ok<AGROW>
                end
            end
            figs = findall(groot, 'Type', 'figure');
            opened = figs(~ismember(figs, mark.figs));
            for k = 1:numel(opened)
                names{end+1} = char(get(opened(k), 'Name')); %#ok<AGROW>
                if isempty(findall(opened(k), 'Type', 'text', 'Tag', 'MessageBox'))
                    seen = G110WorkbookHarness.dialogText(opened(k));
                    if ~isempty(seen)
                        parts{end+1} = seen; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(names)
                names = unique(names, 'stable');
            end
            t = strjoin(parts, newline);
            if closeOwned && ~isempty(opened)
                delete(opened(isvalid(opened)));
            end
        end

        function t = xmlText(file)
            %Every XML part of a workbook, joined (see staleReport).
            xml = G110WorkbookHarness.workbookXml(file);
            t = xml.all;
        end

        function corruptWritten(fn, args)
            %Called by the writecell/writetable seam right after a forwarded
            %write when a corruption mode is armed. Finds the needle in THIS
            %call's payload (first match, row by row), turns its position into
            %a sheet address from the call's own 'Range' origin (A1 when none)
            %plus, for writetable, the variable-name row, and overwrites that
            %one cell of the file just written with the marker. Nothing is
            %raised and the workbook stays valid. The seam folder is off the
            %path while this runs, so the writecell below is the real one.
            %CSV and 'append' writes are left alone.
            s = getappdata(0, G110WorkbookHarness.SeamKey);
            if isempty(s) || ~s.armed || numel(args) < 2 || ~any(strcmp(s.writerMode, ...
                    {'corruptHeader', 'corruptBody', 'corruptNumeric'}))
                return;
            end
            numeric = strcmp(s.writerMode, 'corruptNumeric');
            dest = G110WorkbookHarness.textArg(args{2});
            if ~endsWith(lower(dest), '.xlsx') || ~isfile(dest)
                return;
            end
            origin = 'A1';
            names = strcmp(fn, 'writetable');   % writetable writes variable names first by default
            opts = args(3:end);
            for k = 1:2:numel(opts) - 1
                switch lower(G110WorkbookHarness.textArg(opts{k}))
                    case 'range'
                        origin = G110WorkbookHarness.textArg(opts{k + 1});
                    case 'writevariablenames'
                        names = logical(opts{k + 1});
                    case 'writemode'
                        if strcmpi(G110WorkbookHarness.textArg(opts{k + 1}), 'append')
                            return;
                        end
                end
            end
            if numeric
                [i, j, x] = G110WorkbookHarness.findNumberInPayload(args{1});
                value = x + eps(x);   % the next larger double: equal under any tolerance
                what = sprintf('%.17g -> %.17g', x, value);
            else
                [i, j] = G110WorkbookHarness.findInPayload(args{1}, s.corruptNeedle);
                value = s.corruptMarker;
                what = s.corruptNeedle;
            end
            if isempty(i)
                return;
            end
            [r0, c0] = G110WorkbookHarness.cellRef(origin);
            address = sprintf('%s%d', G110WorkbookHarness.colName(c0 + j - 1), ...
                r0 + double(names) + i - 1);
            writecell({value}, dest, 'Range', address);
            s.corruptApplied = s.corruptApplied + 1;
            s.corruptLog(end+1, :) = {fn, dest, address, what};
            setappdata(0, G110WorkbookHarness.SeamKey, s);
        end

        function args = shiftedArgs(fn, args)
            %Called by the writecell/writetable seam in the 'shift' mode,
            %before it forwards the call: the call's 'Range' origin (A1 when
            %it gives none) moves down Shift(1) rows and right Shift(2)
            %columns, so every write of one save lands the same distance from
            %where it was aimed. The move is logged; nothing else in the call
            %changes. A range such as 'A3:D9' becomes its moved first cell,
            %which is how writecell and writetable place a payload anyway.
            s = getappdata(0, G110WorkbookHarness.SeamKey);
            dest = '';
            if numel(args) >= 2
                dest = G110WorkbookHarness.textArg(args{2});
            end
            origin = 'A1';
            at = 0;   % position of the 'Range' value in ARGS; 0 when absent
            for k = 3:2:numel(args) - 1
                if strcmpi(G110WorkbookHarness.textArg(args{k}), 'range')
                    origin = G110WorkbookHarness.textArg(args{k + 1});
                    at = k + 1;
                end
            end
            [r0, c0] = G110WorkbookHarness.cellRef(origin);
            moved = sprintf('%s%d', G110WorkbookHarness.colName(c0 + s.shift(2)), r0 + s.shift(1));
            if at > 0
                args{at} = moved;
            else
                args = [args, {'Range', moved}];
            end
            s.shiftApplied = s.shiftApplied + 1;
            s.shiftLog(end+1, :) = {fn, dest, origin, moved};
            setappdata(0, G110WorkbookHarness.SeamKey, s);
        end

        function extraWritten(fn, args)
            %Called by the writecell/writetable seam in the 'extra' mode right
            %after a forwarded write into env.work. It adds to the file just
            %written what the Extra option names: 'sheet', a sheet G110_Extra
            %holding one cell (the first sheet is untouched); a cell
            %reference, G110_EXTRA in that cell of the first sheet; 'none',
            %nothing (the control, whose call is still counted). The payload
            %stays where it was aimed and the workbook stays valid. The seam
            %folder is off the path while this runs, so the writecell below is
            %the real one.
            s = getappdata(0, G110WorkbookHarness.SeamKey);
            if isempty(s) || ~s.armed || ~strcmp(s.writerMode, 'extra') || numel(args) < 2
                return;
            end
            dest = G110WorkbookHarness.textArg(args{2});
            if ~endsWith(lower(dest), '.xlsx') || ~isfile(dest)
                return;
            end
            s.extraCalls = s.extraCalls + 1;
            what = '';
            if strcmpi(s.extra, 'sheet')
                state = warning('off', 'MATLAB:xlswrite:AddSheet');
                restore = onCleanup(@() warning(state)); %#ok<NASGU>
                writecell({'G110_EXTRA'}, dest, 'Sheet', 'G110_Extra');
                s.extraSheets = numel(sheetnames(dest));
                what = sprintf('sheet G110_Extra (%d sheets now)', s.extraSheets);
            elseif ~strcmpi(s.extra, 'none')
                writecell({'G110_EXTRA'}, dest, 'Range', s.extra);
                what = sprintf('G110_EXTRA in %s', s.extra);
            end
            if ~isempty(what)
                s.extraApplied = s.extraApplied + 1;
                s.extraLog(end+1, :) = {fn, dest, what};
            end
            setappdata(0, G110WorkbookHarness.SeamKey, s);
        end

        function blk = expectedIdBlock(pd)
            %The ID block a Save IDs export must hold, derived from the
            %fixture's groups: Good IDs and Bad IDs each listed top-down in
            %group order, the group name beside each ID ('none' for the
            %unnamed single group), blank where one list is shorter.
            good = cell(0, 1); goodGroup = cell(0, 1);
            bad = cell(0, 1); badGroup = cell(0, 1);
            for j = 1:numel(pd.relsummary.group)
                g = pd.relsummary.group(j);
                name = char(g.name);
                if isempty(name)
                    name = 'none';
                end
                gi = cellstr(string(g.goodids(:)));
                bi = cellstr(string(g.badids(:)));
                good = [good; gi(:)]; %#ok<AGROW>
                goodGroup = [goodGroup; repmat({name}, numel(gi), 1)]; %#ok<AGROW>
                bad = [bad; bi(:)]; %#ok<AGROW>
                badGroup = [badGroup; repmat({name}, numel(bi), 1)]; %#ok<AGROW>
            end
            n = max(numel(good), numel(bad));
            blk = repmat({''}, n, 4);
            blk(1:numel(good), 1) = good;
            blk(1:numel(bad), 2) = bad;
            blk(1:numel(good), 3) = goodGroup;
            blk(1:numel(bad), 4) = badGroup;
        end

        function [blk, hdr] = idBlock(file)
            %The rows below a Save IDs export's "Good IDs" header row.
            C = G110WorkbookHarness.content(file);
            r = find(strcmp(C(:, 1), 'Good IDs'), 1);
            if isempty(r)
                error('G110Harness:ids', 'No "Good IDs" header row in %s.', file);
            end
            hdr = C(r, :);
            blk = C(r + 1:end, :);
        end

        % ---------------------------------------------------------------
        % Fixtures (the shared factory plus the fields the export headers read)
        % ---------------------------------------------------------------

        function pd = withHeader(pd)
            %The provenance fields the Save headers read, which the
            %summary-only factory data lacks (as TestExportCallbacks adds).
            pd.ver = 'test';
            pd.rel.filename = 'workbook_replacement_fixture.psyrat';
            pd.rel.nchains = 2;
            pd.rel.niter = 600;
            pd.rel.seed = 12345;
            pd.rel.engine = 'cmdstan';
            pd.rel.priors = psyrat_default_priors;
            pd.rel.out.conv.converged = 1;
            pd.rel.out.conv.ndivergent = 0;
        end

        function pd = icData()
            pd = G110WorkbookHarness.withHeader(PsyRATTestDataFactory.makeICSummaryDataWithTrial());
        end

        function pd = trtData()
            pd = G110WorkbookHarness.withHeader(PsyRATTestDataFactory.makeTRTSummaryData());
        end

        function pd = dodData()
            pd = G110WorkbookHarness.withHeader(PsyRATTestDataFactory.makeDoDRelData());
        end

        function map = dodMap(rel)
            %Event indices of the DoD fixture's stored ERP_1..ERP_4 mapping.
            events = cellstr(string(rel.events(:)));
            names = cellstr(string(rel.dod_map(:)));
            [~, map] = ismember(names, events);
            map = map(:)';
        end

        function prefs = criterionPrefs()
            %Criterion table only (no plot), one-facet, cutoff 0.5.
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 0;
            prefs.view.tablecriterion = 1;
            prefs.view.criterioncutoff = 0.5;
            prefs.view.reltype = 1;
            prefs.view.gcoeff = 2;
            prefs.view.depcentmeas = 1;
            prefs.view.relcentmeas = 1;
        end

        function pd = idFixture(spec, shape)
            %Two groups on the route's own fixture. 'long': 10 Good IDs and 4
            %Bad IDs, with s01 and s02 in both groups (the group column is
            %what tells them apart). 'zeroGood': no Good IDs in either group
            %and 4 Bad IDs, so the export is shorter in every column.
            pd = spec.data();
            pd.rel.groups = {'GroupA'; 'GroupB'};
            g = pd.relsummary.group(1);
            switch shape
                case 'long'
                    goodA = {'s01'; 's02'; 's03'; 's04'; 's05'}; badA = {'s06'; 's07'};
                    goodB = {'s01'; 's02'; 's08'; 's09'; 's10'}; badB = {'s11'; 's12'};
                case 'zeroGood'
                    goodA = cell(0, 1); badA = {'s01'; 's02'; 's03'};
                    goodB = cell(0, 1); badB = {'s01'};
                otherwise
                    error('G110Harness:ids', 'Unknown ID fixture shape "%s".', shape);
            end
            A = g; A.name = 'GroupA'; A.goodids = goodA; A.badids = badA; A.event.goodn = numel(goodA);
            B = g; B.name = 'GroupB'; B.goodids = goodB; B.badids = badB; B.event.goodn = numel(goodB);
            pd.relsummary.group = [A B];
        end

        % ---------------------------------------------------------------
        % Small utilities (public so the generated seam gate can use them)
        % ---------------------------------------------------------------

        function t = textArg(x)
            %A path argument as char ('' when it is not text).
            t = '';
            if ischar(x)
                t = x;
            elseif isstring(x) && isscalar(x)
                t = char(x);
            end
        end

        function n = nameOf(p)
            %File name with extension of a path ('' when P is not text).
            [~, base, ext] = fileparts(G110WorkbookHarness.textArg(p));
            n = [base ext];
        end

        function tf = sameFolder(a, b)
            %Folder equality after dropping trailing separators and the
            %macOS /private prefix of the temporary folders.
            tf = strcmp(G110WorkbookHarness.normFolder(a), G110WorkbookHarness.normFolder(b));
        end

        function b = fileBytes(file)
            fid = fopen(file, 'r');
            if fid < 0
                error('G110Harness:read', 'Cannot read %s.', file);
            end
            closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
            b = fread(fid, inf, '*uint8');
        end

        function t = maskTime(t)
            %datestr(clock) text, as the export headers write it
            t = regexprep(t, '\d{1,2}-[A-Za-z]{3}-\d{4} \d{2}:\d{2}:\d{2}', '<TIME>');
        end

        function writeLines(file, lines)
            fid = fopen(file, 'w');
            if fid < 0
                error('G110Harness:write', 'Cannot write %s.', file);
            end
            closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fprintf(fid, '%s\n', lines{:});
        end
    end

    methods (Static, Access = private)

        function fig = openNamed(builder, name)
            delete(findall(groot, 'Type', 'figure', 'Name', name));
            builder();
            figs = findall(groot, 'Type', 'figure', 'Name', name);
            if isempty(figs)
                error('G110Harness:open', 'The "%s" window did not open.', name);
            end
            fig = figs(1);
        end

        function fig = openChild(parent, label)
            %Press LABEL on PARENT and return the one output window it opens.
            before = findall(groot, 'Type', 'figure');
            G110WorkbookHarness.pressButton(parent, label);
            after = findall(groot, 'Type', 'figure');
            opened = after(~ismember(after, before));
            keep = false(size(opened));
            for k = 1:numel(opened)
                keep(k) = strcmp(get(opened(k), 'Tag'), 'psyrat_output');
            end
            opened = opened(keep);
            if numel(opened) ~= 1
                error('G110Harness:open', '"%s" opened %d output windows; expected 1.', ...
                    label, numel(opened));
            end
            fig = opened(1);
        end

        function C = canon(raw)
            C = cell(size(raw));
            for k = 1:numel(raw)
                v = raw{k};
                if isa(v, 'missing') || (isnumeric(v) && isempty(v))
                    C{k} = '';
                elseif ischar(v) || isstring(v)
                    C{k} = G110WorkbookHarness.maskTime(char(v));
                elseif isa(v, 'datetime')
                    C{k} = '<TIME>';
                elseif islogical(v)
                    C{k} = mat2str(v);
                elseif isnumeric(v)
                    C{k} = sprintf('%.15g', v);
                else
                    C{k} = ['<' class(v) '>'];
                end
            end
        end

        function c = column1(T)
            %The first column of a table as a column cell array.
            v = T{:, 1};
            if iscell(v)
                c = v;
            elseif isstring(v) || iscategorical(v)
                c = cellstr(string(v));
            elseif ischar(v)
                c = cellstr(v);
            else
                c = num2cell(v);
            end
            c = c(:);
        end

        function r = columnNameRow(C)
            %The first row with two or more filled cells. Header lines are
            %written one per row in column A, so this is the table's
            %column-name row (or the ID export's "Good IDs" row).
            r = [];
            for k = 1:size(C, 1)
                if nnz(~cellfun(@(x) isempty(strtrim(x)), C(k, :))) >= 2
                    r = k;
                    return;
                end
            end
        end

        function v = namedArg(args, name)
            for k = 1:2:numel(args) - 1
                if ischar(args{k}) && strcmp(args{k}, name)
                    v = args{k + 1};
                    return;
                end
            end
            error('G110Harness:arg', 'The callback carries no "%s" argument.', name);
        end

        function t = joinText(s)
            %A String property (char, char matrix, cellstr or string) as one
            %char row, lines joined by line breaks.
            if ischar(s) && size(s, 1) > 1
                s = cellstr(s);
            end
            if iscell(s) || isstring(s)
                s = strjoin(cellstr(s(:))', newline);
            end
            t = char(s);
        end

        function t = payloadText(v)
            %One payload value as text ('' when it is not text): a 1x1 cell
            %is unwrapped; char, scalar string and scalar categorical count.
            if iscell(v) && isscalar(v)
                v = v{1};
            end
            t = '';
            if ischar(v)
                t = v;
            elseif (isstring(v) || iscategorical(v)) && isscalar(v)
                t = char(string(v));
            end
        end

        function t = firstTextOf(values)
            %The first non-blank text value in a cell row.
            t = '';
            for k = 1:numel(values)
                candidate = G110WorkbookHarness.payloadText(values{k});
                if ~isempty(strtrim(candidate))
                    t = candidate;
                    return;
                end
            end
        end

        function [i, j] = findInPayload(payload, needle)
            %Row and column of the first payload cell whose text is exactly
            %NEEDLE, scanning row by row; [] when absent. A table is scanned by
            %rows and variables (the census tables have one column per
            %variable).
            i = [];
            j = [];
            if iscell(payload)
                nr = size(payload, 1);
                nc = size(payload, 2);
                at = @(r, c) payload{r, c};
            elseif istable(payload)
                nr = height(payload);
                nc = width(payload);
                at = @(r, c) payload{r, c};
            else
                return;
            end
            for r = 1:nr
                for c = 1:nc
                    if strcmp(G110WorkbookHarness.payloadText(at(r, c)), needle)
                        i = r;
                        j = c;
                        return;
                    end
                end
            end
        end

        function [i, j, x] = findNumberInPayload(payload)
            %Row, column and value of the first finite real number (not a
            %logical) in a payload, scanning row by row; [] when there is none.
            i = [];
            j = [];
            x = [];
            if iscell(payload)
                nr = size(payload, 1);
                nc = size(payload, 2);
            elseif istable(payload)
                nr = height(payload);
                nc = width(payload);
            else
                return;
            end
            for r = 1:nr
                for c = 1:nc
                    v = payload{r, c};
                    if iscell(v) && isscalar(v)
                        v = v{1};
                    end
                    if isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v)
                        i = r;
                        j = c;
                        x = double(v);
                        return;
                    end
                end
            end
        end

        function [row, col] = cellRef(ref)
            %Row and column of a cell reference such as 'A22' (the first
            %cell of a range such as 'A22:D40').
            tok = regexp(char(ref), '^\s*([A-Za-z]+)(\d+)', 'tokens', 'once');
            if isempty(tok)
                error('G110Harness:range', 'Cannot read the cell reference "%s".', char(ref));
            end
            col = 0;
            letters = upper(tok{1});
            for k = 1:numel(letters)
                col = col * 26 + (double(letters(k)) - double('A') + 1);
            end
            row = str2double(tok{2});
        end

        function s = colName(c)
            %Spreadsheet column letters of column number C (1 -> A, 27 -> AA).
            s = '';
            while c > 0
                m = mod(c - 1, 26);
                s = [char(double('A') + m) s]; %#ok<AGROW>
                c = floor((c - 1) / 26);
            end
        end

        function xml = workbookXml(file)
            %Unzip a copy of the workbook (byte copy, so no seam is involved)
            %and return every XML part joined, plus xl/workbook.xml.
            tmp = tempname;
            mkdir(tmp);
            cleaner = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
            zipfile = fullfile(tmp, 'book.zip');
            fid = fopen(zipfile, 'w');
            fwrite(fid, G110WorkbookHarness.fileBytes(file), 'uint8');
            fclose(fid);
            parts = unzip(zipfile, fullfile(tmp, 'parts'));
            xml = struct('all', '', 'workbook', '');
            chunks = {};
            for k = 1:numel(parts)
                part = strrep(parts{k}, '\', '/');
                if endsWith(part, '.xml') || endsWith(part, '.rels')
                    text = fileread(parts{k});
                    chunks{end+1} = text; %#ok<AGROW>
                    if endsWith(part, 'xl/workbook.xml')
                        xml.workbook = text;
                    end
                end
            end
            xml.all = strjoin(chunks, newline);
        end

        function p = normFolder(p)
            p = G110WorkbookHarness.textArg(p);
            while numel(p) > 1 && (p(end) == '/' || p(end) == '\')
                p(end) = [];
            end
            if startsWith(p, '/private/var/') || startsWith(p, '/private/tmp/')
                p = p(9:end);
            end
        end

        function dropFromPath(d)
            if ~isempty(d) && any(strcmp(strsplit(path, pathsep), d))
                rmpath(d);
            end
        end

        function removeFolder(d)
            if ~isempty(d) && isfolder(d)
                try
                    rmdir(d, 's');
                catch ME
                    warning('G110Harness:cleanup', 'Could not remove %s: %s', d, ME.message);
                end
            end
        end

        function writeDialogDoubles(d)
            G110WorkbookHarness.writeLines(fullfile(d, 'uiputfile.m'), { ...
                'function [filename, pathname, filterindex] = uiputfile(varargin)'
                '%TEST DOUBLE (G110 TestWorkbookReplacement): answers from a global; 0 is Cancel.'
                'global G110_SAVE_ANSWER G110_SAVE_CALLS'
                'if isempty(G110_SAVE_CALLS), G110_SAVE_CALLS = 0; end'
                'G110_SAVE_CALLS = G110_SAVE_CALLS + 1;'
                'if isempty(G110_SAVE_ANSWER) || isequal(G110_SAVE_ANSWER, 0)'
                '    filename = 0; pathname = 0; filterindex = 0;'
                '    return;'
                'end'
                '[p, n, e] = fileparts(char(G110_SAVE_ANSWER));'
                'filename = [n e];'
                'pathname = [p filesep];'
                'filterindex = 1;'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'questdlg.m'), { ...
                'function answer = questdlg(varargin)'
                '%TEST DOUBLE (G110): never blocks; records the call and answers as a closed dialog.'
                'global G110_QUESTDLG_CALLS'
                'if isempty(G110_QUESTDLG_CALLS), G110_QUESTDLG_CALLS = {}; end'
                'G110_QUESTDLG_CALLS{end+1} = varargin;'
                'answer = '''';'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'uiwait.m'), { ...
                'function uiwait(varargin)'
                '%TEST DOUBLE (G110): never blocks, so an unexpected modal cannot hang a headless run.'
                'global G110_UIWAIT_CALLS'
                'if isempty(G110_UIWAIT_CALLS), G110_UIWAIT_CALLS = 0; end'
                'G110_UIWAIT_CALLS = G110_UIWAIT_CALLS + 1;'
                'end'});
        end

        function lines = writerSeam(fn)
            %Source of the writecell/writetable double named FN.
            lines = { ...
                sprintf('function %s(varargin)', fn)
                sprintf('%%TEST SEAM (G110 TestWorkbookReplacement): logs a %s call made by', fn)
                '%PsyRAT code, fails it on request, and otherwise forwards to the real function.'
                sprintf('action = g110_seam_gate(''writer'', ''%s'', varargin);', fn)
                'if strcmp(action, ''before'')'
                sprintf('    error(''G110Test:injectedWriterFailure'', ''Injected %s failure before writing (G110 test seam).'');', fn)
                'end'
                'if strcmp(action, ''shift'')'
                sprintf('    varargin = G110WorkbookHarness.shiftedArgs(''%s'', varargin);', fn)
                'end'
                'd = fileparts(mfilename(''fullpath''));'
                'rmpath(d);'
                'restore = onCleanup(@() addpath(d)); %#ok<NASGU>'
                sprintf('%s(varargin{:});', fn)
                'if strcmp(action, ''after'')'
                sprintf('    error(''G110Test:injectedWriterFailure'', ''Injected %s failure after writing (G110 test seam).'');', fn)
                'end'
                'if strcmp(action, ''corrupt'')'
                sprintf('    G110WorkbookHarness.corruptWritten(''%s'', varargin);', fn)
                'end'
                'if strcmp(action, ''extra'')'
                sprintf('    G110WorkbookHarness.extraWritten(''%s'', varargin);', fn)
                'end'
                'end'};
        end

        function lines = moverSeam(fn)
            %Source of the movefile/copyfile double named FN.
            lines = { ...
                sprintf('function varargout = %s(varargin)', fn)
                sprintf('%%TEST SEAM (G110 TestWorkbookReplacement): logs a %s call made by PsyRAT', fn)
                '%code, fails the first one onto the watched target on request (the way the real'
                '%function reports a failure), and otherwise forwards to the real function. In the'
                '%wrongButReportsSuccess mode the first one is forwarded, the target is then'
                '%overwritten with other bytes, and the forwarded success is returned.'
                sprintf('action = g110_seam_gate(''replace'', ''%s'', varargin);', fn)
                'if any(strcmp(action, {''damage'', ''remove''}))'
                '    G110WorkbookHarness.damageTarget(varargin{2}, action);'
                '    action = ''fail'';'
                'end'
                'if strcmp(action, ''fail'')'
                '    msg = ''Injected final-replacement failure (G110 test seam).'';'
                '    if nargout == 0'
                '        error(''G110Test:injectedReplaceFailure'', ''%s'', msg);'
                '    end'
                '    out = {false, msg, ''G110Test:injectedReplaceFailure''};'
                '    varargout = out(1:nargout);'
                '    return;'
                'end'
                'd = fileparts(mfilename(''fullpath''));'
                'rmpath(d);'
                'restore = onCleanup(@() addpath(d)); %#ok<NASGU>'
                'if nargout == 0'
                sprintf('    %s(varargin{:});', fn)
                'else'
                sprintf('    [varargout{1:nargout}] = %s(varargin{:});', fn)
                'end'
                'if strcmp(action, ''wrongSuccess'')'
                '    G110WorkbookHarness.damageTarget(varargin{2}, ''damage'');'
                'end'
                'end'};
        end

        function writeSeams(d)
            %The four I/O doubles and their shared gate. Each double asks the
            %gate what to do, then either fails or forwards to the real
            %function with its own folder off the path for that call.
            G110WorkbookHarness.writeLines(fullfile(d, 'writecell.m'), G110WorkbookHarness.writerSeam('writecell'));
            G110WorkbookHarness.writeLines(fullfile(d, 'writetable.m'), G110WorkbookHarness.writerSeam('writetable'));
            G110WorkbookHarness.writeLines(fullfile(d, 'movefile.m'), G110WorkbookHarness.moverSeam('movefile'));
            G110WorkbookHarness.writeLines(fullfile(d, 'copyfile.m'), G110WorkbookHarness.moverSeam('copyfile'));
            G110WorkbookHarness.writeLines(fullfile(d, 'g110_seam_gate.m'), { ...
                'function action = g110_seam_gate(kind, fn, args)'
                '%TEST SEAM (G110 TestWorkbookReplacement): the shared decision for the writecell,'
                '%writetable, movefile and copyfile doubles. Calls made from inside MATLAB''s own'
                '%toolboxes are forwarded uncounted; every other call is logged with its caller.'
                'action = ''forward'';'
                'if ~isappdata(0, ''G110_io_seam'')'
                '    return;'
                'end'
                's = getappdata(0, ''G110_io_seam'');'
                'if ~s.armed'
                '    return;'
                'end'
                'st = dbstack(''-completenames'');'
                '%st(1) is this gate, st(2) the double, st(3) the code that called the double'
                'if numel(st) < 3 || startsWith(st(3).file, matlabroot)'
                '    return;'
                'end'
                '[~, cn, ce] = fileparts(st(3).file);'
                'caller = sprintf(''%s:%d (%s)'', [cn ce], st(3).line, st(3).name);'
                'dest = '''';'
                'if numel(args) >= 2'
                '    dest = G110WorkbookHarness.textArg(args{2});'
                'end'
                '[folder, name, ext] = fileparts(dest);'
                'inwork = G110WorkbookHarness.sameFolder(folder, s.workDir);'
                'if strcmp(kind, ''writer'')'
                '    s.writerLog(end+1, :) = {fn, dest, caller};'
                '    if inwork'
                '        s.writerCalls = s.writerCalls + 1;'
                '        if ~s.writerInjected && s.writerCalls == s.writerFailAt && any(strcmp(s.writerMode, {''before'', ''after''}))'
                '            action = s.writerMode;'
                '            s.writerInjected = true;'
                '        elseif any(strcmp(s.writerMode, {''corruptHeader'', ''corruptBody'', ''corruptNumeric''}))'
                '            action = ''corrupt'';'
                '        elseif strcmp(s.writerMode, ''shift'')'
                '            action = ''shift'';'
                '        elseif strcmp(s.writerMode, ''extra'')'
                '            action = ''extra'';'
                '        end'
                '    end'
                'else'
                '    source = '''';'
                '    if ~isempty(args)'
                '        source = G110WorkbookHarness.textArg(args{1});'
                '    end'
                '    s.moveLog(end+1, :) = {fn, source, dest, caller};'
                '    if inwork && strcmp([name ext], s.targetName)'
                '        s.replaceCalls = s.replaceCalls + 1;'
                '        if strcmp(s.replaceMode, ''failFirst'') && ~s.replaceInjected'
                '            action = ''fail'';'
                '            s.replaceInjected = true;'
                '        elseif any(strcmp(s.replaceMode, {''damageFirst'', ''removeFirst'', ''damageAll''}))'
                '            if ~s.replaceInjected'
                '                action = ''damage'';'
                '                if strcmp(s.replaceMode, ''removeFirst'')'
                '                    action = ''remove'';'
                '                end'
                '                s.replaceInjected = true;'
                '            elseif strcmp(s.replaceMode, ''damageAll'')'
                '                action = ''fail'';'
                '            end'
                '        elseif strcmp(s.replaceMode, ''wrongButReportsSuccess'') && ~s.replaceInjected'
                '            action = ''wrongSuccess'';'
                '            s.replaceInjected = true;'
                '        end'
                '    end'
                'end'
                'setappdata(0, ''G110_io_seam'', s);'
                'end'});
        end

        function writeFaults(d)
            %The fread/ferror/fwrite/delete doubles and their gate. They are
            %built-ins, so each forwards with builtin() (no path changes per
            %call); only calls whose immediate caller is the aimed file are
            %ever changed.
            G110WorkbookHarness.writeLines(fullfile(d, 'fread.m'), { ...
                'function varargout = fread(varargin)'
                '%TEST SEAM (G110 TestWorkbookReplacementHelper): the built-in fread, except that'
                '%the first read of the aimed file by the aimed caller returns half its bytes.'
                'n = max(nargout, 1);'
                '[varargout{1:n}] = builtin(''fread'', varargin{:});'
                'if g110_fault_gate(''read'', varargin)'
                '    keep = floor(numel(varargout{1}) / 2);'
                '    varargout{1} = varargout{1}(1:keep);'
                '    if n >= 2'
                '        varargout{2} = keep;'
                '    end'
                'end'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'ferror.m'), { ...
                'function varargout = ferror(varargin)'
                '%TEST SEAM (G110): reports an error for the file id whose read was cut;'
                '%otherwise the built-in ferror.'
                'msg = g110_fault_gate(''ferror'', varargin);'
                'if ~isempty(msg)'
                '    out = {msg, -1};'
                '    varargout = out(1:max(nargout, 1));'
                '    return;'
                'end'
                '[varargout{1:max(nargout, 1)}] = builtin(''ferror'', varargin{:});'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'fwrite.m'), { ...
                'function varargout = fwrite(varargin)'
                '%TEST SEAM (G110): the built-in fwrite, except that the first write by the aimed'
                '%caller into an aimed file writes half the data and then raises.'
                'if g110_fault_gate(''write'', varargin)'
                '    data = varargin{2};'
                '    builtin(''fwrite'', varargin{1}, data(1:floor(numel(data) / 2)), varargin{3:end});'
                '    error(''G110Test:injectedWriteFailure'', ''Injected write failure after a partial write (G110 test seam).'');'
                'end'
                '[varargout{1:max(nargout, 1)}] = builtin(''fwrite'', varargin{:});'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'delete.m'), { ...
                'function delete(varargin)'
                '%TEST SEAM (G110): the built-in delete, except that deleting an aimed file by the'
                '%aimed caller is refused the way MATLAB refuses it: a warning, and the file stays.'
                'if g110_fault_gate(''delete'', varargin)'
                '    warning(''MATLAB:DELETE:Permission'', ''File not found or permission denied'');'
                '    return;'
                'end'
                'builtin(''delete'', varargin{:});'
                'end'});
            G110WorkbookHarness.writeLines(fullfile(d, 'g110_fault_gate.m'), { ...
                'function out = g110_fault_gate(kind, args)'
                '%TEST SEAM (G110): the decision for the fread, ferror, fwrite and delete doubles.'
                '%Only calls made directly from the aimed caller file, on files in the watched'
                '%folder, are considered; everything else is forwarded unchanged.'
                'out = false;'
                'if strcmp(kind, ''ferror'')'
                '    out = '''';'
                'end'
                'if ~isappdata(0, ''G110_fault_seam'') || isempty(args)'
                '    return;'
                'end'
                's = getappdata(0, ''G110_fault_seam'');'
                'if ~s.armed'
                '    return;'
                'end'
                'st = dbstack(''-completenames'');'
                '%st(1) is this gate, st(2) the double, st(3) the code that called the double'
                'if numel(st) < 3'
                '    return;'
                'end'
                '[~, cn, ce] = fileparts(st(3).file);'
                'if ~strcmpi([cn ce], s.caller)'
                '    return;'
                'end'
                'if strcmp(kind, ''delete'')'
                '    p = G110WorkbookHarness.textArg(args{1});'
                '    [folder, n, e] = fileparts(p);'
                '    if ~isempty(p) && ~isempty(s.refuseDelete) && G110WorkbookHarness.sameFolder(folder, s.workDir) ...'
                '            && any(startsWith([n e], s.refuseDelete))'
                '        out = true;'
                '        s.deleteRefused{end+1} = p;'
                '    end'
                'else'
                '    fid = args{1};'
                '    if ~isnumeric(fid) || ~isscalar(fid) || fid < 3'
                '        return;'
                '    end'
                '    p = fopen(fid);'
                '    [folder, n, e] = fileparts(p);'
                '    aimed = ~isempty(s.cutRead) && ~isempty(p) && strcmp([n e], G110WorkbookHarness.nameOf(s.cutRead)) ...'
                '        && G110WorkbookHarness.sameFolder(folder, fileparts(s.cutRead));'
                '    switch kind'
                '        case ''read'''
                '            if aimed && ~s.readCut'
                '                out = true;'
                '                s.readCut = true;'
                '                s.errorFids(end+1) = fid;'
                '            else'
                '                %a read that was not cut leaves no injected error on this'
                '                %id, as ferror reports the most recent operation; an id'
                '                %closed and reused for any file starts clean'
                '                s.errorFids(s.errorFids == fid) = [];'
                '            end'
                '        case ''ferror'''
                '            s.ferrorCalls = s.ferrorCalls + 1;'
                '            if aimed && any(s.errorFids == fid)'
                '                out = ''Injected read error (G110 test seam).'';'
                '            end'
                '        case ''write'''
                '            if ~s.writeFailed && ~isempty(s.failWrite) && G110WorkbookHarness.sameFolder(folder, s.workDir) ...'
                '                    && startsWith([n e], s.failWrite)'
                '                out = true;'
                '                s.writeFailed = true;'
                '            end'
                '    end'
                'end'
                'setappdata(0, ''G110_fault_seam'', s);'
                'end'});
        end

        function writeRecorders(d)
            %The errordlg recorder: records the call, then forwards to the real
            %errordlg with its own folder off the path for that call (errordlg
            %is a MATLAB file, not a built-in).
            G110WorkbookHarness.writeLines(fullfile(d, 'errordlg.m'), { ...
                'function varargout = errordlg(varargin)'
                '%TEST RECORDER (G110 TestWorkbookReplacementHelper): records every errordlg call'
                '%(caller and arguments), then shows the real dialog. Calls are counted, not'
                '%windows: a modal errordlg reuses an open dialog of the same title.'
                'st = dbstack(''-completenames'');'
                'caller = '''';'
                'if numel(st) >= 2'
                '    [~, cn, ce] = fileparts(st(2).file);'
                '    caller = sprintf(''%s:%d (%s)'', [cn ce], st(2).line, st(2).name);'
                'end'
                'calls = getappdata(0, ''G110_errordlg_calls'');'
                'if ~iscell(calls)'
                '    calls = {};'
                'end'
                'calls{end+1} = {caller, varargin};'
                'setappdata(0, ''G110_errordlg_calls'', calls);'
                'd = fileparts(mfilename(''fullpath''));'
                'rmpath(d);'
                'restore = onCleanup(@() addpath(d)); %#ok<NASGU>'
                'if nargout == 0'
                '    errordlg(varargin{:});'
                'else'
                '    [varargout{1:nargout}] = errordlg(varargin{:});'
                'end'
                'end'});
        end
    end
end
