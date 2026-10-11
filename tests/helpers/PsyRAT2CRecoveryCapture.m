classdef PsyRAT2CRecoveryCapture
    % Opt-in evidence capture for the bounded 2C live recovery cells.
    %
    % PURPOSE. The live recovery tests gate their inferential assertions on
    % convergence assumptions, so a non-converged or otherwise skipped run
    % leaves no evidence behind once the class teardown removes its scratch
    % CmdStan directory. This helper lets a central runner retain, per fit,
    % the REL struct, the planted truth, the simulated data table, and the raw
    % CmdStan output directory BEFORE any convergence assumption executes. Root
    % then adjudicates recovery, coverage and diagnostics from the retained
    % evidence rather than from a pass/fail line alone.
    %
    % ENVIRONMENT CONTRACT.
    %   PSYRAT_2C_RECOVERY_CAPTURE_DIR   absolute path of the capture root.
    %     Unset or empty  -> capture() is a complete no-op (returns
    %                        info.captured = false, touches nothing).
    %     Set             -> must be ABSOLUTE and must not overlap the fit's
    %                        cmdstanOutDir or its parent scratch root (the
    %                        integration cmdstan_artifacts directory) in either
    %                        direction. The root is created if absent; the
    %                        per-run destination <root>/<tag> MUST NOT already
    %                        exist. An existing destination is refused with an
    %                        error (never overwritten, never merged) so a rerun
    %                        cannot silently replace earlier evidence. Choose a
    %                        fresh root (or fresh tags) per execution.
    %   Nothing is ever deleted by this helper. Originals in the CmdStan output
    %   directory are left in place. When the capture env var is set, the
    %   integration classes' teardown SKIPS their scratch cleanup entirely so
    %   the raw scratch survives even a failed capture; root's external
    %   artifact guard archives/restores the owned scratch afterwards.
    %
    % SUCCESS CRITERIA. A capture is reported successful (info.captured = true,
    % meta.capture_status = 'ok') only when ALL of the following hold:
    %   - the source cmdstanOutDir exists;
    %   - at least one chain CSV (CmdStan convention output-*.csv) was copied;
    %   - every copied chain CSV carries a CmdStan header line beginning
    %     'lp__' followed by at least one non-comment sample row;
    %   - every copied file matches its source in size, and every copied
    %     chain CSV matches its source byte for byte.
    % Otherwise the REL MAT / data / truth / meta are still written (they must
    % survive), meta.capture_status = 'failed' with the reasons, and capture()
    % raises an error. A requested raw capture never succeeds on the REL MAT
    % alone. Numeric content of the CSVs is NOT validated (diagnostic and
    % partial outputs are preserved as they are).
    %
    % WHAT IS WRITTEN under <root>/<tag>/
    %   <tag>_rel.mat        rel, truth, dataTbl, meta   (MAT v7; falls back to
    %                        v7.3 only if v7 refuses the payload size)
    %   data.csv             the simulated long table (writetable)
    %   truth.json           planted truth (jsonencode; numeric arrays as lists)
    %   meta.json            tag, timestamps, MATLAB/toolbox versions, the
    %                        route fields of REL, caller-supplied meta, the
    %                        copy manifest with per-file source/copy sizes,
    %                        per-chain-CSV header/row evidence, and
    %                        capture_status ('ok' | 'failed') with reasons
    %   cmdstan_output/      filtered copy of the task-owned CmdStan output
    %                        directory: CSV / Stan / data / metadata / log
    %                        files are kept; executables and build
    %                        intermediates (.hpp, .o, .d, extensionless
    %                        binaries, ...) are skipped and listed in meta.json
    %   README.txt           a copy of this contract for the artifact reader
    %
    % API.
    %   info = PsyRAT2CRecoveryCapture.capture(tag, rel, truth, dataTbl, ...
    %              cmdstanOutDir, 'meta', metaStruct)
    %     tag           unique, filesystem-safe run label (char)
    %     rel           REL struct returned by psyrat_computevarcomp (pre- or
    %                   post-psyrat_checkconv; callers capture pre-assumption)
    %     truth         planted-truth struct from the test simulator
    %     dataTbl       the long table passed to psyrat_computevarcomp
    %     cmdstanOutDir the 'cmdstanoutdir' the fit was given (task-owned)
    %     'meta'        optional struct of caller metadata (settings, design)
    %   info.captured  logical; info.dest, info.matFile, info.copied,
    %   info.skipped, info.chainFiles when captured.
    %
    %   dirOrEmpty = PsyRAT2CRecoveryCapture.captureDir()
    %     Returns the configured root or '' when capture is disabled.
    %
    % FAILURE POLICY. When capture is enabled, any failure to write (refused
    % destination, unwritable root, save error) or any unmet success
    % criterion raises an error rather than returning quietly. A capture that
    % was requested and did not fully happen must not look like one that did.

    properties (Constant)
        EnvVar = 'PSYRAT_2C_RECOVERY_CAPTURE_DIR';
        % file extensions retained from the CmdStan output directory (lower case)
        KeepExtensions = {'.csv', '.stan', '.r', '.json', '.txt', '.log', ...
            '.out', '.err', '.dat', '.data', '.md'};
    end

    methods (Static)

        function d = captureDir()
            % Configured capture root, or '' when capture is disabled.
            d = getenv(PsyRAT2CRecoveryCapture.EnvVar);
            if isempty(d)
                d = '';
            else
                d = strtrim(d);
            end
        end

        function info = capture(tag, rel, truth, dataTbl, cmdstanOutDir, varargin)
            info = struct('captured', false, 'tag', tag, 'dest', '', ...
                'matFile', '', 'copied', {{}}, 'skipped', {{}}, 'chainFiles', {{}});

            root = PsyRAT2CRecoveryCapture.captureDir();
            if isempty(root)
                return;   % disabled: complete no-op
            end

            % ---- optional caller metadata ----
            callerMeta = struct();
            k = 1;
            while k <= numel(varargin)
                if ischar(varargin{k}) && strcmpi(varargin{k}, 'meta') && k < numel(varargin)
                    callerMeta = varargin{k+1};
                    k = k + 2;
                else
                    error('PsyRAT2CRecoveryCapture:badArgs', ...
                        'Unrecognized argument at position %d; only ''meta'', struct is accepted.', k);
                end
            end

            % ---- validate tag, root (absolute, non-overlapping) and destination ----
            if ~ischar(tag) || isempty(tag) || ~isempty(regexp(tag, '[^A-Za-z0-9_\-]', 'once'))
                error('PsyRAT2CRecoveryCapture:badTag', ...
                    'tag must be a nonempty char of [A-Za-z0-9_-]; got "%s".', char(string(tag)));
            end
            if ~ischar(cmdstanOutDir) || isempty(cmdstanOutDir)
                error('PsyRAT2CRecoveryCapture:badSource', ...
                    'cmdstanOutDir must be a nonempty char path.');
            end
            if isempty(regexp(root, '^(/|[A-Za-z]:[\\/]|\\\\)', 'once'))
                error('PsyRAT2CRecoveryCapture:rootNotAbsolute', ...
                    '%s must be an absolute path; got "%s".', PsyRAT2CRecoveryCapture.EnvVar, root);
            end
            % the capture root must not overlap the fit's output directory or its
            % parent scratch root (tests/integration/cmdstan_artifacts) in either
            % direction, so the copy can never nest inside or swallow its source
            scratchRoot = fileparts(PsyRAT2CRecoveryCapture.normalizePath(cmdstanOutDir));
            for guard = {PsyRAT2CRecoveryCapture.normalizePath(cmdstanOutDir), scratchRoot}
                if PsyRAT2CRecoveryCapture.pathsOverlap(PsyRAT2CRecoveryCapture.normalizePath(root), guard{1})
                    error('PsyRAT2CRecoveryCapture:rootOverlapsScratch', ...
                        ['Refusing to capture: %s "%s" overlaps the CmdStan scratch "%s" ', ...
                        '(same path, or one contains the other).'], ...
                        PsyRAT2CRecoveryCapture.EnvVar, root, guard{1});
                end
            end
            if ~isfolder(root)
                [okRoot, msgRoot] = mkdir(root);
                if ~okRoot
                    error('PsyRAT2CRecoveryCapture:rootUnwritable', ...
                        'Cannot create capture root "%s": %s', root, msgRoot);
                end
            end
            dest = fullfile(root, tag);
            if exist(dest, 'dir') || exist(dest, 'file')
                error('PsyRAT2CRecoveryCapture:destinationExists', ...
                    ['Refusing to capture: destination "%s" already exists. This helper ', ...
                    'never overwrites or merges evidence; use a fresh capture root or tag.'], dest);
            end
            [okDest, msgDest] = mkdir(dest);
            if ~okDest
                error('PsyRAT2CRecoveryCapture:destUnwritable', ...
                    'Cannot create capture destination "%s": %s', dest, msgDest);
            end

            % ---- metadata ----
            meta = struct();
            meta.tag = tag;
            meta.captured_local = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
            meta.captured_utc   = char(datetime('now', 'TimeZone', 'UTC', ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
            meta.matlab_version = version();
            meta.computer = computer();
            if exist('psyrat_defineversion', 'file')
                try
                    meta.psyrat_version = psyrat_defineversion();
                catch
                    meta.psyrat_version = 'unavailable';
                end
            else
                meta.psyrat_version = 'unavailable';
            end
            for f = {'analysis', 'family', 'engine', 'gammascale'}
                if isstruct(rel) && isfield(rel, f{1})
                    meta.(['rel_' f{1}]) = rel.(f{1});
                end
            end
            if isstruct(rel) && isfield(rel, 'out') && isstruct(rel.out) ...
                    && isfield(rel.out, 'conv') && isstruct(rel.out.conv)
                meta.rel_conv_fields = fieldnames(rel.out.conv)';
                if isfield(rel.out.conv, 'converged')
                    meta.rel_conv_converged = rel.out.conv.converged;
                end
            end
            meta.n_rows = height(dataTbl);
            meta.data_variables = dataTbl.Properties.VariableNames;
            meta.cmdstan_source_dir = cmdstanOutDir;
            meta.caller = callerMeta;

            % ---- REL / truth / data as MAT ----
            matFile = fullfile(dest, [tag '_rel.mat']);
            try
                save(matFile, 'rel', 'truth', 'dataTbl', 'meta', '-v7');
            catch errV7
                % v7 has a 2 GB per-variable cap; fall back rather than lose the run
                warning('PsyRAT2CRecoveryCapture:v73fallback', ...
                    'MAT v7 save failed (%s); retrying with -v7.3.', errV7.message);
                save(matFile, 'rel', 'truth', 'dataTbl', 'meta', '-v7.3');
            end

            % ---- data as CSV, truth as JSON (external readers) ----
            writetable(dataTbl, fullfile(dest, 'data.csv'));
            PsyRAT2CRecoveryCapture.writeText(fullfile(dest, 'truth.json'), ...
                jsonencode(truth, 'PrettyPrint', true));

            % ---- raw CmdStan output (filtered copy) and success criteria ----
            % REL/data/truth are already on disk above, so every failure below
            % still leaves them; the raw evidence is what decides success.
            failures = {};
            copied = {}; skipped = {}; sizes = struct([]); chainEvidence = struct([]);
            meta.cmdstan_source_existed = isfolder(cmdstanOutDir);
            if ~meta.cmdstan_source_existed
                failures{end+1} = sprintf('source CmdStan directory "%s" does not exist', cmdstanOutDir);
            else
                [copied, skipped, sizes, copyErrs] = PsyRAT2CRecoveryCapture.copyFiltered( ...
                    cmdstanOutDir, fullfile(dest, 'cmdstan_output'));
                failures = [failures, copyErrs];
                % chain CSVs (CmdStan convention output-*.csv) must exist, carry a
                % header + at least one sample row, and match the source bytes
                isChain = ~cellfun(@isempty, regexp(copied, '(^|[\\/])output-[^\\/]*\.csv$', 'once'));
                chainFiles = copied(isChain);
                if isempty(chainFiles)
                    failures{end+1} = 'no chain CSV (output-*.csv) was copied from the source directory';
                end
                for i = 1:numel(chainFiles)
                    src = fullfile(cmdstanOutDir, chainFiles{i});
                    dst = fullfile(dest, 'cmdstan_output', chainFiles{i});
                    ev = PsyRAT2CRecoveryCapture.chainEvidence(src, dst);
                    ev.file = chainFiles{i};
                    if isempty(chainEvidence), chainEvidence = ev; else, chainEvidence(end+1) = ev; end %#ok<AGROW>
                    if ~ev.header_ok
                        failures{end+1} = sprintf('%s: no CmdStan header line beginning lp__', chainFiles{i}); %#ok<AGROW>
                    end
                    if ev.n_sample_rows < 1
                        failures{end+1} = sprintf('%s: no sample rows after the header', chainFiles{i}); %#ok<AGROW>
                    end
                    if ~ev.bytes_identical
                        failures{end+1} = sprintf('%s: copied bytes differ from the source', chainFiles{i}); %#ok<AGROW>
                    end
                end
                info.chainFiles = chainFiles;
            end
            meta.cmdstan_copied  = copied;
            meta.cmdstan_skipped = skipped;
            meta.cmdstan_file_sizes = sizes;
            meta.chain_csv_evidence = chainEvidence;
            meta.capture_failures = failures;
            if isempty(failures)
                meta.capture_status = 'ok';
            else
                meta.capture_status = 'failed';
            end

            PsyRAT2CRecoveryCapture.writeText(fullfile(dest, 'meta.json'), ...
                jsonencode(meta, 'PrettyPrint', true));
            PsyRAT2CRecoveryCapture.writeText(fullfile(dest, 'README.txt'), ...
                PsyRAT2CRecoveryCapture.readmeText());

            info.dest = dest;
            info.matFile = matFile;
            info.copied = copied;
            info.skipped = skipped;
            if ~isempty(failures)
                error('PsyRAT2CRecoveryCapture:rawCaptureIncomplete', ...
                    ['Raw CmdStan capture for "%s" is INCOMPLETE (REL/data/truth were written to "%s"; ', ...
                    'meta.json records capture_status=failed). Reasons: %s'], ...
                    tag, dest, strjoin(failures, '; '));
            end
            info.captured = true;
            fprintf('[2C recovery capture] %s -> %s (%d CmdStan files copied, %d skipped, %d chain CSVs verified)\n', ...
                tag, dest, numel(copied), numel(skipped), numel(info.chainFiles));
        end

    end

    methods (Static, Access = private)

        function [copied, skipped, sizes, errs] = copyFiltered(srcDir, dstDir)
            % Recursive filtered copy: keep CSV/Stan/data/metadata/log files,
            % skip executables and build intermediates. Records both lists
            % relative to srcDir, the source/copy byte sizes of every copied
            % file, and a list of copy/size failures (the caller decides the
            % outcome; nothing is silently dropped). Never modifies srcDir.
            % The caller has already verified that srcDir exists.
            copied = {};
            skipped = {};
            sizes = struct([]);
            errs = {};
            listing = dir(fullfile(srcDir, '**', '*'));
            listing = listing(~[listing.isdir]);
            if isempty(listing)
                errs{end+1} = sprintf('source CmdStan directory "%s" contains no files', srcDir);
                return;
            end
            keepExt = PsyRAT2CRecoveryCapture.KeepExtensions;
            for i = 1:numel(listing)
                srcFile = fullfile(listing(i).folder, listing(i).name);
                relPath = PsyRAT2CRecoveryCapture.relativeTo(srcDir, srcFile);
                [~, ~, ext] = fileparts(listing(i).name);
                if any(strcmpi(ext, keepExt))
                    dstFile = fullfile(dstDir, relPath);
                    dstFolder = fileparts(dstFile);
                    if ~isfolder(dstFolder)
                        mkdir(dstFolder);
                    end
                    [okCopy, msgCopy] = copyfile(srcFile, dstFile);
                    if ~okCopy
                        errs{end+1} = sprintf('could not copy "%s": %s', relPath, msgCopy); %#ok<AGROW>
                        continue;
                    end
                    d = dir(dstFile);
                    s = struct('file', relPath, 'source_bytes', listing(i).bytes, ...
                        'copy_bytes', d(1).bytes);
                    if isempty(sizes), sizes = s; else, sizes(end+1) = s; end %#ok<AGROW>
                    if d(1).bytes ~= listing(i).bytes
                        errs{end+1} = sprintf('"%s": copy is %d bytes, source is %d bytes', ...
                            relPath, d(1).bytes, listing(i).bytes); %#ok<AGROW>
                    end
                    copied{end+1} = relPath; %#ok<AGROW>
                else
                    skipped{end+1} = relPath; %#ok<AGROW>
                end
            end
        end

        function ev = chainEvidence(srcFile, dstFile)
            % Header/sample evidence for one CmdStan chain CSV plus a byte
            % comparison against its source. CmdStan CSVs start with '#'
            % comment lines, then a header beginning 'lp__', then sample rows
            % (further '#' lines, e.g. adaptation/timing, are skipped).
            % Numeric content is deliberately NOT validated.
            ev = struct('file', '', 'header_ok', false, 'n_columns', 0, ...
                'n_sample_rows', 0, 'bytes_identical', false, 'source_bytes', 0, 'copy_bytes', 0);
            srcBytes = PsyRAT2CRecoveryCapture.readBytes(srcFile);
            dstBytes = PsyRAT2CRecoveryCapture.readBytes(dstFile);
            ev.source_bytes = numel(srcBytes);
            ev.copy_bytes = numel(dstBytes);
            ev.bytes_identical = isequal(srcBytes, dstBytes);
            lines = regexp(char(dstBytes(:)'), '\r?\n', 'split');
            headerSeen = false;
            for i = 1:numel(lines)
                ln = strtrim(lines{i});
                if isempty(ln) || ln(1) == '#'
                    continue;
                end
                if ~headerSeen
                    ev.header_ok = strncmp(ln, 'lp__', 4);
                    if ~ev.header_ok
                        return;   % first non-comment line is not a CmdStan header
                    end
                    ev.n_columns = numel(strsplit(ln, ','));
                    headerSeen = true;
                else
                    ev.n_sample_rows = ev.n_sample_rows + 1;
                end
            end
        end

        function b = readBytes(fileName)
            fid = fopen(fileName, 'r');
            if fid < 0
                b = uint8([]);
                return;
            end
            cleaner = onCleanup(@() fclose(fid));
            b = fread(fid, Inf, '*uint8');
        end

        function p = normalizePath(p)
            % Absolute, separator-normalized path without a trailing separator.
            p = char(p);
            if isempty(regexp(p, '^(/|[A-Za-z]:[\\/]|\\\\)', 'once'))
                p = fullfile(pwd, p);
            end
            p = regexprep(p, '[\\/]+', filesep);
            p = regexprep(p, '[\\/]+$', '');
            if ispc, p = lower(p); end
        end

        function tf = pathsOverlap(a, b)
            % True when the two normalized paths are equal or one contains the other.
            tf = strcmp(a, b) || strncmp([a filesep], [b filesep], numel(a) + 1) ...
                || strncmp([b filesep], [a filesep], numel(b) + 1);
        end

        function rel = relativeTo(baseDir, absPath)
            baseDir = char(baseDir); absPath = char(absPath);
            if strncmp(absPath, baseDir, numel(baseDir))
                rel = absPath(numel(baseDir)+1:end);
                rel = regexprep(rel, '^[\\/]+', '');
            else
                [~, n, e] = fileparts(absPath);
                rel = [n e];
            end
        end

        function writeText(fileName, txt)
            fid = fopen(fileName, 'w');
            if fid < 0
                error('PsyRAT2CRecoveryCapture:writeFailed', 'Cannot open "%s" for writing.', fileName);
            end
            cleaner = onCleanup(@() fclose(fid));
            fwrite(fid, txt, 'char');
        end

        function txt = readmeText()
            txt = sprintf([ ...
                'PsyRAT 2C recovery capture\n', ...
                '==========================\n', ...
                'Written by tests/helpers/PsyRAT2CRecoveryCapture.m when the environment\n', ...
                'variable PSYRAT_2C_RECOVERY_CAPTURE_DIR was set during a live recovery test.\n', ...
                'The capture is taken immediately after psyrat_computevarcomp returns and\n', ...
                'BEFORE the test''s convergence assumptions, so a skipped (Incomplete) run\n', ...
                'still leaves its evidence here.\n\n', ...
                'Contents:\n', ...
                '  <tag>_rel.mat     rel (REL struct), truth (planted), dataTbl (input), meta\n', ...
                '  data.csv          the simulated long table\n', ...
                '  truth.json        planted truth (numeric arrays as JSON lists)\n', ...
                '  meta.json         run/route metadata, the CmdStan copy manifest with sizes,\n', ...
                '                    per-chain-CSV header/sample-row evidence, and capture_status\n', ...
                '                    (''ok'' only if chain CSVs were copied and verified; ''failed''\n', ...
                '                    with reasons otherwise; a failed capture still keeps REL/data/truth)\n', ...
                '  cmdstan_output/   filtered copy of the raw CmdStan output directory\n', ...
                '                    (CSV/Stan/data/metadata/log kept; executables and\n', ...
                '                    build intermediates skipped and listed in meta.json)\n\n', ...
                'This directory is evidence, not a claim: convergence, coverage and\n', ...
                'recovery are adjudicated from it separately. Do not delete or overwrite.\n']);
        end

    end
end
