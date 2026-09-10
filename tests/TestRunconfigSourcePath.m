classdef TestRunconfigSourcePath < PsyRATTestBase
    %Audit 2026-09-05: the sidecar's input-file hash was computed from
    %raw.filename, the bare file name both producers store for display, so
    %exist() resolved it only when the file sat in the current folder or on
    %the MATLAB path, and could resolve a same-named file elsewhere. Both
    %producers (psyrat_run, psyrat_startproc) now also record raw.sourcepath,
    %the full path, and psyrat_write_runconfig hashes from it and records it
    %as data.source_path. A struct without sourcepath keeps the old
    %best-effort behavior (the existing TestCliEntryPoint pins pass a full
    %path in raw.filename and are unaffected).
    %
    %RED against the pre-fix code: testHashesFromTheSourcePath (the file is
    %in a temp folder that is neither pwd nor on the path, so the old code
    %left file_sha256 empty and had no source_path key).

    methods (Test)

        function testHashesFromTheSourcePath(testCase)
            tdir = tempname; mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
            csvpath = fullfile(tdir, 'scores.csv');
            tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.2; 0.9; 1.0], ...
                'VariableNames', {'id','meas'});
            writetable(tbl, csvpath);
            %guard the premise: the bare name must NOT resolve from here
            testCase.assertNotEqual(pwd, tdir);
            testCase.assertEqual(exist('scores.csv', 'file'), 0, ...
                'a file named scores.csv must not be reachable by bare name, or this proves nothing');

            data = localData(tdir, tbl);
            data.raw.filename = 'scores.csv';          % bare, as both producers store it
            data.raw.sourcepath = csvpath;             % full, as both producers now add
            sidecar = psyrat_write_runconfig(data, psyrat_defaults());
            cfg = jsondecode(fileread(sidecar));

            testCase.verifyEqual(cfg.data.source_file, 'scores.csv');
            testCase.verifyEqual(cfg.data.source_path, csvpath);
            testCase.verifyTrue(localIsSha256(cfg.data.file_sha256), ...
                'the file hash must be computed from the full path');
            testCase.verifyEqual(cfg.data.file_sha256, psyrat_hash_file(csvpath));
        end

        function testWithoutSourcePathKeepsTheOldBehavior(testCase)
            tdir = tempname; mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
            tbl = table({'s1';'s2'}, [1.1; 0.9], 'VariableNames', {'id','meas'});
            data = localData(tdir, tbl);
            data.raw.filename = 'not_reachable_by_bare_name.csv';
            sidecar = psyrat_write_runconfig(data, psyrat_defaults());
            cfg = jsondecode(fileread(sidecar));
            testCase.verifyEqual(cfg.data.source_file, 'not_reachable_by_bare_name.csv');
            testCase.verifyEqual(cfg.data.source_path, '');
            testCase.verifyEqual(cfg.data.file_sha256, '');
        end

        function testPsyratRunRecordsTheSourcePath(testCase)
            %the producer side: psyrat_run must carry the full path into raw
            src = fileread(which('psyrat_run'));
            src = regexprep(src, '(?m)%[^\n]*', '');
            testCase.verifyTrue(contains(src, 'rawsourcepath = rawfile;'), ...
                'psyrat_run must record the full input path for the sidecar');
            testCase.verifyTrue(contains(src, 'psyrat_data.raw.sourcepath = rawsourcepath;'));
            gui = fileread(which('psyrat_startproc'));
            gui = regexprep(gui, '(?m)%[^\n]*', '');
            testCase.verifyTrue(contains(gui, 'psyrat_data.raw.sourcepath = fullfile(pathpart,filepart);'), ...
                'the GUI must record the full input path for the sidecar');
        end

    end
end

function data = localData(tdir, tbl)
data = struct();
data.proc.savepath = tdir;
data.proc.savename = 'demo.psyrat';
data.proc.idheader = 'id';
data.proc.measheader = 'meas';
data.proc.eventheader = '';
data.proc.groupheader = '';
data.proc.timeheader = '';
data.proc.whichevents = '';
data.proc.whichgroups = '';
data.proc.whichtimes = '';
data.rel.analysis = 'ic';
data.raw.data = tbl;
end

function tf = localIsSha256(s)
tf = ischar(s) && numel(s) == 64 && all(ismember(lower(s), '0123456789abcdef'));
end
