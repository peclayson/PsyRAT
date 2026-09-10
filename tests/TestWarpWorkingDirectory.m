classdef TestWarpWorkingDirectory < matlab.unittest.TestCase
    %TESTWARPWORKINGDIRECTORY Pins psyrat_computevarcompwarp's working-directory
    % contract: the run must end in the directory it started in, whether
    % estimation succeeds or throws.
    %
    % WHY THIS EXISTS
    % The warp cd's into a per-run Stan temp directory before estimation
    % (origdir = cd(...)) and cd's back roughly three hundred lines later. Nothing
    % spanned the psyrat_computevarcomp call in between, so when estimation threw
    % - a rejected option, a preflight failure, a CmdStan error - the restore
    % never ran and MATLAB was left INSIDE the temp directory.
    %
    % The failure is delayed and misleading, which is what makes it worth a test.
    % The temp directory is removed after the run, so the session is left pointing
    % at a path that no longer exists. Nothing goes wrong until the NEXT run,
    % which records that dangling path as its own origdir and dies in cd() with
    % 'MATLAB:cd:NonExistentFolder' naming a temp directory from an unrelated
    % earlier analysis. A user would have no way to connect that to the run that
    % actually failed.
    %
    % The trigger below (dispersion=2 outside the group-level concurrent Gamma
    % difference design) is chosen only because it is the cheapest engine-side
    % error that needs no CmdStan and no special data. Any error raised by
    % psyrat_computevarcomp exercises the same path; this is a test of the warp's
    % cleanup contract, not of the dispersion option.

    methods (Test)

        function testWorkingDirectoryRestoredWhenEstimationErrors(testCase)
            tdir = localTempDir(testCase);
            tbl = localTable();

            before = pwd;
            testCase.verifyError(@() psyrat_run('data', tbl, 'dispersion', 2, ...
                'savepath', tdir, 'savename', 'x.psyrat'), 'varargin:dispersion');

            testCase.verifyEqual(pwd, before, ...
                ['The warp must restore the working directory when estimation ' ...
                'errors. Leaving the session in the run''s temp directory makes ' ...
                'the NEXT run fail in cd() with an unrelated path.']);
        end

        % ON WHY THERE IS ONLY ONE TEST HERE.
        % Two further tests were written and both were DELETED after they passed
        % against the unfixed code, i.e. they proved nothing:
        %   - asserting pwd still EXISTS right after the failed run passes,
        %     because the warp aborts before its own rmdir, so the temp directory
        %     is still there at that instant;
        %   - a failed run, then removing its tree, then a SECOND failed run also
        %     passes, because run 2 aborts before reaching its own cd(origdir).
        % The delayed failure needs a later run that SUCCEEDS far enough to reach
        % the restore. That path is already covered, in the full-suite ordering,
        % by TestInstallationAndWorkflowValidation - which is precisely where this
        % bug surfaced, as 'MATLAB:cd:NonExistentFolder' naming a temp directory
        % from an unrelated earlier test. Reproducing it here would mean driving a
        % complete estimation, and the single assertion above already pins the
        % root cause directly. A test that cannot fail is worse than no test.

    end
end

function tbl = localTable()
%Minimal valid one-facet long table; the run never reaches estimation.
tbl = table({'s1';'s1';'s2';'s2'}, [1.1; 1.4; 0.9; 1.0], ...
    'VariableNames', {'id','meas'});
end

function d = localTempDir(testCase)
d = tempname;
mkdir(d);
testCase.addTeardown(@() localRemoveDir(d));
end

function localRemoveDir(d)
if exist(d, 'dir') == 7
    try, rmdir(d, 's'); catch, end
end
end
