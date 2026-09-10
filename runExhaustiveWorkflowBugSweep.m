function [allPassed, report] = runExhaustiveWorkflowBugSweep(varargin)
%RUNEXHAUSTIVEWORKFLOWBUGSWEEP Exhaustive end-to-end workflow bug sweep.
%
% This script runs CmdStan estimation from test_data and then executes
% exhaustive relsummary/output combinations for each analysis family.
% Difference-score workflows sweep both absolute and relative decisions
% via diffgcoeff values [1 2].
%
% Usage:
%   [allPassed, report] = runExhaustiveWorkflowBugSweep();
%
% Optional name-value inputs:
%   'Chains'            - Stan chains (default: 1)
%   'Iter'              - Stan iterations (default: 30)
%   'Verbose'           - print progress (default: true)
%   'IncludePlots'      - run plotting outputs (default: true)
%   'ScenarioFilter'    - scenario names to run (default: all)
%   'CmdStanOutputDir'  - output dir for CmdStan artifacts
%                         (default: tests/integration/cmdstan_artifacts)
%   'SaveReport'        - write csv+mat report files (default: true)
%   'ReportPrefix'      - report file prefix (default: workflow_bug_sweep)
%   'DiffWpCov'         - difference-score residual covariance mode
%                         1 = force wp cov = 0 (default), 2 = estimate/use wp cov
%
% Report fields:
%   report.summary
%   report.results (table with per-step pass/fail details)
%   report.options

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

parser = inputParser;
parser.addParameter('Chains', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('Iter', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('IncludePlots', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('ScenarioFilter', {}, @(x) iscell(x) || isstring(x) || ischar(x));
parser.addParameter('CmdStanOutputDir', '', @(x) ischar(x) || isstring(x));
parser.addParameter('SaveReport', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('ReportPrefix', 'workflow_bug_sweep', @(x) ischar(x) || isstring(x));
parser.addParameter('DiffWpCov', 1, @(x) isnumeric(x) && isscalar(x) && any(x == [1 2]));
parser.parse(varargin{:});
opts = parser.Results;

opts.Verbose = logical(opts.Verbose);
opts.IncludePlots = logical(opts.IncludePlots);
opts.SaveReport = logical(opts.SaveReport);
opts.ScenarioFilter = cellstr(string(opts.ScenarioFilter));
opts.ReportPrefix = char(string(opts.ReportPrefix));

projectRoot = fileparts(mfilename('fullpath'));
bundledStan = fullfile(projectRoot, 'bundled_dependents', 'MatlabStan-2.15.1.0');
if exist(bundledStan, 'dir')
    addpath(genpath(bundledStan), '-begin');
    rehash;
    % Force MatlabStan class reload so stale StanModel definitions from
    % prior paths cannot persist across repeated sweeps.
    clear StanModel StanFit;
    rehash;
    stanModelPath = which('StanModel.m');
    if isempty(stanModelPath) || ~contains(stanModelPath, bundledStan)
        warning('workflow:stanmodelpath', ...
            'StanModel.m is not resolving to bundled MatlabStan: %s', stanModelPath);
    end
end
addpath(genpath(fullfile(projectRoot, 'subroutines')), '-begin');
addpath(fullfile(projectRoot, 'tests', 'helpers'), '-begin');

if isempty(strtrim(getenv('PSYRAT_TEST_CMDSTAN_PATH')))
    %G51: consult the shared resolver (env vars > toolbox-managed home globs
    %> MATLAB path, validity-checked; see psyrat_locate_cmdstan/G43) instead
    %of a path-only which() probe, so an off-path CmdStan the startup gate
    %accepts is found here too. The resolver already enforces the
    %makefile+bin criterion the old inline check duplicated. The subroutines
    %path was prepended above, so the resolver is reachable here.
    cmdstanHome = psyrat_locate_cmdstan();
    if ~isempty(cmdstanHome)
        setenv('PSYRAT_TEST_CMDSTAN_PATH', cmdstanHome);
    end
end

if isempty(which('PsyRATCmdStanIntegrationFactory'))
    error('workflow:missingfactory', ...
        'PsyRATCmdStanIntegrationFactory is unavailable. Ensure tests/helpers is on path.');
end

if isempty(which('processManager'))
    warning('workflow:missingprocessmanager', ...
        ['processManager is unavailable. Set PSYRAT_TEST_MPM_PATH or add ', ...
        'MatlabProcessManager to path.']);
end
if isempty(which('stan'))
    warning('workflow:missingstan', ...
        ['stan is unavailable. Set PSYRAT_TEST_MATLABSTAN_PATH or add ', ...
        'MatlabStan to path.']);
end

if strlength(string(opts.CmdStanOutputDir)) == 0
    outDir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
else
    outDir = char(string(opts.CmdStanOutputDir));
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

scenarios = PsyRATCmdStanIntegrationFactory.scenarios();
scenarioNames = fieldnames(scenarios);
if ~isempty(opts.ScenarioFilter)
    keep = ismember(scenarioNames, opts.ScenarioFilter);
    scenarioNames = scenarioNames(keep);
end

if isempty(scenarioNames)
    error('workflow:noscenarios', 'No scenarios selected to run.');
end

if opts.Verbose
    fprintf('\nRunning exhaustive workflow bug sweep\n');
    fprintf('Scenarios: %d\n', numel(scenarioNames));
    fprintf('CmdStan output dir: %s\n\n', outDir);
end

allRows = {};
for s = 1:numel(scenarioNames)
    scenarioName = scenarioNames{s};
    cfg = scenarios.(scenarioName);
    if opts.Verbose
        fprintf('--- Scenario %d/%d: %s ---\n', s, numel(scenarioNames), scenarioName);
    end
    
    tStart = tic;
    [okEst, msgEst, psyrat_data] = localRunEstimation(projectRoot, cfg, opts, outDir);
    allRows(end+1,:) = localRow('estimate', scenarioName, cfg.expectedAnalysis, ...
        'estimate_cmdstan', okEst, toc(tStart), msgEst); %#ok<AGROW>
    
    if ~okEst
        if opts.Verbose
            fprintf('  Estimation failed: %s\n\n', msgEst);
        end
        continue;
    end
    
    runs = psyrat_workflow_runs(psyrat_data.rel.analysis);
    for r = 1:numel(runs)
        runCfg = runs{r};
        [rowsOut, passCount, failCount] = localRunWorkflow(psyrat_data, scenarioName, ...
            runCfg, opts);
        allRows = [allRows; rowsOut]; %#ok<AGROW>
        
        if opts.Verbose
            fprintf('  %s: %d passed, %d failed\n', runCfg.analysis, passCount, failCount);
        end
    end
    
    if opts.Verbose
        fprintf('\n');
    end
end

results = localRowsToTable(allRows);
nFail = sum(~results.passed);
nPass = sum(results.passed);
allPassed = nFail == 0;

summary = struct();
summary.total = height(results);
summary.passed = nPass;
summary.failed = nFail;
summary.scenarios = scenarioNames;

report = struct();
report.summary = summary;
report.results = results;
report.options = opts;

if opts.SaveReport
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    csvFile = fullfile(projectRoot, sprintf('%s_%s.csv', opts.ReportPrefix, ts));
    matFile = fullfile(projectRoot, sprintf('%s_%s.mat', opts.ReportPrefix, ts));
    writetable(results, csvFile);
    save(matFile, 'report');
    if opts.Verbose
        fprintf('Saved report:\n  %s\n  %s\n', csvFile, matFile);
    end
end

if opts.Verbose
    fprintf('\nSweep complete. Passed: %d, Failed: %d, Total: %d\n', ...
        summary.passed, summary.failed, summary.total);
end

end

function [ok, msg, psyrat_data] = localRunEstimation(projectRoot, cfg, opts, outDir)
ok = false;
msg = '';
psyrat_data = struct();
try
    dataTbl = PsyRATCmdStanIntegrationFactory.makeScenarioTable(projectRoot, cfg);
    estimatorFn = 'psyrat_computevarcomp';
    if isempty(which(estimatorFn))
        estimatorFn = 'psyrat_computerel';
    end
    rel = feval(estimatorFn, ...
        'data', dataTbl, ...
        'chains', opts.Chains, ...
        'iter', opts.Iter, ...
        'verbose', 1, ...
        'showgui', 1, ...
        'cmdstanoutdir', outDir, ...
        'sserrvar', cfg.sserrvar, ...
        'diffest', cfg.diffest, ...
        'diffwpcov', opts.DiffWpCov);
    psyrat_data.rel = rel;
    ok = true;
catch ME
    msg = localErrMsg(ME);
end
end

function [rowsOut, passCount, failCount] = localRunWorkflow(psyrat_data, scenarioName, runCfg, opts)
rowsOut = {};
passCount = 0;
failCount = 0;

switch runCfg.analysis
    case {'sing', 'sing_sserr', 'sing_diff', 'sing_diff_sserr'}
        if isfield(runCfg,'diffgcoeff') && ~isempty(runCfg.diffgcoeff)
            diffgcoeffValues = runCfg.diffgcoeff(:)';
        else
            diffgcoeffValues = 1;
        end
        
        for meascutoff = 1:3
            for depcentmeas = 1:2
                for diffgcoeff = diffgcoeffValues
                    if contains(runCfg.analysis,'sing_diff')
                        comboName = sprintf('meascutoff=%d;depcent=%d;diffg=%d', ...
                            meascutoff, depcentmeas, diffgcoeff);
                        tStart = tic;
                        [ok, msg, pdat] = localRunRelsummary( ...
                            'psyrat_data', psyrat_data, ...
                            'analysis', runCfg.analysis, ...
                            'depcutoff', runCfg.depcutoff, ...
                            'meascutoff', meascutoff, ...
                            'depcentmeas', depcentmeas, ...
                            'diffgcoeff', diffgcoeff, ...
                            'CI', 0.95);
                    else
                        comboName = sprintf('meascutoff=%d;depcent=%d', meascutoff, depcentmeas);
                        tStart = tic;
                        [ok, msg, pdat] = localRunRelsummary( ...
                            'psyrat_data', psyrat_data, ...
                            'analysis', runCfg.analysis, ...
                            'depcutoff', runCfg.depcutoff, ...
                            'meascutoff', meascutoff, ...
                            'depcentmeas', depcentmeas, ...
                            'CI', 0.95);
                    end
                    
                    rowsOut(end+1,:) = localRow('relsummary', scenarioName, runCfg.analysis, ...
                        comboName, ok, toc(tStart), msg); %#ok<AGROW>
                    [passCount, failCount] = localCount(ok, passCount, failCount);
                    if ~ok
                        continue;
                    end
                    if contains(msg, 'nogooddata=1')
                        continue;
                    end
                    
                    [rowsOut, passCount, failCount] = localRunOutputs( ...
                        rowsOut, passCount, failCount, pdat, scenarioName, runCfg.analysis, ...
                        comboName, opts);
                end
            end
        end
        
    case 'trt'
        for gcoeff = 1:2
            for reltype = 1:3
                for meascutoff = 1:3
                    for relcentmeas = 1:2
                        comboName = sprintf('g=%d;reltype=%d;meascutoff=%d;relcent=%d', ...
                            gcoeff, reltype, meascutoff, relcentmeas);
                        tStart = tic;
                        [ok, msg, pdat] = localRunRelsummary( ...
                            'psyrat_data', psyrat_data, ...
                            'analysis', 'trt', ...
                            'gcoeff', gcoeff, ...
                            'reltype', reltype, ...
                            'relcutoff', runCfg.relcutoff, ...
                            'meascutoff', meascutoff, ...
                            'relcentmeas', relcentmeas, ...
                            'CI', 0.95);
                        rowsOut(end+1,:) = localRow('relsummary', scenarioName, runCfg.analysis, ...
                            comboName, ok, toc(tStart), msg); %#ok<AGROW>
                        [passCount, failCount] = localCount(ok, passCount, failCount);
                        if ~ok
                            continue;
                        end
                        if contains(msg, 'nogooddata=1')
                            continue;
                        end
                        
                        [rowsOut, passCount, failCount] = localRunOutputs( ...
                            rowsOut, passCount, failCount, pdat, scenarioName, runCfg.analysis, ...
                            comboName, opts);
                    end
                end
            end
        end

    case 'trt_sserr'
        %subject-level (per-participant) crossed test-retest (trt_sserrvar):
        %test-retest coefficient selection (gcoeff/reltype/meascutoff) with the
        %subject-level participant-inclusion cutoff (depcutoff, not relcutoff).
        for gcoeff = 1:2
            for reltype = 1:3
                for meascutoff = 1:3
                    comboName = sprintf('g=%d;reltype=%d;meascutoff=%d', ...
                        gcoeff, reltype, meascutoff);
                    tStart = tic;
                    [ok, msg, pdat] = localRunRelsummary( ...
                        'psyrat_data', psyrat_data, ...
                        'analysis', 'trt_sserr', ...
                        'gcoeff', gcoeff, ...
                        'reltype', reltype, ...
                        'depcutoff', runCfg.depcutoff, ...
                        'meascutoff', meascutoff, ...
                        'CI', 0.95);
                    rowsOut(end+1,:) = localRow('relsummary', scenarioName, runCfg.analysis, ...
                        comboName, ok, toc(tStart), msg); %#ok<AGROW>
                    [passCount, failCount] = localCount(ok, passCount, failCount);
                    if ~ok
                        continue;
                    end
                    if contains(msg, 'nogooddata=1')
                        continue;
                    end
                    [rowsOut, passCount, failCount] = localRunOutputs( ...
                        rowsOut, passCount, failCount, pdat, scenarioName, runCfg.analysis, ...
                        comboName, opts);
                end
            end
        end

    case 'trt_diff'
        if isfield(runCfg,'diffgcoeff') && ~isempty(runCfg.diffgcoeff)
            diffgcoeffValues = runCfg.diffgcoeff(:)';
        else
            diffgcoeffValues = 1;
        end
        for gcoeff = 1:2
            for reltype = 1:3
                for meascutoff = 1:3
                    for relcentmeas = 1:2
                        for diffgcoeff = diffgcoeffValues
                            comboName = sprintf(...
                                'g=%d;reltype=%d;meascutoff=%d;relcent=%d;diffg=%d', ...
                                gcoeff, reltype, meascutoff, relcentmeas, diffgcoeff);
                            tStart = tic;
                            [ok, msg, pdat] = localRunRelsummary( ...
                                'psyrat_data', psyrat_data, ...
                                'analysis', 'trt_diff', ...
                                'gcoeff', gcoeff, ...
                                'diffgcoeff', diffgcoeff, ...
                                'reltype', reltype, ...
                                'relcutoff', runCfg.relcutoff, ...
                                'meascutoff', meascutoff, ...
                                'relcentmeas', relcentmeas, ...
                                'CI', 0.95);
                            rowsOut(end+1,:) = localRow('relsummary', scenarioName, runCfg.analysis, ...
                                comboName, ok, toc(tStart), msg); %#ok<AGROW>
                            [passCount, failCount] = localCount(ok, passCount, failCount);
                            if ~ok
                                continue;
                            end
                            if contains(msg, 'nogooddata=1')
                                continue;
                            end

                            [rowsOut, passCount, failCount] = localRunOutputs( ...
                                rowsOut, passCount, failCount, pdat, scenarioName, runCfg.analysis, ...
                                comboName, opts);
                        end
                    end
                end
            end
        end

    otherwise
        rowsOut(end+1,:) = localRow('workflow', scenarioName, runCfg.analysis, ...
            'unsupported_analysis', false, 0, ...
            sprintf('Unsupported analysis type: %s', runCfg.analysis)); %#ok<AGROW>
        failCount = failCount + 1;
end
end

function [rowsOut, passCount, failCount] = localRunOutputs( ...
    rowsOut, passCount, failCount, pdat, scenarioName, analysis, comboName, opts)

isTRT = strcmp(analysis, 'trt');
isSS = isfield(pdat, 'rel') && isfield(pdat.rel, 'analysis') && ...
    any(strcmpi(pdat.rel.analysis, {'ic_sserrvar', 'trt_sserrvar', 'ic_diff_sserrvar'}));

steps = {};
if any(strcmp(analysis, {'sing','sing_diff'}))
    steps{end+1} = {'table_depcutoff', @() psyrat_depcutofft('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    steps{end+1} = {'table_depoverall', @() psyrat_depoverallt('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    steps{end+1} = {'table_variance', @() psyrat_variancet('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
elseif any(strcmp(analysis, {'sing_sserr','sing_diff_sserr'}))
    steps{end+1} = {'table_depoverall', @() psyrat_depoverallt('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    steps{end+1} = {'table_variance', @() psyrat_variancet('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    if strcmp(analysis, 'sing_diff_sserr')
        %G60: the real viewer opens the trial-cutoff table on the DIFF
        %subject-level route (psyrat_relfigures' sing_diff_sserr case; the
        %diff row shows the observed-design coefficient with '---' for the
        %count), and this was the one rendered surface the sweep never
        %exercised. sing_sserr stays out: the viewer does not call it there.
        steps{end+1} = {'table_depcutoff', @() psyrat_depcutofft('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    end
elseif strcmp(analysis, 'trt_sserr')
    %subject-level test-retest: the group-level SD/SEM/ICC table (the
    %per-subject coefficients ride on the isSS ss-plot steps below).
    steps{end+1} = {'table_variance', @() psyrat_variancet('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
elseif isTRT
    steps{end+1} = {'table_trt_cutoff', @() psyrat_trt_relcutofft('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    steps{end+1} = {'table_trt_overall', @() psyrat_trt_reloverallt('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
    steps{end+1} = {'table_variance', @() psyrat_variancet('psyrat_data', pdat, 'gui', 0)}; %#ok<AGROW>
end

if opts.IncludePlots
    if any(strcmp(analysis, {'sing','sing_diff'}))
        steps{end+1} = {'plot_dep_trials', @() psyrat_depvtrialsplot( ...
            'psyrat_data', pdat, 'trials', [1 30], 'depline', 1, 'depcutoff', 0.5)}; %#ok<AGROW>
        steps{end+1} = {'plot_icc', @() psyrat_ptintervalplot('psyrat_data', pdat, 'stat', 'icc')}; %#ok<AGROW>
        steps{end+1} = {'plot_bet', @() psyrat_ptintervalplot('psyrat_data', pdat, 'stat', 'bet')}; %#ok<AGROW>
    elseif isTRT
        steps{end+1} = {'plot_trt_trials', @() psyrat_trt_relvtrialsplot( ...
            'psyrat_data', pdat, 'trials', [1 30], 'relline', 1, 'relcutoff', 0.5)}; %#ok<AGROW>
        steps{end+1} = {'plot_icc', @() psyrat_ptintervalplot('psyrat_data', pdat, 'stat', 'icc')}; %#ok<AGROW>
        steps{end+1} = {'plot_bet', @() psyrat_ptintervalplot('psyrat_data', pdat, 'stat', 'bet')}; %#ok<AGROW>
    end
    
    if isSS
        steps{end+1} = {'plot_ss_dep', @() psyrat_ssrelplot('psyrat_data', pdat, 'stat', 'dep')}; %#ok<AGROW>
        steps{end+1} = {'plot_ss_icc', @() psyrat_ssrelplot('psyrat_data', pdat, 'stat', 'icc')}; %#ok<AGROW>
    end
end

if strcmp(analysis,'sing')
    steps{end+1} = {'criterion_outputs', @() psyrat_criterionfigures( ...
        'psyrat_data', pdat, 'analysis', 'sing')}; %#ok<AGROW>
elseif isTRT
    steps{end+1} = {'criterion_outputs', @() psyrat_criterionfigures( ...
        'psyrat_data', pdat, 'analysis', 'trt')}; %#ok<AGROW>
end

for i = 1:numel(steps)
    stepName = steps{i}{1};
    fn = steps{i}{2};
    tStart = tic;
    [ok, msg] = localCall(fn);
    close all;
    rowsOut(end+1,:) = localRow(stepName, scenarioName, analysis, comboName, ...
        ok, toc(tStart), msg); %#ok<AGROW>
    [passCount, failCount] = localCount(ok, passCount, failCount);
end
end

function [ok, msg, pdat] = localRunRelsummary(varargin)
ok = false;
msg = '';
pdat = [];
try
    [pdat, relerr] = psyrat_relsummary(varargin{:});
    if isfield(relerr, 'nogooddata') && relerr.nogooddata == 1
        msg = 'psyrat_relsummary returned nogooddata=1 (no subjects met cutoff)';
        ok = true;
        return;
    end
    ok = true;
catch ME
    msg = localErrMsg(ME);
end
end

function [ok, msg, out] = localCall(fn)
ok = false;
msg = '';
out = [];
try
    fn();
    ok = true;
catch ME
    msg = localErrMsg(ME);
end
end

function [passCount, failCount] = localCount(ok, passCount, failCount)
if ok
    passCount = passCount + 1;
else
    failCount = failCount + 1;
end
end

function row = localRow(step, scenario, analysis, combo, passed, elapsed, message)
row = {char(string(step)), char(string(scenario)), char(string(analysis)), ...
    char(string(combo)), logical(passed), elapsed, char(string(message))};
end

function t = localRowsToTable(rows)
if isempty(rows)
    t = table;
    return;
end
t = cell2table(rows, 'VariableNames', ...
    {'step','scenario','analysis','combo','passed','elapsed_sec','message'});
end

function msg = localErrMsg(ME)
parts = {};
if ~isempty(ME.identifier)
    parts{end+1} = ['[' ME.identifier ']']; %#ok<AGROW>
end
if ~isempty(ME.message)
    parts{end+1} = ME.message; %#ok<AGROW>
else
    parts{end+1} = class(ME); %#ok<AGROW>
end
if ~isempty(ME.stack)
    st = ME.stack(1);
    parts{end+1} = sprintf('(%s:%d)',st.name,st.line); %#ok<AGROW>
end
msg = strjoin(parts,' ');
if numel(msg) > 500
    msg = msg(1:500);
end
end
