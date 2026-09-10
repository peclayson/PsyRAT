function out_path = psyrat_make_ep_real_subset_fixture(varargin)
%PSYRAT_MAKE_EP_REAL_SUBSET_FIXTURE  One-time generator for the real-EP .mat regression fixture.
%
%   Derives a small slice of a REAL EP Toolkit single-trial export so that the
%   committed unit test (testImportRealEpToolkitMatExactScoring) exercises the
%   ".mat" import + scoring path on genuine recorded data. This closes the EP
%   half of GAP-2 (recorded in the GUI review notes removed in commit fa4e23d):
%   EP Toolkit files are plain MATLAB
%   structs, so the ".ept" extension is cosmetic and the same data carried as a
%   literal ".mat" must route to the EP adapter and score exactly.
%
%   The full source .ept (~260 MB) is intentionally NOT committed to the repo;
%   only the tiny derived .mat fixture (tens of KB) is. Re-run this generator if
%   the fixture ever needs to be regenerated. It is fully deterministic: fixed
%   channel/sample/trial selection, no randomness, no hidden state.
%
%   out_path = psyrat_make_ep_real_subset_fixture()
%   out_path = psyrat_make_ep_real_subset_fixture('source', PATH, 'output', PATH, ...)
%
%   Name/value options:
%     'source'      Path to a real EP Toolkit single-trial export (.ept/.mat).
%                   Required: no default is committed (the source is real
%                   recorded data). The generator errors if it is left unset.
%     'output'      Destination .mat path (default: <root>/test_data/ep_toolkit_real_subset.mat).
%     'n_channels'  Number of leading channels to keep (default 8; preserves E6).
%     'n_trials'    Number of leading single trials to keep (default 5).
%     'time_lo_ms'  Lower edge of the retained time window (default -100 ms).
%     'time_hi_ms'  Upper edge of the retained time window (default 200 ms).
%
%   The retained window must span the 0-100 ms scoring window with margin; the
%   defaults keep E1..E8 (so the canonical FCz=E6 scoring channel survives at
%   index 6), the first 5 single trials, and -100..200 ms.

opts = struct( ...
    'source', '', ...
    'output', '', ...
    'n_channels', 8, ...
    'n_trials', 5, ...
    'time_lo_ms', -100, ...
    'time_hi_ms', 200);

if mod(numel(varargin), 2) ~= 0
    error('psyrat_make_ep_real_subset_fixture:badArgs', ...
        'Inputs are incomplete. Provide name/value pairs.');
end
for i = 1:2:numel(varargin)
    key = lower(char(string(varargin{i})));
    if ~isfield(opts, key)
        error('psyrat_make_ep_real_subset_fixture:unknownOption', ...
            'Unknown option ''%s''.', key);
    end
    opts.(key) = varargin{i+1};
end

% Default output sits in the root-level test_data/ dir (same convention as
% test_data/testdata.csv). mfilename resolves tests/helpers -> tests -> root.
if isempty(opts.output)
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    opts.output = fullfile(projectRoot, 'test_data', 'ep_toolkit_real_subset.mat');
end

if exist(opts.source, 'file') ~= 2
    error('psyrat_make_ep_real_subset_fixture:sourceMissing', ...
        ['Source EP file not found: %s\n', ...
        'How to fix: this one-time generator needs a real EP Toolkit single-trial ', ...
        'export (not committed to the repo). Point ''source'' at a real .ept/.mat ', ...
        'EP file, e.g. an EP Toolkit single-trial flanker export.'], opts.source);
end

% Load the source and locate the EP struct using the same rule as the adapter
% (psyrat_adapter_epmat_read): prefer EPdata, else the first struct that carries
% a data field plus chanNames/timeNames.
S = load('-mat', opts.source);
if isfield(S, 'EPdata') && isstruct(S.EPdata)
    EP = S.EPdata;
else
    EP = [];
    fn = fieldnames(S);
    for i = 1:numel(fn)
        cand = S.(fn{i});
        if isstruct(cand) && isfield(cand, 'data') && ...
                (isfield(cand, 'chanNames') || isfield(cand, 'timeNames'))
            EP = cand;
            break;
        end
    end
    if isempty(EP)
        error('psyrat_make_ep_real_subset_fixture:noEPStruct', ...
            'Could not locate an EP Toolkit-like struct in %s', opts.source);
    end
end

if ~isfield(EP, 'data') || isempty(EP.data)
    error('psyrat_make_ep_real_subset_fixture:noData', ...
        'EP struct in %s has no data field.', opts.source);
end

% Real single-trial EP data are chan x time x trial (the subject dimension is
% singleton for these single-subject single-trial files).
data = double(EP.data);
nChan = size(data, 1);
nSamp = size(data, 2);
nWave = size(data, 3);

chanIdx = 1:min(opts.n_channels, nChan);
trialIdx = 1:min(opts.n_trials, nWave);

if ~isfield(EP, 'timeNames') || isempty(EP.timeNames)
    error('psyrat_make_ep_real_subset_fixture:noTime', ...
        'EP struct has no timeNames; cannot select a time window.');
end
timeNames = double(EP.timeNames(:).');
timeIdx = find(timeNames >= opts.time_lo_ms & timeNames <= opts.time_hi_ms);
if isempty(timeIdx)
    error('psyrat_make_ep_real_subset_fixture:emptyWindow', ...
        'No samples fell within %g..%g ms (timeNames range %g..%g).', ...
        opts.time_lo_ms, opts.time_hi_ms, min(timeNames), max(timeNames));
end

% Build a minimal EPdata struct holding only the fields the EP adapter consumes,
% but with the REAL recorded values. Keeping it minimal (no eloc/ced/history/...)
% is what holds the committed fixture to tens of KB while still being genuine
% recorded data, not synthetic.
EPdata = struct();
EPdata.dataType  = 'single_trial';
EPdata.Fs        = double(EP.Fs);
EPdata.data      = data(chanIdx, timeIdx, trialIdx);     % chan x time x trial
EPdata.timeNames = timeNames(timeIdx);
EPdata.chanNames = localNameSubset(EP, 'chanNames', chanIdx, 'Ch');
EPdata.cellNames = localNameSubset(EP, 'cellNames', trialIdx, 'cell');
EPdata.subNames  = localNameSubset(EP, 'subNames', 1:min(1, max(1, numel_safe(EP, 'subNames'))), 'sub');

outDir = fileparts(opts.output);
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
save(opts.output, 'EPdata', '-v7');
out_path = opts.output;

info = dir(opts.output);
fprintf('Wrote real-EP .mat fixture: %s\n', opts.output);
fprintf('  source           : %s\n', opts.source);
fprintf('  data dims (c/t/tr): %d x %d x %d\n', numel(chanIdx), numel(timeIdx), numel(trialIdx));
fprintf('  time window       : %g..%g ms (%d samples), Fs=%g\n', ...
    EPdata.timeNames(1), EPdata.timeNames(end), numel(timeIdx), EPdata.Fs);
fprintf('  chanNames(1:end)  : %s\n', strjoin(EPdata.chanNames(:)', ', '));
fprintf('  cellNames(1:end)  : %s\n', strjoin(EPdata.cellNames(:)', ', '));
fprintf('  file size         : %.1f KB\n', info.bytes/1024);

end

function n = numel_safe(EP, field)
if isfield(EP, field) && ~isempty(EP.(field))
    n = numel(EP.(field));
else
    n = 1;
end
end

function names = localNameSubset(EP, field, idx, defaultPrefix)
%Robustly subset an EP name list (char/string/cell/numeric) to a cell array.
if isfield(EP, field) && ~isempty(EP.(field))
    vals = EP.(field);
    if ischar(vals)
        vals = cellstr(vals);
    elseif isstring(vals)
        vals = cellstr(vals(:));
    elseif ~iscell(vals)
        vals = cellstr(string(vals(:)));
    else
        vals = cellfun(@(x) char(string(x)), vals(:), 'UniformOutput', false);
    end
    idx = idx(idx <= numel(vals));
    if ~isempty(idx)
        names = vals(idx);
        return;
    end
end
names = cell(numel(idx), 1);
for i = 1:numel(idx)
    names{i} = sprintf('%s_%d', defaultPrefix, idx(i));
end
end
