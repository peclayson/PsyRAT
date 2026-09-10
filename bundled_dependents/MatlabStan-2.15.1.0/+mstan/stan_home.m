function d = stan_home()
% Resolve CmdStan root directory for MatlabStan.

d = '';
cands = {};

% Discovery precedence (the first candidate that is a valid CmdStan root wins):
%   1. Environment variables (PSYRAT_TEST_CMDSTAN_PATH, STAN_HOME,
%      CMDSTAN_HOME, CMDSTAN) - an explicit user/CI override.
%   2. Toolbox-managed install dirs under the user home
%      (PsyRATDependents/cmdstan-*, then ~/cmdstan-*).
%   3. Whatever is discoverable on the MATLAB path via which().
% Env vars are enumerated into cands FIRST so an explicit override always
% takes priority; the home-dir install is intentionally preferred over a path
% entry so the version psyrat_installdependents built is used in preference to
% an unrelated CmdStan that merely happens to be on the MATLAB path.
envnames = {'PSYRAT_TEST_CMDSTAN_PATH','STAN_HOME','CMDSTAN_HOME','CMDSTAN'};
for i = 1:length(envnames)
    v = strtrim(getenv(envnames{i}));
    if ~isempty(v)
        cands{end+1} = v; %#ok<AGROW>
    end
end

homeDir = char(java.lang.System.getProperty('user.home'));
if ~isempty(homeDir)
    dirPats = { ...
        fullfile(homeDir,'Documents','MATLAB','PsyRATDependents','cmdstan-*') ...
        fullfile(homeDir,'Documents','PsyRATDependents','cmdstan-*') ...
        fullfile(homeDir,'MATLAB','PsyRATDependents','cmdstan-*') ...
        fullfile(homeDir,'cmdstan-*') ...
        };
    for p = 1:length(dirPats)
        dd = dir(dirPats{p});
        dd = dd([dd.isdir]);
        % G43 ruling (2026-08-29): within each managed location, prefer the
        % NEWEST parseable cmdstan-X.Y.Z install rather than dir()'s
        % lexicographic order, under which cmdstan-2.38.0 beats
        % cmdstan-2.40.0 (digit-by-digit, '3' < '4'). Names that do not
        % parse keep their dir() order after the parseable ones.
        dd = local_sort_by_cmdstan_version(dd);
        for j = 1:length(dd)
            cands{end+1} = fullfile(dd(j).folder,dd(j).name); %#ok<AGROW>
        end
    end
end

d = local_first_cmdstan_home(cands);
if ~isempty(d)
    return;
end

v = which('runCmdStanTests.py');
if ~isempty(v)
    cands{end+1} = v; %#ok<AGROW>
end

v = which('stanc');
if ~isempty(v)
    cands{end+1} = v; %#ok<AGROW>
end

v = which('makefile');
if ~isempty(v)
    cands{end+1} = v; %#ok<AGROW>
end

pth = strsplit(path,pathsep);
for i = 1:length(pth)
    cand = fullfile(pth{i},'runCmdStanTests.py');
    if exist(cand,'file') == 2
        cands{end+1} = cand; %#ok<AGROW>
    end
end

d = local_first_cmdstan_home(cands);
end

function d = local_first_cmdstan_home(cands)
d = '';
normcands = cell(size(cands));
for i = 1:length(cands)
    normcands{i} = local_normalize_cmdstan_home(cands{i});
end
normcands = normcands(~cellfun(@isempty,normcands));
if ~isempty(normcands)
    [~,ia] = unique(normcands,'stable');
    normcands = normcands(sort(ia));
end

for i = 1:length(normcands)
    if local_is_cmdstan_home(normcands{i})
        d = normcands{i};
        return;
    end
end
end

function tf = local_is_cmdstan_home(p)
tf = exist(fullfile(p,'makefile'),'file') == 2 && ...
    exist(fullfile(p,'bin'),'dir') == 7;
end

function dd = local_sort_by_cmdstan_version(dd)
% Newest parseable cmdstan-X.Y.Z first; unparseable names after, in their
% original dir() order (G43).
ver = nan(numel(dd),3);
for k = 1:numel(dd)
    tok = regexp(dd(k).name,'^cmdstan-(\d+)\.(\d+)\.(\d+)$','tokens','once');
    if ~isempty(tok)
        ver(k,:) = cellfun(@str2double,tok);
    end
end
parseable = find(~any(isnan(ver),2));
unparseable = find(any(isnan(ver),2));
[~,ord] = sortrows(ver(parseable,:),[-1 -2 -3]);
dd = [dd(parseable(ord)); dd(unparseable)];
end

function p = local_normalize_cmdstan_home(p)
if isstring(p)
    p = char(p);
end
if isempty(p) || ~ischar(p)
    p = '';
    return;
end
p = strtrim(p);
if isempty(p)
    return;
end

if exist(p,'file') == 2
    [pp,nm,ex] = fileparts(p);
    token = lower([nm ex]);
    if strcmp(token,'runcmdstantests.py')
        p = pp;
    elseif strcmp(token,'stanc') || strcmp(token,'stanc.exe')
        p = fileparts(pp);
    elseif strcmp(token,'makefile')
        p = pp;
    else
        p = pp;
    end
elseif exist(p,'dir') == 7
    [pp,nm] = fileparts(p);
    if strcmpi(nm,'bin') && exist(fullfile(pp,'makefile'),'file') == 2
        p = pp;
    end
else
    p = '';
end
end
