function info = psyrat_checkgithubreleases(varargin)
%Check installed PsyRAT version against GitHub releases.
%
% info = psyrat_checkgithubreleases
% info = psyrat_checkgithubreleases('currentversion','0.1.0', ...
%     'repo','peclayson/PsyRAT','verbose',true)
%
% Optional inputs:
% currentversion - semantic version for current install (default: psyrat_defineversion)
% repo - GitHub owner/repo to check (default: peclayson/PsyRAT)
% verbose - whether to print/notify the user (default: true)
% timeout - request timeout in seconds (default: 8)
% webreadfcn - function handle used to fetch API data (default: @webread)
%
% Output fields:
% statusCode - 0: update available, 1: current stable, 2: beta/dev ahead,
%              3: no stable release available, -1: unable to check
% status - text label for statusCode
% currentVersion - installed version string
% latestReleaseVersion - latest stable release semantic version (if found)
% latestReleaseTag - latest stable release tag (if found)
% latestReleaseURL - GitHub releases page for the repository
% latestReleaseAssetURL - direct download URL of the PsyRAT.zip asset attached to
%              the latest stable release ('' when there is no stable release or the
%              release carries no such asset); what psyrat_updatepsyrat downloads
% hasStableRelease - true when at least one stable release exists
% isBeta - true when running an unreleased/beta build
% isGitCheckout - true when toolbox appears to run from a git checkout
% message - what happened
% howto - actionable next step
%
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
%

if mod(length(varargin), 2)
    error('varargin:incomplete', ...
        ['Input arguments are incomplete. Provide name/value pairs. ', ...
        'How to fix: call psyrat_checkgithubreleases with paired inputs, ', ...
        'for example ''repo'',''peclayson/PsyRAT''.']);
end

opts = local_defaults();

for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};

    if ~(ischar(key) || (isstring(key) && isscalar(key)))
        error('varargin:keytype', ...
            'Input names must be character vectors or scalar strings.');
    end

    switch lower(char(key))
        case 'currentversion'
            opts.currentversion = char(string(val));
        case 'repo'
            opts.repo = char(string(val));
        case 'verbose'
            opts.verbose = logical(val);
        case 'timeout'
            opts.timeout = double(val);
        case 'webreadfcn'
            opts.webreadfcn = val;
        otherwise
            error('varargin:unknown', ...
                ['Unknown input ''%s''. How to fix: remove the input or use one ', ...
                'of: currentversion, repo, verbose, timeout, webreadfcn.'], ...
                char(string(key)));
    end
end

if isempty(opts.currentversion)
    opts.currentversion = psyrat_defineversion;
end

if ~(isnumeric(opts.timeout) && isscalar(opts.timeout) && isfinite(opts.timeout) && opts.timeout > 0)
    error('varargin:timeout', ...
        ['Input ''timeout'' must be a positive numeric scalar. ', ...
        'How to fix: pass a value such as ''timeout'', 8.']);
end

if ~isa(opts.webreadfcn, 'function_handle')
    error('varargin:webreadfcn', ...
        ['Input ''webreadfcn'' must be a function handle. ', ...
        'How to fix: pass @webread or a compatible mock function handle for tests.']);
end

info = local_empty_info(opts);

try
    releases = local_fetch_releases(opts);
catch me
    if local_is_http_404(me)
        info.statusCode = 3;
        info.status = 'no_stable_release';
        info.hasStableRelease = false;
        info.isBeta = true;

        sourceHint = ['GitHub returned 404 for the releases API. This usually means ', ...
            'the repository is private/inaccessible from this MATLAB session, ', ...
            'or no publicly accessible release metadata is available yet.'];
        info.message = sprintf(['Could not verify a stable GitHub release for ''%s''. %s ', ...
            'Treat this installation as a potentially unstable beta/development build.'], ...
            opts.repo, sourceHint);
        info.howto = sprintf(['If this is a private beta repo, this warning is expected. ', ...
            'For stable update checks, publish an accessible GitHub release at %s.'], ...
            info.latestReleaseURL);
        if opts.verbose
            warning('psyrat:betaReleaseApi404', '%s\nHow to fix: %s', ...
                info.message, info.howto);
        end
        return;
    else
        info.statusCode = -1;
        info.status = 'unable_to_check';
        info.message = sprintf(['Could not check GitHub releases for ''%s''. ', ...
            'Reason: %s.'], opts.repo, me.message);
        info.howto = ['Check your internet connection/firewall and rerun ', ...
            'psyrat_checkgithubreleases.'];
        if opts.verbose
            warning('psyrat:releaseCheckUnavailable', '%s\nHow to fix: %s', ...
                info.message, info.howto);
        end
        return;
    end
end

stable = local_select_stable_releases(releases);
if isempty(stable)
    info.statusCode = 3;
    info.status = 'no_stable_release';
    info.hasStableRelease = false;
    info.isBeta = true;

    %Wording (audit 2026-09-05): when releases exist but every one is a
    %prerelease, the user has most likely installed the official beta, so say
    %that rather than "unreleased". The "unreleased" hint is kept for a
    %repository with no releases at all.
    if ~isempty(releases)
        sourceHint = ['Only prerelease (beta) versions have been published so far, ', ...
            'so this is a beta build; there is no stable release to compare it with yet.'];
    else
        sourceHint = 'This installation should be treated as an unreleased beta build.';
    end
    if info.isGitCheckout
        sourceHint = ['This appears to be a git checkout from the repository ', ...
            '(not an official release package).'];
    end

    info.message = sprintf(['No official stable GitHub releases were found for ''%s''. ', ...
        '%s'], opts.repo, sourceHint);
    info.howto = sprintf(['Use this build with caution and monitor new releases at %s.'], ...
        info.latestReleaseURL);

    if opts.verbose
        warning('psyrat:betaNoStableRelease', '%s\nHow to fix: %s', ...
            info.message, info.howto);
    end
    return;
end

latest = stable(1);
info.hasStableRelease = true;
info.latestReleaseTag = local_char_field(latest, 'tag_name');
if isempty(info.latestReleaseTag)
    info.latestReleaseTag = local_char_field(latest, 'name');
end
info.latestReleaseVersion = local_extract_semver(info.latestReleaseTag);

releaseURL = local_char_field(latest, 'html_url');
if ~isempty(releaseURL)
    info.latestReleaseURL = releaseURL;
end

%The PsyRAT.zip asset release/PUBLISHING.md attaches to every release: the
%in-app updater downloads it through this URL instead of scraping the
%releases web page (audit finding G87, owner ruling 2026-09-09). Stable
%releases only: the prerelease filter above already dropped the betas.
info.latestReleaseAssetURL = local_find_asset_url(latest, 'PsyRAT.zip');

if isempty(info.latestReleaseVersion)
    info.statusCode = -1;
    info.status = 'unable_to_parse_release_version';
    info.message = sprintf(['Latest stable release tag ''%s'' does not contain a semantic version.'], ...
        info.latestReleaseTag);
    info.howto = ['Use release tags that include semantic versions (for example v0.2.0), ', ...
        'or pass currentversion manually for custom workflows.'];
    if opts.verbose
        warning('psyrat:releaseVersionParse', '%s\nHow to fix: %s', ...
            info.message, info.howto);
    end
    return;
end

cmp = local_compare_semver(info.currentVersion, info.latestReleaseVersion);
switch cmp
    case -1
        info.statusCode = 0;
        info.status = 'update_available';
        info.message = sprintf(['A newer PsyRAT release is available (installed: %s, latest: %s).'], ...
            info.currentVersion, info.latestReleaseVersion);
        info.howto = sprintf('Download/install the latest stable release: %s', info.latestReleaseURL);
        if opts.verbose
            fprintf('\n%s\nHow to fix: %s\n', info.message, info.howto);
        end

    case 0
        info.statusCode = 1;
        info.status = 'current_stable';
        info.message = sprintf('You are running the latest stable PsyRAT release (%s).', ...
            info.currentVersion);
        info.howto = 'No action needed.';
        if opts.verbose
            fprintf('\n%s\n', info.message);
        end

    case 1
        info.statusCode = 2;
        info.status = 'beta_ahead_of_latest_release';
        info.isBeta = true;

        sourceHint = '';
        if info.isGitCheckout
            sourceHint = [' This appears to be a repository checkout (not an official ', ...
                'release artifact).'];
        end

        info.message = sprintf(['Installed PsyRAT version (%s) is newer than the latest ', ...
            'stable release (%s). This is likely a potentially unstable beta/development build.%s'], ...
            info.currentVersion, info.latestReleaseVersion, sourceHint);
        info.howto = sprintf(['For maximum stability, use the latest official release: %s. ', ...
            'If you stay on this build, validate outputs carefully before publication.'], ...
            info.latestReleaseURL);

        if opts.verbose
            warning('psyrat:betaBuild', '%s\nHow to fix: %s', ...
                info.message, info.howto);
        end
end

end

function opts = local_defaults()
opts = struct();
opts.currentversion = '';
opts.repo = 'peclayson/PsyRAT';
opts.verbose = true;
opts.timeout = 8;
opts.webreadfcn = @webread;
end

function info = local_empty_info(opts)
info = struct();
info.statusCode = -1;
info.status = 'not_checked';
info.currentVersion = opts.currentversion;
info.latestReleaseVersion = '';
info.latestReleaseTag = '';
info.latestReleaseURL = sprintf('https://github.com/%s/releases', opts.repo);
info.latestReleaseAssetURL = '';
info.hasStableRelease = false;
info.isBeta = false;
info.isGitCheckout = local_is_git_checkout();
info.message = '';
info.howto = '';
end

function releases = local_fetch_releases(opts)
url = sprintf('https://api.github.com/repos/%s/releases?per_page=25', opts.repo);

if isequal(opts.webreadfcn, @webread)
    wopts = weboptions( ...
        'Timeout', opts.timeout, ...
        'HeaderFields', {'User-Agent', 'PsyRAT-ReleaseCheck'}, ...
        'ContentType', 'json');
    raw = opts.webreadfcn(url, wopts);
else
    raw = opts.webreadfcn(url);
end

if isempty(raw)
    releases = struct([]);
elseif isstruct(raw)
    releases = raw;
elseif iscell(raw) && all(cellfun(@isstruct, raw))
    releases = [raw{:}];
else
    error('psyrat:releaseformat', ...
        ['Unexpected GitHub releases response format. ', ...
        'How to fix: verify GitHub API access and return JSON release objects.']);
end
end

function stable = local_select_stable_releases(releases)
if isempty(releases)
    stable = releases;
    return;
end

keep = true(size(releases));
for i = 1:numel(releases)
    isDraft = local_bool_field(releases(i), 'draft');
    isPre = local_bool_field(releases(i), 'prerelease');
    keep(i) = ~(isDraft || isPre);
end

stable = releases(keep);
end

function out = local_bool_field(s, fieldname)
out = false;
if isfield(s, fieldname)
    val = s.(fieldname);
    if islogical(val) && isscalar(val)
        out = val;
    elseif isnumeric(val) && isscalar(val)
        out = (val ~= 0);
    end
end
end

function out = local_char_field(s, fieldname)
out = '';
if isfield(s, fieldname)
    raw = s.(fieldname);
    if ischar(raw)
        out = strtrim(raw);
    elseif isstring(raw) && isscalar(raw)
        out = strtrim(char(raw));
    end
end
end

function url = local_find_asset_url(release, assetname)
%Direct download URL of the asset named ASSETNAME on one release, '' when
%the release has no assets field, no assets, or none by that name. Tolerates
%the shapes webread returns for a JSON array: a struct array (same-shaped
%objects), a cell of structs (objects whose fields differ), or empty.
url = '';
if ~isstruct(release) || ~isfield(release, 'assets')
    return;
end

assets = release.assets;
if isempty(assets)
    return;
end

if iscell(assets)
    assets = assets(cellfun(@isstruct, assets));
elseif isstruct(assets)
    assets = num2cell(assets);
else
    return;
end

for i = 1:numel(assets)
    if strcmp(local_char_field(assets{i}, 'name'), assetname)
        url = local_char_field(assets{i}, 'browser_download_url');
        return;
    end
end
end

function ver = local_extract_semver(raw)
ver = '';
if isempty(raw)
    return;
end

token = regexp(char(raw), '\d+\.\d+\.\d+', 'match', 'once');
if ~isempty(token)
    ver = token;
end
end

function cmp = local_compare_semver(a, b)
aparts = sscanf(char(a), '%d.%d.%d')';
bparts = sscanf(char(b), '%d.%d.%d')';

if numel(aparts) ~= 3
    error('psyrat:versionparse', ...
        ['Current version ''%s'' is not semantic x.y.z. ', ...
        'How to fix: update psyrat_defineversion to a semantic version.'], char(a));
end

if numel(bparts) ~= 3
    error('psyrat:versionparse', ...
        ['Release version ''%s'' is not semantic x.y.z. ', ...
        'How to fix: use semantic version tags in GitHub releases.'], char(b));
end

cmp = 0;
for i = 1:3
    if aparts(i) < bparts(i)
        cmp = -1;
        return;
    elseif aparts(i) > bparts(i)
        cmp = 1;
        return;
    end
end
end

function tf = local_is_git_checkout()
tf = false;
thisFile = mfilename('fullpath');
if isempty(thisFile)
    return;
end

toolboxRoot = fileparts(fileparts(fileparts(thisFile)));
tf = isfolder(fullfile(toolboxRoot, '.git'));
end

function tf = local_is_http_404(me)
tf = false;
if isempty(me)
    return;
end

if isprop(me, 'identifier') && ~isempty(me.identifier)
    if ~isempty(regexp(me.identifier, '404', 'once'))
        tf = true;
        return;
    end
end

if isprop(me, 'message') && ~isempty(me.message)
    if ~isempty(regexp(me.message, '\b404\b', 'once'))
        tf = true;
        return;
    end
end
end
