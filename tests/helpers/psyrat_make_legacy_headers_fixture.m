function T = psyrat_make_legacy_headers_fixture(outfile)
%PSYRAT_MAKE_LEGACY_HEADERS_FIXTURE Generate the non-canonical-header fixture.
%
%   T = PSYRAT_MAKE_LEGACY_HEADERS_FIXTURE() regenerates
%   test_data/legacy_headers_fixture.csv and returns it as a table.
%
% PURPOSE. A small synthetic file whose HEADERS, not values, are the point:
% subjid, group, event, ern. None of those match the canonical id/meas
% contract, and 'ern' matches no Measurement candidate literal, so this file
% exercises the column auto-detect fallback in psyrat_startproc (take the
% last numeric column no other role claimed). It replaces the retired
% development fixture testdata.csv in tests/TestStartprocColumnAutodetect.m:
% that file no longer ships (owner ruling 2026-08-29), and this one does, so
% the regression stays runnable from a public release tree.
%
% DESIGN. Deliberately minimal and fully synthetic: 8 participants (4-digit
% ids 5001-5008) split between two independent groups, two events (cor, err),
% 5 trials per cell. Values are ERN-plausible microvolts (errors more
% negative than corrects) but nothing downstream depends on them.
%
% REPRODUCIBILITY. rng(20260829, 'twister'); measurements rounded to 6
% decimals before writing. Welded to the committed CSV by
% tests/TestManualDatasets.m. Standalone base MATLAB.

SEED   = 20260829;
GROUPS = {'grpA', 'grpB'};
EVENTS = {'cor', 'err'};
MU     = [-1.8, -4.6];            % per event (cor, err)
SIGP   = 1.8;                     % person SD
SIGE   = 5.0;                     % trial SD
NPERGRP = 4;
NTRL    = 5;

if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'legacy_headers_fixture.csv');
end

rng(SEED, 'twister');

subjid = []; group = {}; event = {}; ern = [];
nextid = 5001;
for g = 1:numel(GROUPS)
    for p = 1:NPERGRP
        pid = nextid; nextid = nextid + 1;
        peff = SIGP * randn;
        for e = 1:numel(EVENTS)
            scores = MU(e) + peff + SIGE * randn(NTRL, 1);
            subjid = [subjid; repmat(pid, NTRL, 1)];         %#ok<AGROW>
            group  = [group;  repmat(GROUPS(g), NTRL, 1)];   %#ok<AGROW>
            event  = [event;  repmat(EVENTS(e), NTRL, 1)];   %#ok<AGROW>
            ern    = [ern;    scores];                       %#ok<AGROW>
        end
    end
end

ern = round(ern, 6);

T = table(subjid, group, event, ern, ...
    'VariableNames', {'subjid', 'group', 'event', 'ern'});

writetable(T, outfile);
fprintf('psyrat_make_legacy_headers_fixture: wrote %d rows to %s\n', ...
    height(T), outfile);
end
