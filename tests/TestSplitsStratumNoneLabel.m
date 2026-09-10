classdef TestSplitsStratumNoneLabel < PsyRATTestBase
    % B26: psyrat_splits_summary's facet-absent test must key on the CARRIER
    % TYPE. The engine's facet-absent sentinel is always the bare char 'none'
    % (a scalar string is accepted defensively); a CELL REL.groups/REL.events
    % carries real labels by construction. The pre-fix local_isnone unwrapped
    % x{1} before comparing, and unique() sorts labels, so a group/event label
    % set whose alphabetically-first entry is literally 'none' (for example
    % {'none','patients'}) read as "facet absent" and local_stratum_rows
    % returned true(n,1): every stratum's observed design (nbar_i, n_s) was
    % silently recomputed from ALL rows while the output stayed labeled
    % per-stratum.
    %
    % The discriminator: two groups whose items-per-split weights sit around
    % 10 vs 40. Correct subsetting recovers each group's own harmonic mean;
    % the pooled arm reports one intermediate value for both strata. The
    % fixtures mirror TestSplitsNonparallel's localTrtGroupREL.

    methods (Test)

        function noneSortedFirstStillSubsetsStrata(testCase)
            % {'none','patients'}: 'none' is a REAL group label sorted first.
            REL = localTwoGroupREL({'none','patients'});
            s = psyrat_splits_summary(REL);
            testCase.verifyEqual(numel(s.strata), 2);
            for k = 1:2
                [~, nbar_i] = localGroupDesign(REL, REL.out.labels{k});
                testCase.verifyEqual(s.strata(k).nbar_i, nbar_i, ...
                    'AbsTol', 1e-12, sprintf(['Stratum ''%s'' must average ' ...
                    'only its own weights; a pooled value means the ' ...
                    'facet-absent test misread the cellstr carrier (B26).'], ...
                    REL.out.labels{k}));
            end
        end

        function noneSortedSecondKeepsSubsetting(testCase)
            % {'controls','none'}: the alphabetical trap's control arm - the
            % pre-fix code subset this correctly, and the fix must keep
            % treating a cellstr containing 'none' anywhere as real labels.
            REL = localTwoGroupREL({'controls','none'});
            s = psyrat_splits_summary(REL);
            for k = 1:2
                [~, nbar_i] = localGroupDesign(REL, REL.out.labels{k});
                testCase.verifyEqual(s.strata(k).nbar_i, nbar_i, ...
                    'AbsTol', 1e-12);
            end
        end

        function charSentinelStillReadsFacetAbsent(testCase)
            % The engine's real sentinel (bare char 'none') must still mean
            % "no group facet": one stratum, deliberately pooling all rows.
            pd = PsyRATTestDataFactory.makeSplitsTRTSummaryData();
            s = psyrat_splits_summary(pd.rel);
            testCase.verifyEqual(numel(s.strata), 1);
            w = double(pd.rel.data.weight);
            testCase.verifyEqual(s.strata(1).nbar_i, ...
                numel(w) / sum(1 ./ w), 'AbsTol', 1e-12);
        end
    end
end

function REL = localTwoGroupREL(glabels)
% Two-group trt_splits REL (mirrors TestSplitsNonparallel's localTrtGroupREL):
% each group 2 persons x 2 occasions x 3 splits; the first group's weights sit
% around 10, the second's around 40, so the per-stratum harmonic means differ
% from each other and from the pooled value.
gA = localGroupBlock(glabels{1}, [1 2], 10);
gB = localGroupBlock(glabels{2}, [3 4], 40);
REL.analysis = 'trt_splits';
REL.data = [gA; gB];
REL.groups = glabels(:);          % cellstr = REAL labels (engine contract)
REL.events = 'none';              % event facet genuinely absent (sentinel)
REL.time = {'t1','t2'};
p = struct('sig_id',[2.0 2.4],'sig_occ',[0.8 0.6],'sig_trl',[1.0 1.2],...
    'sig_trlxid',[1.5 1.1],'sig_occxid',[1.2 0.9],'sig_trlxocc',[0.5 0.4],...
    'sig_posxid',[0.9 0.7],'sig_err',[3.0 2.5]);
flds = fieldnames(p);
for k = 1:numel(flds)
    REL.out.(flds{k}) = p.(flds{k});   % 1 x 2 (one column per stratum)
end
REL.out.labels = glabels(:)';
end

function tbl = localGroupBlock(gname, ids, wbase)
% 2 persons x 2 occ x 3 splits for one group; weights vary around wbase
id = repelem(ids(:), 6);
time = repmat(repelem((1:2)', 3), 2, 1);
weight = wbase + (0:11)';
meas = (1:12)';
group = repmat({gname}, 12, 1);
tbl = table(id, meas, weight, time, group);
end

function [nsplit,nbar_i] = localGroupDesign(REL, label)
% The per-group observed design computed independently of the summary.
rows = strcmp(string(REL.data.group), string(label));
w = double(REL.data.weight(rows));
nbar_i = numel(w) / sum(1 ./ w);
key = strcat(string(REL.data.id(rows)), '_;_', string(REL.data.time(rows)));
[~,~,g] = unique(key);
nsplit = median(accumarray(g, 1));
end
