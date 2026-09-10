classdef PsyRATCmdStanIntegrationFactory
    % Builds deterministic integration tables from /test_data/testdata.csv
    % for CmdStan-backed PsyRAT model execution.

    methods (Static)
        function scenarioMap = scenarios()
            cfgs = { ...
                localCfg('ic_base', false, false, false, 1, 1, 'ic'); ...
                localCfg('ic_group', true, false, false, 1, 1, 'ic'); ...
                localCfg('ic_event', false, true, false, 1, 1, 'ic'); ...
                localCfg('ic_group_event', true, true, false, 1, 1, 'ic'); ...
                localCfg('trt_time', false, false, true, 1, 1, 'trt'); ...
                localCfg('trt_group_time', true, false, true, 1, 1, 'trt'); ...
                localCfg('trt_event_time', false, true, true, 1, 1, 'trt'); ...
                localCfg('trt_group_event_time', true, true, true, 1, 1, 'trt'); ...
                localCfg('trt_group_event_time_sserr', true, true, true, 2, 1, ...
                    'trt_sserrvar'); ...
                localCfg('ic_sserr_base', false, false, false, 2, 1, 'ic_sserrvar'); ...
                localCfg('ic_sserr_group', true, false, false, 2, 1, 'ic_sserrvar'); ...
                localCfg('ic_sserr_event', false, true, false, 2, 1, 'ic_sserrvar'); ...
                localCfg('ic_sserr_group_event', true, true, false, 2, 1, 'ic_sserrvar'); ...
                localCfg('ic_diff_event', false, true, false, 1, 2, 'ic_diff'); ...
                localCfg('ic_diff_group_event', true, true, false, 1, 2, 'ic_diff'); ...
                localCfg('ic_diff_group_event_sserr', true, true, false, 2, 2, 'ic_diff_sserrvar', '', 2) ...
                };

            scenarioMap = struct();
            for i = 1:numel(cfgs)
                scenarioMap.(cfgs{i}.name) = cfgs{i};
            end
        end

        function tbl = makeScenarioTable(projectRoot, cfg)
            raw = readtable(fullfile(projectRoot, 'test_data', 'testdata.csv'));
            raw = table( ...
                cellstr(string(raw.subjid)), ...
                cellstr(string(raw.group)), ...
                cellstr(string(raw.event)), ...
                raw.ern, ...
                'VariableNames', {'id', 'group', 'event', 'meas'});

            balanced = PsyRATCmdStanIntegrationFactory.buildBalancedSubset(raw);

            % Preserve unique subject identities when group is omitted.
            if ~cfg.includeGroup
                balanced.id = strcat(balanced.id, '__', balanced.group);
            end

            if cfg.includeTime
                balanced = PsyRATCmdStanIntegrationFactory.addOccasionColumn(balanced);
            end

            % The source data reuses subject IDs across groups.
            % When group is omitted, make IDs globally unique to avoid
            % cross-group collisions in grouping/pairing logic.
            if ~cfg.includeGroup
                balanced.id = cellstr(strcat(string(balanced.group), "__", string(balanced.id)));
            end

            cols = {'id', 'meas'};
            if cfg.includeGroup
                cols{end+1} = 'group'; %#ok<AGROW>
            end
            if cfg.includeEvent
                cols{end+1} = 'event'; %#ok<AGROW>
            end
            if cfg.includeTime
                cols{end+1} = 'time'; %#ok<AGROW>
            end

            tbl = balanced(:, cols);
            tbl.Properties.Description = cfg.name;
        end
    end

    methods (Static, Access = private)
        function balanced = buildBalancedSubset(raw)
            groups = unique(raw.group, 'stable');
            events = unique(raw.event, 'stable');
            nIdsPerGroup = 2;
            nTrialsPerEvent = 8;

            keep = false(height(raw), 1);

            for g = 1:numel(groups)
                grp = groups{g};
                ids = unique(raw.id(strcmp(raw.group, grp)), 'stable');
                picked = 0;

                for i = 1:numel(ids)
                    if picked >= nIdsPerGroup
                        break;
                    end

                    idv = ids{i};
                    ok = true;
                    for e = 1:numel(events)
                        evt = events{e};
                        n = sum(strcmp(raw.group, grp) & strcmp(raw.id, idv) & strcmp(raw.event, evt));
                        if n < nTrialsPerEvent
                            ok = false;
                            break;
                        end
                    end

                    if ~ok
                        continue;
                    end

                    for e = 1:numel(events)
                        evt = events{e};
                        idx = find(strcmp(raw.group, grp) & strcmp(raw.id, idv) & strcmp(raw.event, evt));
                        keep(idx(1:nTrialsPerEvent)) = true;
                    end
                    picked = picked + 1;
                end
            end

            balanced = raw(keep, :);
            balanced = sortrows(balanced, {'group', 'id', 'event'});
        end

        function out = addOccasionColumn(tbl)
            out = tbl;
            out.time = repmat({''}, height(out), 1);

            keys = strcat(string(out.id), "__", string(out.group), "__", string(out.event));
            uk = unique(keys, 'stable');

            for i = 1:numel(uk)
                idx = find(keys == uk(i));
                n = numel(idx);
                half = floor(n / 2);
                out.time(idx(1:half)) = {'t1'};
                out.time(idx(half+1:end)) = {'t2'};
            end
        end
    end
end

function cfg = localCfg(name, includeGroup, includeEvent, includeTime, ...
    sserrvar, diffest, expectedAnalysis, expectedErrorId, diffwpcov)
if nargin < 8
    expectedErrorId = '';
end
if nargin < 9
    diffwpcov = 1;
end

cfg = struct();
cfg.name = name;
cfg.includeGroup = includeGroup;
cfg.includeEvent = includeEvent;
cfg.includeTime = includeTime;
cfg.sserrvar = sserrvar;
cfg.diffest = diffest;
cfg.diffwpcov = diffwpcov;
cfg.diffrescor = diffwpcov;
cfg.expectedAnalysis = expectedAnalysis;
cfg.expectedErrorId = expectedErrorId;
end
