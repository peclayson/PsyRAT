classdef TestSplitsNonparallel < PsyRATTestBase
    % Deterministic unit tests for the nonparallel data-splits observed-design
    % summary (psyrat_splits_summary; Rocha Tables 4/5, analysis cases 23/24). No
    % CmdStan/MatlabStan required: the summary is checked against the independent
    % closed forms in PsyRATAccuracyOracle (splitsSingle/splitsTrt + the SEM
    % oracles). The live known-truth CmdStan recovery is a separate slow test.
    %
    % The summary recomputes the observed design (n_s splits per person, n_o
    % occasions, harmonic-mean items-per-split n_i) from REL.data; the tests build
    % a synthetic REL with a KNOWN design and verify both that the summary
    % recovered that design and that the coefficients/SEMs match the oracle fed
    % the same design. Single-draw inputs collapse the credible interval to the
    % point estimate for exact checks; one test uses a draw vector to exercise the
    % interval path.

    methods (Test)

        % ---- single occasion (ic_splits, Rocha Table 4) ----------------------

        function testSingleOccasionMatchesOracle(testCase)
            p = localSinglePars();
            REL = localSingleREL(p);
            [nsplit,nbar_i] = localSingleDesign(REL);
            s = psyrat_splits_summary(REL, 'CI', .95);

            % the summary recovered the known observed design
            testCase.verifyEqual(s.strata(1).nsplit, nsplit);
            testCase.verifyEqual(s.strata(1).nocc, 1);
            testCase.verifyEqual(s.strata(1).nbar_i, nbar_i, 'AbsTol', 1e-12);

            c = s.strata(1).coef(1);
            D = row3(@() PsyRATAccuracyOracle.splitsSingle(1, p.sig_u, ...
                p.sig_trl, p.sig_ps, p.sig_err, nsplit, nbar_i, .95));
            G = row3(@() PsyRATAccuracyOracle.splitsSingle(2, p.sig_u, ...
                p.sig_trl, p.sig_ps, p.sig_err, nsplit, nbar_i, .95));
            testCase.verifyEqual(c.D, D, 'AbsTol', 1e-10);
            testCase.verifyEqual(c.G, G, 'AbsTol', 1e-10);

            semabs = row3(@() PsyRATAccuracyOracle.splitsSingleSem('abs', ...
                p.sig_trl, p.sig_ps, p.sig_err, nsplit, nbar_i, .95));
            semrel = row3(@() PsyRATAccuracyOracle.splitsSingleSem('rel', ...
                p.sig_trl, p.sig_ps, p.sig_err, nsplit, nbar_i, .95));
            testCase.verifyEqual(c.SEM_abs, semabs, 'AbsTol', 1e-10);
            testCase.verifyEqual(c.SEM_rel, semrel, 'AbsTol', 1e-10);
            testCase.verifyTrue(c.gbias);
        end

        function testSingleOccasionDependabilityNotAboveGeneralizability(testCase)
            % Reported D includes the split main effect in the error and G does
            % not, so D <= G for the reported (observed-design) values.
            REL = localSingleREL(localSinglePars());
            s = psyrat_splits_summary(REL);
            c = s.strata(1).coef(1);
            testCase.verifyLessThanOrEqual(c.D(2), c.G(2) + 1e-12);
        end

        function testSingleOccasionIntervalPath(testCase)
            % A draw vector: the summary's [ll pt ul] must match the oracle's
            % over the same draws (exercises mean + quantile, not just a scalar).
            p = localSinglePars();
            p.sig_u  = 1.8 + 0.4 * randn(300,1);
            p.sig_trl = abs(0.9 + 0.3 * randn(300,1));
            p.sig_ps = abs(1.4 + 0.3 * randn(300,1));
            p.sig_err = abs(2.8 + 0.5 * randn(300,1));
            REL = localSingleREL(p);
            [nsplit,nbar_i] = localSingleDesign(REL);
            s = psyrat_splits_summary(REL, 'CI', .90);
            c = s.strata(1).coef(1);
            G = row3(@() PsyRATAccuracyOracle.splitsSingle(2, p.sig_u, ...
                p.sig_trl, p.sig_ps, p.sig_err, nsplit, nbar_i, .90));
            testCase.verifyEqual(c.G, G, 'AbsTol', 1e-12);
        end

        % ---- multi occasion (trt_splits, Rocha Table 5) ----------------------

        function testMultiOccasionMatchesOracle(testCase)
            p = localTrtPars();
            REL = localTrtREL(p);
            [nsplit,nbar_i,nocc] = localTrtDesign(REL);
            s = psyrat_splits_summary(REL, 'CI', .95);

            testCase.verifyEqual(s.strata(1).nsplit, nsplit);
            testCase.verifyEqual(s.strata(1).nocc, nocc);
            testCase.verifyEqual(s.strata(1).nbar_i, nbar_i, 'AbsTol', 1e-12);
            testCase.verifyNumElements(s.strata(1).coef, 3); % CE/CS/CES

            for r = 1:3
                c = s.strata(1).coef(r);
                D = row3(@() PsyRATAccuracyOracle.splitsTrt(1, r, p.sig_id, ...
                    p.sig_occ, p.sig_trl, p.sig_trlxid, p.sig_occxid, ...
                    p.sig_trlxocc, p.sig_posxid, p.sig_err, nsplit, nocc, ...
                    nbar_i, .95));
                G = row3(@() PsyRATAccuracyOracle.splitsTrt(2, r, p.sig_id, ...
                    p.sig_occ, p.sig_trl, p.sig_trlxid, p.sig_occxid, ...
                    p.sig_trlxocc, p.sig_posxid, p.sig_err, nsplit, nocc, ...
                    nbar_i, .95));
                testCase.verifyEqual(c.D, D, 'AbsTol', 1e-10);
                testCase.verifyEqual(c.G, G, 'AbsTol', 1e-10);
                testCase.verifyLessThanOrEqual(c.D(2), c.G(2) + 1e-12);

                semabs = row3(@() PsyRATAccuracyOracle.splitsTrtSem('abs', r, ...
                    p.sig_occ, p.sig_trl, p.sig_trlxid, p.sig_occxid, ...
                    p.sig_trlxocc, p.sig_posxid, p.sig_err, nsplit, nocc, ...
                    nbar_i, .95));
                testCase.verifyEqual(c.SEM_abs, semabs, 'AbsTol', 1e-10);
                testCase.verifyTrue(c.gbias);
            end
        end

        function testMultiOccasionStrataSubsetting(testCase)
            % Two groups with different per-split item counts: the summary must
            % subset REL.data per stratum and use each stratum's own n_i / n_s.
            p = localTrtGroupPars();
            REL = localTrtGroupREL(p);
            s = psyrat_splits_summary(REL, 'CI', .95);
            testCase.verifyNumElements(s.strata, 2);

            for gi = 1:2
                [nsplit,nbar_i,nocc] = localTrtGroupDesign(REL, ...
                    REL.out.labels{gi});
                testCase.verifyEqual(s.strata(gi).nbar_i, nbar_i, 'AbsTol', 1e-12);
                testCase.verifyEqual(s.strata(gi).nsplit, nsplit);
                G = row3(@() PsyRATAccuracyOracle.splitsTrt(2, 3, p.sig_id(:,gi), ...
                    p.sig_occ(:,gi), p.sig_trl(:,gi), p.sig_trlxid(:,gi), ...
                    p.sig_occxid(:,gi), p.sig_trlxocc(:,gi), p.sig_posxid(:,gi), ...
                    p.sig_err(:,gi), nsplit, nocc, nbar_i, .95));
                testCase.verifyEqual(s.strata(gi).coef(3).G, G, 'AbsTol', 1e-10);
            end
            % the two strata have different harmonic-mean item counts
            testCase.verifyNotEqual(s.strata(1).nbar_i, s.strata(2).nbar_i);
        end

        % ---- options + structure --------------------------------------------

        function testObsAndNoccOverrides(testCase)
            REL = localTrtREL(localTrtPars());
            s = psyrat_splits_summary(REL, 'obs', 10, 'nocc', 5);
            testCase.verifyEqual(s.strata(1).nsplit, 10);
            testCase.verifyEqual(s.strata(1).nocc, 5);
        end

        function testNoccObservedKeyword(testCase)
            REL = localTrtREL(localTrtPars());
            [~,~,nocc] = localTrtDesign(REL);
            s = psyrat_splits_summary(REL, 'nocc', 'observed');
            testCase.verifyEqual(s.strata(1).nocc, nocc);
        end

        function testVarcompTableShape(testCase)
            REL = localSingleREL(localSinglePars());
            s = psyrat_splits_summary(REL);
            T = s.strata(1).varcomp;
            testCase.verifyEqual(T.Properties.VariableNames, ...
                {'component','symbol','estimate','ci_lower','ci_upper'});
            % single-occasion design reports 4 components (p, s, ps, residual)
            testCase.verifyEqual(height(T), 4);
            testCase.verifyTrue(any(strcmp(T.symbol, 'sigma^2_p')));

            RELt = localTrtREL(localTrtPars());
            St = psyrat_splits_summary(RELt);
            testCase.verifyEqual(height(St.strata(1).varcomp), 8);
        end

        function testRejectsBadAnalysis(testCase)
            REL.analysis = 'ic_base';
            testCase.verifyError(@() psyrat_splits_summary(REL), ...
                'splits:badinput');
        end

        function testRejectsNonScalarObs(testCase)
            % A range obs would silently malform the [ll pt ul] coefficients;
            % the observed-design summary takes a scalar only.
            REL = localSingleREL(localSinglePars());
            testCase.verifyError(@() psyrat_splits_summary(REL, 'obs', [10 50]), ...
                'splits:badobs');
            RELt = localTrtREL(localTrtPars());
            testCase.verifyError(@() psyrat_splits_summary(RELt, 'nocc', [2 4]), ...
                'splits:badnocc');
        end

    end
end

% =============================== local helpers ===============================

function p = localSinglePars()
% scalar component "draws" (single draw -> CI collapses to the point estimate)
p.sig_u   = 2.0;   % sigma_p
p.sig_trl = 1.0;   % sigma_s (split main)
p.sig_ps  = 1.5;   % sigma_ps (person x split)
p.sig_err = 3.0;   % per-item residual
end

function REL = localSingleREL(p)
% ic_splits, one stratum, 3 persons x 4 splits, varying items-per-split (n_i)
id = repelem((1:3)', 4);
meas = (1:12)';                                  % arbitrary; unused by the calc
weight = [10;20;30;40; 15;25;35;45; 12;22;32;42]; % n_i variation
REL.analysis = 'ic_splits';
REL.data = table(id, meas, weight);
REL.groups = 'none';
REL.events = 'none';
REL.out.sig_u        = p.sig_u(:);
REL.out.sig_trl      = p.sig_trl(:);
REL.out.sig_splitxid = p.sig_ps(:);
REL.out.sig_err      = p.sig_err(:);
REL.out.labels = 'measure';
end

function [nsplit,nbar_i] = localSingleDesign(REL)
w = double(REL.data.weight);
nbar_i = numel(w) / sum(1 ./ w);
nsplit = median(localCounts(REL.data.id));
end

function p = localTrtPars()
p.sig_id       = 2.0;
p.sig_occ      = 0.8;
p.sig_trl      = 1.0;
p.sig_trlxid   = 1.5;
p.sig_occxid   = 1.2;
p.sig_trlxocc  = 0.5;
p.sig_posxid   = 0.9;
p.sig_err      = 3.0;
end

function REL = localTrtREL(p)
% trt_splits, one stratum (no group/event), 3 persons x 2 occ x 4 splits
id = repelem((1:3)', 8);
time = repmat(repelem((1:2)', 4), 3, 1);
meas = (1:24)';
weight = (10:33)';                               % n_i variation across rows
REL.analysis = 'trt_splits';
REL.data = table(id, meas, weight, time);
REL.groups = 'none';
REL.events = 'none';
REL.time = {'t1','t2'};
REL.out.sig_id       = p.sig_id(:);
REL.out.sig_occ      = p.sig_occ(:);
REL.out.sig_trl      = p.sig_trl(:);
REL.out.sig_trlxid   = p.sig_trlxid(:);
REL.out.sig_occxid   = p.sig_occxid(:);
REL.out.sig_trlxocc  = p.sig_trlxocc(:);
REL.out.sig_posxid   = p.sig_posxid(:);
REL.out.sig_err      = p.sig_err(:);
REL.out.labels = {'none'};
end

function [nsplit,nbar_i,nocc] = localTrtDesign(REL)
w = double(REL.data.weight);
nbar_i = numel(w) / sum(1 ./ w);
key = strcat(string(REL.data.id), '_;_', string(REL.data.time));
nsplit = median(localCounts(key));
nocc = numel(unique(REL.data.time));
end

function p = localTrtGroupPars()
% two strata (columns): different magnitudes so coefficients differ by group
p.sig_id      = [2.0 2.4];
p.sig_occ     = [0.8 0.6];
p.sig_trl     = [1.0 1.2];
p.sig_trlxid  = [1.5 1.1];
p.sig_occxid  = [1.2 0.9];
p.sig_trlxocc = [0.5 0.4];
p.sig_posxid  = [0.9 0.7];
p.sig_err     = [3.0 2.5];
end

function REL = localTrtGroupREL(p)
% two groups, each 2 persons x 2 occ x 3 splits; group B has larger n_i
gA = localGroupBlock('A', [1 2], 10);
gB = localGroupBlock('B', [3 4], 40);
REL.analysis = 'trt_splits';
REL.data = [gA; gB];
REL.groups = {'A','B'};
REL.events = 'none';
REL.time = {'t1','t2'};
flds = {'sig_id','sig_occ','sig_trl','sig_trlxid','sig_occxid',...
    'sig_trlxocc','sig_posxid','sig_err'};
for k = 1:numel(flds)
    REL.out.(flds{k}) = p.(flds{k});   % 1 x 2 (one column per stratum)
end
REL.out.labels = {'A','B'};
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

function [nsplit,nbar_i,nocc] = localTrtGroupDesign(REL, label)
rows = strcmp(string(REL.data.group), string(label));
w = double(REL.data.weight(rows));
nbar_i = numel(w) / sum(1 ./ w);
key = strcat(string(REL.data.id(rows)), '_;_', string(REL.data.time(rows)));
nsplit = median(localCounts(key));
nocc = numel(unique(REL.data.time(rows)));
end

function counts = localCounts(key)
[~,~,g] = unique(key);
counts = accumarray(g, 1);
end

function v = row3(fn)
% Capture the oracle's three outputs [ll,pt,ul] into a [ll pt ul] row vector to
% compare against the summary's coefficient/SEM fields (stored as [ll pt ul]).
[a,b,c] = fn();
v = [a b c];
end
