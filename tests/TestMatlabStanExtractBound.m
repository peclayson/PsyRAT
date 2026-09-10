classdef TestMatlabStanExtractBound < PsyRATTestBase
    %TESTMATLABSTANEXTRACTBOUND Regression tests for the bundled MatlabStan's
    % permuted extract (audit finding S19).
    %
    % Upstream mcmc.m indexed the draw permutation with 1:max(sz), where sz is
    % the STORED size of the parameter array [ndraws x dim1 x ...]. Because
    % permute_index has length equal to the total post-warmup draws, any
    % parameter with a single dimension larger than the draw count (the case-20
    % gamma arm's r_tid is matrix[NTID,2] with NTID = nsub*ntrl) crashed the
    % extraction with an index error AFTER the full CmdStan fit had completed
    % and been paid for. The local patch indexes 1:sz(1) - the draw dimension -
    % which is provably bit-identical whenever the upstream expression
    % succeeded: a permutation of 1..L supplies max(sz) DISTINCT row indices,
    % so the original line could only ever succeed when max(sz) == sz(1).
    %
    % No CmdStan needed: mcmc is a plain handle class, and its public
    % append_helper populates the sample store directly.

    methods (Test)

        function testWideParameterExtractsInsteadOfErroring(testCase)
            %The S19 crash shape: a parameter whose second dimension exceeds
            %the total draw count must extract, with the draw dimension
            %permuted and every draw's slice kept intact.
            [m, chains, ind] = localTwoChainFixture(testCase);
            nd = size(chains{1}.r_big,1) + size(chains{2}.r_big,1);

            out = m.extract('names','r_big');

            allr = cat(1, chains{1}.r_big, chains{2}.r_big);
            testCase.verifySize(out.r_big, size(allr));
            for k = 1:nd
                testCase.verifyEqual(squeeze(out.r_big(k,:,:)), ...
                    squeeze(allr(ind(k),:,:)), 'AbsTol', 0, sprintf( ...
                    'Row %d must be draw %d of the concatenated store.', ...
                    k, ind(k)));
            end
        end

        function testSmallParameterPermutationIsTheDocumentedOne(testCase)
            %Behavior preservation: for a parameter where the upstream
            %expression worked (every dimension <= ndraws), the returned
            %permutation must be exactly randperm(ndraws) evaluated at the
            %object's construction seed - the same vector the upstream code
            %produced, so previously working extractions are bit-identical.
            [m, chains, ind] = localTwoChainFixture(testCase);

            out = m.extract('names','b');

            allb = cat(1, chains{1}.b, chains{2}.b);
            testCase.verifyEqual(out.b, allb(ind), 'AbsTol', 0);
        end

        function testNamedDiagnosticExtractTouchesOnlyThatField(testCase)
            %The production divergence counter (local_count_divergences in
            %psyrat_computevarcomp) must request 'pars','divergent__' rather
            %than extracting everything: a no-names extract permutes and
            %copies EVERY stored parameter - ~2x the whole sample store,
            %post-fit - to read one [ndraws x 1] vector (2026-08-11 code
            %review). This pins the mechanism that fix relies on: a
            %diagnostics-named field is extractable by name, alone, and its
            %sum is permutation-invariant. (StanFit.extract forwards its
            %'pars' argument to this mcmc-level 'names' parameter.)
            [m, chains, ~] = localTwoChainFixture(testCase);

            out = m.extract('names','divergent__');

            testCase.verifyEqual(fieldnames(out), {'divergent__'}, ...
                'A named extract must return only the requested field.');
            testCase.verifyEqual(sum(out.divergent__(:)), ...
                sum(chains{1}.divergent__) + sum(chains{2}.divergent__), ...
                'The divergence count is permutation-invariant.');
        end

        function testNoNamesExtractCoversEveryParameter(testCase)
            %The divergence-count diagnostics path calls extract with NO names,
            %which iterates every stored parameter. Before the patch that call
            %crashed on the first wide parameter (swallowed by a try/catch in
            %production, silently reporting the divergence count as NaN for
            %every Gaussian two-facet fit at realistic scale). It must now
            %return every field.
            [m, chains, ~] = localTwoChainFixture(testCase);

            out = m.extract();

            testCase.verifyEqual(sort(fieldnames(out)), ...
                sort(fieldnames(chains{1})));
            testCase.verifySize(out.r_big, ...
                size(cat(1, chains{1}.r_big, chains{2}.r_big)));
            testCase.verifySize(out.b, [6 1]);
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function [m, chains, ind] = localTwoChainFixture(testCase)
%LOCALTWOCHAINFIXTURE A two-chain mcmc object whose store holds one narrow
% parameter (b, vector) and one wide one (r_big, stored [draws x 10 x 2] with
% 10 > 6 total draws - the S19 shape). Values are sequential integers so a
% permutation defect cannot pass by coincidence. Returns the expected
% permutation, reproduced from the object's construction seed exactly the way
% mcmc.permute_index does (rng restored afterward).
root = testCase.projectRoot();
addpath(genpath(fullfile(root, 'bundled_dependents', 'MatlabStan-2.15.1.0')));

seed = 4242;
m = mcmc(seed);

ndPerChain = 3;
N = 10;
chains = cell(1,2);
for c = 1:2
    base = (c-1) * 1000;
    s = struct();
    s.b = base + (1:ndPerChain)';
    s.r_big = base + reshape(1:(ndPerChain*N*2), ndPerChain, N, 2);
    %sampler diagnostic column, stored in samples exactly like a parameter
    %(read_stan_csv keeps the whole CSV header) - what the production
    %divergence counter extracts by name.
    s.divergent__ = double((1:ndPerChain)' == c);
    chains{c} = s;
end
%append_helper with an empty store installs the whole 1x2 struct array.
store = [chains{1}, chains{2}];
m.append_helper('samples', store, 1);

%the expected permutation: mcmc(seed) stores rng state right after rng(seed),
%and permute_index restores that state before randperm(total draws).
curr = rng;
rng(seed);
ind = randperm(2 * ndPerChain);
rng(curr);
end
