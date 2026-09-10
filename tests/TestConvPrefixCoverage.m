classdef TestConvPrefixCoverage < PsyRATTestBase
    % Pins the convergence-monitoring parameter-name lists in psyrat_convprefix
    % against the parameters the corresponding Stan models actually declare.
    %
    % WHY THIS EXISTS
    % psyrat_storeconv keeps only those rows of the CmdStan summary whose
    % parameter label passes psyrat_hasprefix(label, psyrat_convprefix(trt)).
    % psyrat_hasprefix is a PREFIX match (strncmp), so a list entry that names no
    % real parameter silently monitors nothing, and a real parameter that no entry
    % prefixes is silently never checked for R-hat / n_eff. Neither failure mode
    % is visible from the outside: estimation still runs, the convergence table is
    % still produced, and the design still reports "converged" - just on fewer
    % parameters than the user is told.
    %
    % That is exactly how the analysis-1 list drifted: it read
    % {'lp__' 'mu_' 'sig_u_' 'sig_e_'} while ic_base.stan declares mu, sig_u,
    % sig_trl, sig_e. strncmp('mu','mu_',3) is false, so the toolbox's most-used
    % design reported its verdict on the joint log-density alone.
    %
    % HOW IT WORKS
    % psyrat_convprefix is a LOCAL function of psyrat_computevarcomp.m and cannot
    % be called from outside that file, so this test reads the literal lists out
    % of the source text and compares them to the parameter names declared in the
    % committed Stan goldens (tests/golden/stan/*.stan). Both artifacts are
    % already pinned elsewhere - the goldens by TestStanModelSnapshots - so this
    % test adds the missing link BETWEEN them. No CmdStan, no MATLAB toolboxes.
    %
    % Only the non-threaded goldens are read: the within-chain threaded twin
    % differs only in the model block (reduce_sum), never in the declarations.

    methods (Test)
        function testEveryPrefixNamesARealParameter(testCase)
            % Direction 1: no dead entries. For every (trt branch, golden) pair,
            % each prefix must match at least one name declared in that model's
            % parameters / transformed parameters / generated quantities blocks
            % (stansummary reports all three), unless the pair is listed in
            % localKnownUnmatched below with a stated reason.
            prefixLists = localParseConvPrefix(testCase);
            goldenMap = localGoldenMap();
            exempt = localKnownUnmatched();

            branches = cell2mat(keys(prefixLists));
            for b = branches
                goldens = goldenMap(b);
                for g = 1:numel(goldens)
                    name = goldens{g};
                    declared = localDeclaredNames(testCase, name);
                    prefixes = prefixLists(b);

                    dead = {};
                    for p = 1:numel(prefixes)
                        if strcmp(prefixes{p}, 'lp__')
                            continue;   % sampler quantity, never a model declaration
                        end
                        if ~any(strncmp(declared, prefixes{p}, numel(prefixes{p})))
                            dead{end+1} = prefixes{p}; %#ok<AGROW>
                        end
                    end

                    expected = localLookupExemption(exempt, b, name);
                    testCase.verifyEqual(sort(dead(:)), sort(expected(:)), ...
                        sprintf(['psyrat_convprefix(%d) vs golden %s.stan: the set of ', ...
                        'prefixes matching NO declared parameter changed. A prefix that ', ...
                        'matches nothing monitors nothing (psyrat_hasprefix is strncmp). ', ...
                        'If you fixed the list, delete the matching entry from ', ...
                        'localKnownUnmatched; if you added a design, add its reason there.'], ...
                        b, name));
                end
            end
        end

        function testEveryScalarVarianceComponentIsMonitored(testCase)
            % Direction 2: no unmonitored scalar. Every free parameter declared as
            % a bare scalar `real` (the grand means, the variance components, the
            % scalar correlations - the quantities a user would expect a
            % convergence verdict on) must be matched by some prefix.
            %
            % Restricted to scalar `real` on purpose: the vector/matrix free
            % parameters in these models are the non-centered auxiliaries (z_*,
            % *_raw, L_*, gro_effs_stndzd), which every branch deliberately leaves
            % out of the table because they are NSUB/NTRL-dimensional. Asserting
            % over them would be asserting a convention, not a defect.
            prefixLists = localParseConvPrefix(testCase);
            goldenMap = localGoldenMap();
            known = localKnownCoverageGaps();

            branches = cell2mat(keys(prefixLists));
            for b = branches
                goldens = goldenMap(b);
                for g = 1:numel(goldens)
                    name = goldens{g};
                    scalars = localScalarRealParameters(testCase, name);
                    prefixes = prefixLists(b);

                    gaps = {};
                    for s = 1:numel(scalars)
                        hit = false;
                        for p = 1:numel(prefixes)
                            if strncmp(scalars{s}, prefixes{p}, numel(prefixes{p}))
                                hit = true;
                                break;
                            end
                        end
                        if ~hit
                            gaps{end+1} = scalars{s}; %#ok<AGROW>
                        end
                    end

                    expected = localLookupExemption(known, b, name);
                    testCase.verifyEqual(sort(gaps(:)), sort(expected(:)), ...
                        sprintf(['psyrat_convprefix(%d) vs golden %s.stan: the set of ', ...
                        'scalar `real` parameters matched by NO prefix changed. An ', ...
                        'unmatched variance component never gets an R-hat / n_eff check, ', ...
                        'so the reported convergence verdict silently excludes it. If you ', ...
                        'closed a gap, delete the matching entry from ', ...
                        'localKnownCoverageGaps.'], b, name));
                end
            end
        end

        function testEveryBranchIsMappedToAGolden(testCase)
            % Guard: the two tests above are only as complete as localGoldenMap.
            % Assert the map covers exactly the branches psyrat_convprefix defines,
            % so adding a design without deciding which golden pins it fails here
            % instead of quietly going unchecked.
            prefixLists = localParseConvPrefix(testCase);
            goldenMap = localGoldenMap();

            testCase.verifyEqual(sort(cell2mat(keys(prefixLists))), ...
                sort(cell2mat(keys(goldenMap))), ...
                ['psyrat_convprefix branches and localGoldenMap disagree. Every trt ', ...
                'code with a prefix list needs a golden mapping (or an explicit ', ...
                'entry documenting why it cannot be mapped).']);
        end

        function testEveryCoefficientBlockIsMonitored(testCase)
            % Direction 3: no unmonitored coefficient block. The population-level
            % effect vectors (b, b_sigma, b_nu, and the event-crossed *_dim
            % variants) are the mean and dispersion SLOPES - for the dynamic /
            % conditional-reliability designs they are the estimand the analysis
            % exists to report.
            %
            % WHY THIS IS A THIRD TEST AND NOT PART OF DIRECTION 2. Direction 2 is
            % restricted to scalar `real` parameters, deliberately, because the
            % other vector/matrix free parameters are non-centered auxiliaries
            % (z_*, *_raw, L_*) that no branch monitors. The coefficient blocks are
            % vectors too, so they fell through BOTH existing directions: S15 found
            % code 3 omitting b / b_sigma while its two-facet sibling code 7 listed
            % them, and neither test could see it. Direction 1 pins them only while
            % they are present - deleting them again would fail nothing. This test
            % is what makes that re-drift fail.
            %
            % No exemption table: every mapped golden that declares one of these
            % blocks has a branch that monitors it. If that ever stops being true,
            % the right response is almost certainly to fix the list.
            prefixLists = localParseConvPrefix(testCase);
            goldenMap = localGoldenMap();

            branches = cell2mat(keys(prefixLists));
            for b = branches
                goldens = goldenMap(b);
                for g = 1:numel(goldens)
                    name = goldens{g};
                    blocks = localCoefficientBlockParameters(testCase, name);
                    prefixes = prefixLists(b);

                    for c = 1:numel(blocks)
                        hit = any(cellfun(@(p) strncmp(blocks{c}, p, numel(p)), ...
                            prefixes));
                        testCase.verifyTrue(hit, sprintf( ...
                            ['psyrat_convprefix(%d) vs golden %s.stan: the ', ...
                            'coefficient block ''%s'' is declared but matched by ', ...
                            'no prefix, so it never gets an R-hat / n_eff check. ', ...
                            'For a dynamic-reliability design these slopes are the ', ...
                            'reported estimand. Add it to the list for code %d.'], ...
                            b, name, blocks{c}, b));
                    end
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% trt code -> the committed Stan goldens that code's convergence list serves.
%
% Derived from the psyrat_storeconv(fit, <code>) call sites in
% psyrat_computevarcomp.m and the scenario definitions in
% tests/helpers/psyrat_snapshot_scenarios.m + PsyRATCmdStanIntegrationFactory.m.
% Several codes are reused by more than one analysis case, which is why some
% entries mix families - that reuse is the main source of the drift these tests
% detect, so it is represented faithfully rather than tidied away.
% -------------------------------------------------------------------------
function m = localGoldenMap()
pairs = { ...
    ...% Gaussian one-facet family, analyses 1-4 (psyrat_storeconv(fit) default)
    0,  {'ic_base', 'ic_base_priors', 'ic_group', 'ic_event', ...
         'ic_event_priors', 'ic_group_event'}; ...
    ...% Gaussian crossed test-retest, analysis 5
    1,  {'trt_time', 'trt_time_priors', 'trt_group_time', 'trt_event_time', ...
         'trt_group_event_time'}; ...
    ...% Gaussian subject-level error variance (analysis 6) AND the Gaussian
    ...% one-facet dynamic/conditional model (analysis 11), which reuses code 3
    3,  {'ic_sserr_base', 'ic_sserr_group', 'ic_sserr_event', ...
         'ic_sserr_group_event', 'dynrel_onefacet'}; ...
    ...% Gaussian two-event difference (analysis 7) AND the four-event
    ...% difference-of-differences (analysis 9), which reuses code 4
    4,  {'ic_diff_event', 'ic_diff_event_priors', 'ic_diff_group_event', ...
         'ic_diff_rescor_event', 'dod_event', 'dod_event_priors'}; ...
    ...% Gaussian subject-level two-event difference, analysis 8. Both goldens
    ...% capture the CONCURRENT (diffwpcov=2) form; the non-concurrent case-8
    ...% model has no golden, so its prefixes are unchecked here.
    5,  {'ic_diff_group_event_sserr', 'ic_diff_group_event_sserr_priors'}; ...
    ...% Gaussian two-facet difference (analysis 10) AND the one-facet dynamic
    ...% difference designs (analyses 12, 13), which reuse code 6
    6,  {'trt_diff_event_time', 'trt_diff_event_time_priors', ...
         'trt_diff_rescor_event_time', 'dynrel_diff_onefacet', ...
         'dynrel_diff_sserr_onefacet'}; ...
    7,  {'dynrel_trt_onefacet'}; ...                       % analysis 14
    8,  {'dynrel_diff_trt_onefacet'}; ...                  % analysis 15
    9,  {'dynrel_diff_sserr_trt_onefacet'}; ...            % analysis 16
    10, {'dynrel_diff_rescor_onefacet'}; ...               % analysis 17
    11, {'dynrel_diff_rescor_trt_onefacet'}; ...           % analysis 18
    12, {'dynrel_diff_sserr_rescor_onefacet'}; ...         % analysis 19
    13, {'dynrel_diff_sserr_rescor_trt_onefacet'}; ...     % analysis 20
    14, {'dynrel_diff_sserr_ssrescor_onefacet'}; ...       % analysis 21
    15, {'dynrel_diff_sserr_ssrescor_trt_onefacet'}; ...   % analysis 22
    23, {'splits_iins_single'}; ...                        % analysis 23
    24, {'splits_iins_trt'}; ...                           % analysis 24
    25, {'trt_group_event_time_sserr'}; ...                % analysis 25 (Gaussian)
    ...% Gamma / scaled-chi-square family. Analysis 4 (group x event) has no
    ...% gamma golden, so code 30 is pinned by the three that exist.
    30, {'ic_gamma_base', 'ic_gamma_group', 'ic_gamma_event'}; ...
    31, {'ic_gamma_diff_event'}; ...                       % analysis 7  (gamma)
    32, {'ic_gamma_trt'}; ...                              % analysis 5  (gamma)
    33, {'ic_gamma_diff_trt'}; ...                         % analysis 10 (gamma)
    34, {'ic_gamma_dodiff'}; ...                           % analysis 9  (gamma)
    35, {'dynrel_gamma_onefacet'}; ...                     % analysis 11 (gamma)
    36, {'dynrel_gamma_trt_onefacet'}; ...                 % analysis 14 (gamma)
    37, {'ic_gamma_diff_sserr'}; ...                       % analysis 8  (gamma)
    38, {'ic_gamma_diff_copula', 'ic_gamma_diff_copula_fixnu'}; ...        % 7 concurrent
    39, {'ic_gamma_diff_trt_copula', 'ic_gamma_diff_trt_copula_fixnu'}; ...% 10 concurrent
    40, {'ic_gamma_diff_sserr_copula'}; ...                % analysis 8  concurrent
    41, {'dynrel_gamma_diff_onefacet'}; ...                % analysis 12 (gamma)
    42, {'ic_gamma_trt_sserr'}; ...                        % analysis 25 (gamma)
    43, {'dynrel_gamma_sserr_onefacet'}; ...               % analysis 26 (gamma)
    44, {'dynrel_gamma_sserr_trt_onefacet'}; ...           % analysis 27 (gamma)
    ...% Gamma LOCATION-SCALE dynrel (gammascale=2): the log-nu scale submodel is
    ...% replaced by a direct log residual-SD submodel, so these need their own
    ...% codes - 'Intercept_nu'/'b_nu' prefix nothing in these goldens.
    45, {'dynrel_gamma_ls_onefacet'}; ...                  % analysis 11 (gamma, LS)
    46, {'dynrel_gamma_ls_sserr_onefacet'}; ...            % analysis 26 (gamma, LS)
    ...% The STATIC subject-level location-scale designs (owner ruling H,
    ...% 2026-08-04). Same argument as 45/46: the log-nu codes 30 and 42 name
    ...% 'Intercept_nu', which prefixes nothing in these goldens, so reusing them
    ...% would monitor the scale submodel not at all while still reporting a
    ...% verdict. These lists carry no b/b_sigma because the static designs have no
    ...% dimension slopes.
    47, {'ic_gamma_ls_sserr'}; ...                         % analysis 6  (gamma, LS)
    48, {'ic_gamma_ls_trt_sserr'}; ...                     % analysis 25 (gamma, LS)
    ...% The STATIC two-event DIFFERENCE location-scale designs (same ruling).
    ...% Code 31 names 'b_nu', which prefixes nothing in these goldens, so it
    ...% cannot be reused. Unlike 47/48 these DO carry b/b_sigma - the per-event
    ...% intercept vectors - but still no dimension slopes. 50's golden is
    ...% byte-identical to 49's by delegation, exactly as 37's is to 31's.
    49, {'ic_gamma_ls_diff_event'}; ...                    % analysis 7  (gamma, LS)
    50, {'ic_gamma_ls_diff_sserr'}; ...                    % analysis 8  (gamma, LS)
    ...% The STATIC TWO-FACET difference location-scale design (same ruling).
    ...% Code 33 names 'b_nu', which prefixes nothing in this golden. Unlike 49/50
    ...% this list also carries the five crossed-factor SDs and the six cross-event
    ...% correlations, because the two-facet mean submodel is unchanged from 33.
    ...% There is no delegation twin: analysis 10 has no subject-level variant.
    51, {'ic_gamma_ls_diff_trt'}; ...                      % analysis 10 (gamma, LS)
    ...% The four-event DIFFERENCE-OF-DIFFERENCES location-scale design (same
    ...% ruling), which completes the static family. Code 34 names 'b_nu', which
    ...% prefixes nothing in this golden. Like 34 it excludes the cross-cell
    ...% correlation matrices, whose constant unit diagonals have no R-hat, so the
    ...% list is 34's with b_nu -> b_sigma. There is no delegation twin: analysis 9
    ...% is selected by diffest==3 && sserrvar==1 and has neither a subject-level
    ...% nor a concurrent variant.
    52, {'ic_gamma_ls_dodiff'}; ...                        % analysis 9  (gamma, LS)
    ...% The TWO-FACET DYNAMIC location-scale designs (same ruling), which open the
    ...% dynamic family. Codes 36/44 name 'Intercept_nu'/'b_nu', which prefix
    ...% nothing in these goldens. These lists are 45/46's with the five crossed
    ...% log-mean SDs added - equivalently 36/44's with the two nu names swapped for
    ...% their sigma counterparts. 54's golden is byte-identical to 53's by the same
    ...% shared-case-block delegation that makes 46's identical to 45's.
    53, {'dynrel_gamma_ls_trt_onefacet'}; ...              % analysis 14 (gamma, LS)
    54, {'dynrel_gamma_ls_sserr_trt_onefacet'}; ...        % analysis 27 (gamma, LS)
    ...% The DYNAMIC DIFFERENCE location-scale design (increment 8b), which closes
    ...% the one-facet gamma LS family. Code 41 names 'b_nu'/'b_nu_dim', which
    ...% prefix nothing in this golden - the same silent-narrowing failure codes
    ...% 45-54 exist to avoid. This list is 41's with the dispersion pair swapped
    ...% for their sigma counterparts, equivalently 49's (static analysis 7 LS)
    ...% with the two dimension-slope blocks added.
    55, {'dynrel_gamma_ls_diff_onefacet'}; ...             % analysis 12 (gamma, LS)
    ...% The SUBJECT-LEVEL dynamic difference (increment 8c). Code 56's list is 55's
    ...% VERBATIM, because the emitted model is 55's verbatim - the builder is a
    ...% one-line delegation and the golden is byte-identical, pinned by
    ...% TestStanModelSnapshots. The per-subject effects the extraction reads (r_id)
    ...% are transformed parameters, deliberately not listed, exactly as for 46<-45
    ...% and 50<-49. The separate code exists so a run's convergence record names the
    ...% analysis it belongs to, NOT to close a coverage hole - 55 would cover every
    ...% label here, since it is the same label set.
    ...% NOTE this is the only entry in the chain with no log-nu counterpart:
    ...% analysis 13 is location-scale only. (Codes 57/58 below exist, but they
    ...% belong to the modular cut-copula designs, not to a log-nu analysis 13.)
    56, {'dynrel_gamma_ls_diff_sserr_onefacet'}; ...       % analysis 13 (gamma, LS)
    ...% The MODULAR CUT-COPULA stage-1 margins designs (analyses 19/20 under
    ...% family='gamma'). 57's list is 56's verbatim - same free-parameter
    ...% families - and 58 adds the four extra crossed facet SDs and their
    ...% cross-event correlations.
    ...%
    ...% WHAT IS DELIBERATELY ABSENT, and why it is not a coverage hole. The
    ...% concurrent gamma codes 38/39/40 all list 'rho_e', because in those
    ...% FULL-JOINT models the copula correlation is a sampled parameter. Here it
    ...% is not: stage 1 fits the margins only and rho_e is estimated afterwards
    ...% on a grid, conditionally on each retained draw. Listing it would name a
    ...% parameter the model never declares, which is exactly the dead-prefix
    ...% failure this test exists to catch. The stage-2 posterior carries its own
    ...% diagnostics in REL.out.copula instead.
    57, {'dynrel_gamma_ls_diff_sserr_rescor_onefacet'}; ...     % analysis 19 (gamma, LS)
    58, {'dynrel_gamma_ls_diff_sserr_rescor_trt_onefacet'}; ... % analysis 20 (gamma, LS)
    ...% The person-specific dynamic NONCONCURRENT difference-of-differences
    ...% (analysis 28, gamma LS arm). Four margins in one model: the same
    ...% free-parameter families as 57 at four cells instead of two, with the
    ...% FULL correlation matrices (cor_id 8x8, cor_trl 4x4) listed because
    ...% stage 3 consumes every off-diagonal. rho_e is absent by DEFINITION,
    ...% not by stage separation: nothing estimates residual coupling anywhere
    ...% in this design (the six residual covariances are fixed to zero).
    59, {'dynrel_gamma_ls_dodiff_sserr_onefacet'}; ...          % analysis 28 (gamma, LS)
    ...% Code 59 plus the four crossed facet blocks (SDs and full 4x4
    ...% correlation matrices), the same 57->58 relationship at four cells.
    60, {'dynrel_gamma_ls_dodiff_sserr_trt_onefacet'}; ...      % analysis 29 (gamma, LS)
    ...% The GAUSSIAN DoD arms (owner ruling 2026-08-20). The monitored lists
    ...% are identical to 59/60 - the builders share every parameter name by
    ...% construction - but the codes are separate because psyrat_convprefix
    ...% records WHICH design was fit, and reusing 59/60 would erase the
    ...% family from the convergence record.
    61, {'dynrel_dodiff_sserr_onefacet'}; ...                   % analysis 28 (Gaussian)
    62, {'dynrel_dodiff_sserr_trt_onefacet'}};                  % analysis 29 (Gaussian)

m = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:size(pairs, 1)
    m(pairs{k, 1}) = pairs{k, 2};
end
end

% -------------------------------------------------------------------------
% Prefixes that legitimately (or, where flagged, knowingly) match nothing in a
% given golden. Each row is {trt, goldenName, {prefixes}}.
% -------------------------------------------------------------------------
function e = localKnownUnmatched()
% CONDITIONAL, by design. Codes 4 and 6 serve both the non-concurrent and the
% concurrent (diffrescor=2) form of their analysis. Lrescor / rescor /
% err_varcov exist only in the concurrent model, so they are dead against the
% non-concurrent goldens and alive against the *_rescor_* ones. One shared list
% covering both forms is the intended arrangement.
rescorOnly = {'Lrescor', 'rescor', 'err_varcov'};

% CONDITIONAL, by design (S15, 2026-08-02). Codes 3 and 4 each serve two
% families whose builders declare different names, so closing the S15 coverage
% gaps required listing BOTH spellings and accepting that each is dead against
% the other family's goldens. Same shape as rescorOnly above.
%
% dimSlopesOnly: the analysis-11 dynamic model declares the dimension slopes
% b / b_sigma; the four analysis-6 ic_sserr_* models that share code 3 do not.
% trialSdSpellings: the analysis-7 difference builders declare sd_trl, the
% analysis-9 DoD builder declares sd_trial, and neither prefix-matches the
% other (strncmp('sd_trial','sd_trl',6) is false). Each row below names the
% spelling that is dead in that golden.
dimSlopesOnly = {'b', 'b_sigma'};

% OPEN FINDING (reported, not fixed - changing it would change what gets
% reported for a design, which is maintainer-adjudicated work):
%
% (a) code 4, 'Intercept' and 'sigma_b': dead in EVERY golden code 4 serves.
%     No case-7 or case-9 Gaussian model declares either name.
% (b) code 6 / dynrel_diff_*_onefacet: analyses 12 and 13 are ONE-facet but
%     reuse code 6, which was written for the two-facet analysis 10. The
%     occasion-related entries are dead there.
deadInAllOf4 = {'Intercept', 'sigma_b'};
twoFacetOnly = {'sd_occ', 'sd_tid', 'sd_oid', 'sd_to', ...
    'occ_varcov', 'tid_varcov', 'oid_varcov', 'to_varcov'};

e = { ...
    3, 'ic_sserr_base',              dimSlopesOnly; ...
    3, 'ic_sserr_group',             dimSlopesOnly; ...
    3, 'ic_sserr_event',             dimSlopesOnly; ...
    3, 'ic_sserr_group_event',       dimSlopesOnly; ...
    4, 'ic_diff_event',              [deadInAllOf4, rescorOnly, {'sd_trial'}]; ...
    4, 'ic_diff_event_priors',       [deadInAllOf4, rescorOnly, {'sd_trial'}]; ...
    4, 'ic_diff_group_event',        [deadInAllOf4, rescorOnly, {'sd_trial'}]; ...
    4, 'ic_diff_rescor_event',       [deadInAllOf4, {'sd_trial'}]; ...
    4, 'dod_event',                  [deadInAllOf4, rescorOnly, {'sd_trl'}]; ...
    4, 'dod_event_priors',           [deadInAllOf4, rescorOnly, {'sd_trl'}]; ...
    6, 'trt_diff_event_time',        rescorOnly; ...
    6, 'trt_diff_event_time_priors', rescorOnly; ...
    6, 'dynrel_diff_onefacet',       [twoFacetOnly, rescorOnly]; ...
    6, 'dynrel_diff_sserr_onefacet', [twoFacetOnly, rescorOnly]};
end

% -------------------------------------------------------------------------
% Scalar `real` free parameters knowingly left out of a branch's list. Each row
% is {trt, goldenName, {parameterNames}}. Every entry here is an OPEN FINDING:
% the named component is estimated but never convergence-checked.
% -------------------------------------------------------------------------
function e = localKnownCoverageGaps()
% EMPTY, and that is the assertion. Every scalar `real` free parameter in every
% mapped golden is now matched by some prefix in its branch's list.
%
% It has not always been empty. S15 (2026-08-02) filed two gaps that lived here
% until the maintainer ruled to close them:
%   - code 3, 'sig_trl': the trial main-effect SD (sigma_i), declared by all
%     five code-3 goldens but missing from the list since S10 added it.
%   - code 5, 'sd_rescor': the SD of the per-person residual-correlation
%     hierarchy in the analysis-8 concurrent model. The list carried 'rescor',
%     which prefix-matches 'rescor_mu' but not 'sd_rescor'.
% Both were closed in psyrat_convprefix rather than documented further, so the
% rows were deleted rather than annotated.
%
% Keep the {trt, goldenName, {parameterNames}} shape when adding a row: a new
% gap belongs here WITH a stated reason, so that leaving a component
% unmonitored is a deliberate, reviewable act rather than a silent one.
e = cell(0, 3);
end

function names = localLookupExemption(tbl, trt, goldenName)
names = {};
for k = 1:size(tbl, 1)
    if tbl{k, 1} == trt && strcmp(tbl{k, 2}, goldenName)
        names = tbl{k, 3};
        return;
    end
end
end

% -------------------------------------------------------------------------
% Source parsing: pull the literal prefix lists out of psyrat_convprefix.
% psyrat_convprefix is a local function of psyrat_computevarcomp.m, so it is not
% reachable by a normal call; reading the source is the only way to pin it.
% -------------------------------------------------------------------------
function lists = localParseConvPrefix(testCase)
srcPath = fullfile(testCase.projectRoot(), 'subroutines', 'estimation', ...
    'psyrat_computevarcomp.m');
testCase.assertTrue(isfile(srcPath), ...
    sprintf('Cannot find psyrat_computevarcomp.m at %s', srcPath));

lines = regexp(fileread(srcPath), '\r\n|\n|\r', 'split');

startIdx = find(~cellfun(@isempty, regexp(lines, ...
    '^\s*function\s+inp2check\s*=\s*psyrat_convprefix\(', 'once')), 1);
testCase.assertNotEmpty(startIdx, ...
    'psyrat_convprefix was not found in psyrat_computevarcomp.m.');

stopIdx = startIdx + find(~cellfun(@isempty, regexp( ...
    lines(startIdx+1:end), '^\s*function\s', 'once')), 1) - 1;
if isempty(stopIdx)
    stopIdx = numel(lines);
end

lists = containers.Map('KeyType', 'double', 'ValueType', 'any');
currentTrt = [];
buffer = '';
for k = startIdx+1:stopIdx
    ln = strtrim(lines{k});
    if isempty(ln) || startsWith(ln, '%')
        continue;   % blank line or comment
    end

    tok = regexp(ln, '^(?:else)?if\s+trt\s*==\s*(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        currentTrt = str2double(tok{1});
        continue;
    end

    if startsWith(ln, 'inp2check')
        buffer = ln;
    elseif ~isempty(buffer)
        buffer = [buffer ' ' ln]; %#ok<AGROW>
    else
        continue;
    end
    buffer = strrep(buffer, '...', ' ');

    if contains(buffer, '};')
        if ~isempty(currentTrt)
            % The trailing `else` arm has no trt code and is skipped.
            quoted = regexp(buffer, '''([^'']*)''', 'tokens');
            lists(currentTrt) = cellfun(@(c) c{1}, quoted, ...
                'UniformOutput', false);
        end
        currentTrt = [];
        buffer = '';
    end
end

testCase.assertNotEmpty(keys(lists), ...
    'Parsed no prefix lists out of psyrat_convprefix - the parser is stale.');
end

% -------------------------------------------------------------------------
% Stan golden parsing.
% -------------------------------------------------------------------------
function names = localDeclaredNames(testCase, goldenName)
% Every name declared at the top level of the blocks stansummary reports on:
% sampled parameters, transformed parameters, and generated quantities.
src = localReadGolden(testCase, goldenName);
names = {};
for blk = {'parameters', 'transformed parameters', 'generated quantities'}
    body = localBlockBody(src, blk{1});
    stmts = localTopLevelStatements(body);
    for s = 1:numel(stmts)
        nm = localDeclaredName(stmts{s});
        if ~isempty(nm)
            names{end+1} = nm; %#ok<AGROW>
        end
    end
end
end

function names = localCoefficientBlockParameters(testCase, goldenName)
% Free parameters that are population-level effect (slope) blocks, matched by
% EXACT name against the set the builders use. An exact-name set rather than a
% b* pattern on purpose: dod_event declares b_cell / b_sigma_cell and analysis
% 21/22 declare per-subject blocks, and sweeping those in would turn this test
% into an assertion about naming conventions instead of about coverage.
coeffNames = {'b', 'b_sigma', 'b_nu', 'b_dim', 'b_sigma_dim', 'b_nu_dim'};
src = localReadGolden(testCase, goldenName);
stmts = localTopLevelStatements(localBlockBody(src, 'parameters'));
names = {};
for s = 1:numel(stmts)
    nm = localDeclaredName(stmts{s});
    if ~isempty(nm) && any(strcmp(nm, coeffNames))
        names{end+1} = nm; %#ok<AGROW>
    end
end
end

function names = localScalarRealParameters(testCase, goldenName)
% Free parameters declared as a bare scalar `real` (constraints allowed).
src = localReadGolden(testCase, goldenName);
stmts = localTopLevelStatements(localBlockBody(src, 'parameters'));
names = {};
for s = 1:numel(stmts)
    stripped = strtrim(regexprep(stmts{s}, '<[^<>]*>', ''));
    tok = regexp(stripped, '^real\s+([A-Za-z_]\w*)$', 'tokens', 'once');
    if ~isempty(tok)
        names{end+1} = tok{1}; %#ok<AGROW>
    end
end
end

function src = localReadGolden(testCase, goldenName)
p = fullfile(testCase.projectRoot(), 'tests', 'golden', 'stan', ...
    [goldenName '.stan']);
testCase.assertTrue(isfile(p), sprintf('Missing Stan golden: %s', p));
src = fileread(p);
% Strip comments before any brace matching - a comment may contain braces.
src = regexprep(src, '(?s)/\*.*?\*/', '');
src = regexprep(src, '//[^\r\n]*', '');
end

function body = localBlockBody(src, blockName)
% Text between the braces of a top-level Stan block ('' when absent).
body = '';
pat = ['(?m)^' regexprep(blockName, '\s+', '\\s+') '\s*\{'];
openIdx = regexp(src, pat, 'end', 'once');
if isempty(openIdx)
    return;
end
tail = src(openIdx:end);
depth = cumsum((tail == '{') - (tail == '}'));
closeIdx = find(depth == 0, 1);
if isempty(closeIdx)
    return;
end
body = tail(2:closeIdx-1);
end

function stmts = localTopLevelStatements(body)
% Split on semicolons that sit outside any nested braces (loops in transformed
% parameters / generated quantities), so only declarations survive intact.
stmts = {};
if isempty(body)
    return;
end
depth = cumsum((body == '{') - (body == '}'));
cuts = find(body == ';' & depth == 0);
prev = 1;
for k = 1:numel(cuts)
    stmts{end+1} = strtrim(body(prev:cuts(k)-1)); %#ok<AGROW>
    prev = cuts(k) + 1;
end
end

function nm = localDeclaredName(stmt)
% Variable name declared by a Stan declaration statement ('' if not one).
nm = '';
s = regexprep(stmt, '<[^<>]*>', '');    % drop constraint clauses (they hold '=')
eq = strfind(s, '=');
if ~isempty(eq)
    s = s(1:eq(1)-1);                   % drop any initializer
end
s = regexprep(s, '^\s*array\s*\[[^\]]*\]\s*', '');   % leading array qualifier
types = ['real|vector|row_vector|matrix|simplex|ordered|positive_ordered|' ...
    'unit_vector|cholesky_factor_corr|cholesky_factor_cov|corr_matrix|' ...
    'cov_matrix|int|complex'];
tok = regexp(s, ['^\s*(?:' types ')\s*(?:\[[^\]]*\])?\s*([A-Za-z_]\w*)'], ...
    'tokens', 'once');
if ~isempty(tok)
    nm = tok{1};
end
end
