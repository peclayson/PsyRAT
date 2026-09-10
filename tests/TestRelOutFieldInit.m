classdef TestRelOutFieldInit < PsyRATTestBase
    % Guards psyrat_computevarcomp against growing a REL.out field that the
    % matching init branch never created.
    %
    % WHY THIS EXISTS. Every analysis case in psyrat_computevarcomp initializes
    % its REL.out slots once (inside `if i == 1`, after `REL.out = []`) and then
    % appends one column per stratum with `(:,end+1)` or `{end+1}`. MATLAB cannot
    % grow a field that does not exist yet: `s.newfield(:,end+1) = v` throws
    % MATLAB:nonExistentField rather than creating it. Verified in R2025b:
    %
    %   o = struct(); o.pop_lognu = [];
    %   o.pop_sdlog(:,end+1) = [1;2;3]
    %   -> Unrecognized field name "pop_sdlog".
    %
    % So an init branch and an extraction branch that disagree on the slot names
    % is a hard crash, and it lands AFTER the first fit has run - the most
    % expensive possible place to discover it. That is exactly what happened when
    % `gammascale = 2` (the gamma location-scale dynrel) was added: the case
    % {11,26} init block still branched on the family alone and created
    % pop_lognu/b_nu/ind_nu, while the new gammascale==2 extraction branch wrote
    % pop_sdlog/b_sigma/ind_sdlog.
    %
    % Nothing in the offline lane could catch it. The extraction block only runs
    % against a real CmdStan fit, so unit tests, the accuracy suite and the Stan
    % snapshot goldens were all green while a gammascale=2 run could not
    % complete. This test reads the source instead of executing it, so it reports
    % the same answer with or without CmdStan installed.
    %
    % WHY BRANCH-AWARE. A scanner that simply unions every init in a case block
    % does NOT catch this: case {11,26} also carries a Gaussian branch that
    % initializes pop_sdlog/b_sigma/ind_sdlog, so the missing slots look present.
    % Measured against the pre-fix source, the coarse form reported 0 problems
    % across all 25 case blocks. The scanner therefore pairs an init with a grow
    % only when both sit under the SAME family-guard text (or the init is
    % unconditional within the case block). Against that same pre-fix source the
    % branch-aware form reports all nine slots, pop_sdlog included.
    %
    % SCOPE - what this proves and what it does not. It pins WIRING: that each
    % branch creates the slots it later fills. It says nothing about whether the
    % values in those slots are right; only a live CmdStan recovery run can
    % establish that.
    %
    % Three deliberate gaps, none of which should be read as coverage:
    %   1. It keys on guards mentioning REL.family, because that is where the
    %      parameterization branching lives. A future dispatch on some other
    %      condition needs its own key added to FamilyGuardToken below.
    %   2. It scans only the analysis switch. Appends past the switch's `end`
    %      (the case-8 psyrat_compute_diff_sserr helpers, which carry their own
    %      init blocks) are not examined.
    %   3. It resolves field names STATICALLY, so appends written through a
    %      dynamic field name - REL.out.(measname)(:,end+1) - are skipped
    %      entirely. testGrowAndInitMatchersIgnoreLookalikes pins that they are
    %      skipped rather than silently matched to the wrong field, so this is a
    %      known hole and not an accident; a bug introduced at one of those sites
    %      would not be caught here.
    %
    % All three gaps are in the permissive direction: they can hide a real defect,
    % never invent one. That is the right way round for a source scanner, but it
    % means a green run is not a proof that every append is wired.

    properties (Constant)
        % Substring identifying the if/elseif chains that select a likelihood
        % family or parameterization. Chains opening on any other condition are
        % treated as unguarded, which is the permissive direction: it can produce
        % a missed finding, never a false failure.
        FamilyGuardToken = 'REL.family';
    end

    methods (Test)

        function testEveryRelOutGrowHasAMatchingBranchInit(testCase)
            % The real scan over the shipped estimator.
            root = testCase.projectRoot();
            file = fullfile(root, 'subroutines', 'estimation', ...
                'psyrat_computevarcomp.m');
            testCase.assertTrue(isfile(file), ...
                'psyrat_computevarcomp.m not found - the scan would pass vacuously.');

            lines = localReadLines(file);
            [problems, stats] = localScanSource(lines, ...
                TestRelOutFieldInit.FamilyGuardToken);

            % Sanity checks on the scanner's own reach. Without these a regex
            % typo that matched nothing would report a clean tree, which is the
            % failure mode this whole file exists to prevent.
            testCase.assertGreaterThanOrEqual(stats.nCases, 20, ...
                ['Found implausibly few top-level case blocks. The case-boundary ', ...
                'matcher is probably broken, which would make this scan vacuous.']);
            testCase.assertGreaterThanOrEqual(stats.nGrows, 100, ...
                ['Found implausibly few REL.out grow sites. The grow matcher is ', ...
                'probably broken, which would make this scan vacuous.']);
            testCase.assertGreaterThanOrEqual(stats.nFamilyChains, 2, ...
                ['Found no family-guarded if/elseif chains. The chain finder is ', ...
                'probably broken, so every branch would look unguarded and the ', ...
                'mutually-exclusive-branch case would go unchecked.']);

            testCase.verifyEmpty(problems, localFormatReport(problems));
        end

        function testScannerFlagsAnInitInAMutuallyExclusiveBranch(testCase)
            % NON-VACUITY PIN. The fixture reproduces the exact defect shape:
            % a family chain whose init branches and extraction branches
            % disagree, with the missing slot present in a SIBLING branch. A
            % scanner that unions inits across the case block reports this
            % clean; this one must not.
            broken = localFixture(false);
            [problems, stats] = localScanSource(broken, ...
                TestRelOutFieldInit.FamilyGuardToken);

            testCase.verifyEqual(stats.nFamilyChains, 2, ...
                'Fixture must expose exactly two family-guarded chains.');
            testCase.verifyNumElements(problems, 1, ...
                'Exactly the one uninitialized slot must be reported.');
            testCase.verifyEqual(problems{1}.field, 'pop_sdlog');
            testCase.verifyEqual(problems{1}.guard, 'strcmp(REL.family,''new'')');

            % The same fixture with the init branch corrected must come back
            % clean, so the finding tracks the defect and not the fixture.
            fixed = localFixture(true);
            testCase.verifyEmpty(localScanSource(fixed, ...
                TestRelOutFieldInit.FamilyGuardToken), ...
                'The corrected fixture must produce no finding.');
        end

        function testChainFinderPairsBranchesByGuardText(testCase)
            % Pins the part most likely to go subtly wrong. An `end` that closes
            % a NESTED block must not be mistaken for the chain's own `end`,
            % because that would truncate the branch ranges and silently drop
            % the grows that follow.
            lines = {
                '    case 1'
                '        if REL.family == 1'
                '            if a == 2'
                '                x = 1;'
                '            end'
                '            y = 2;'
                '        else'
                '            z = 3;'
                '        end'
                '    case 2'
                };
            chains = localFamilyChains(lines, 1, numel(lines), 'REL.family');

            testCase.verifyNumElements(chains, 1);
            testCase.verifyEqual(chains(1).guard, {'REL.family == 1', 'else'});
            % branch 1 spans the `if` through the line before `else`; branch 2
            % spans `else` through the line before the chain's own `end`.
            testCase.verifyEqual(chains(1).range, [2 6; 7 8]);
        end

        function testGrowAndInitMatchersIgnoreLookalikes(testCase)
            % The two regexes decide everything. Pinned against the forms that
            % appear in the estimator so a rewrite cannot quietly narrow them.
            grows = {
                'REL.out.mu(:,end+1) = measvalue;'                  % numeric append
                'REL.out.b(:,end+1) = {num2cell(measvalue)};'       % cell append
                'REL.out.ssinfo{end+1} = local_dynrel_ssinfo(d);'   % brace append
                'REL.out.conv.data{end+1} = psyrat_storeconv(f,3);' % nested field
                'REL.out.ntrials(:,end + 1) = median(x);'           % spaced +
                };
            for k = 1:numel(grows)
                testCase.verifyNotEmpty(localGrownFields(grows{k}), ...
                    sprintf('Grow matcher missed: %s', grows{k}));
            end

            notGrows = {
                'bsig = cell2mat(REL.out.b_sigma{:,gloc});'   % a read, not a grow
                'REL.out.labels = darray.names;'              % a plain assignment
                'x = REL.out.mu(:,end);'                      % reads the last col
                'REL.out.(measname)(:,end+1) = {v};'          % dynamic name
                };
            for k = 1:numel(notGrows)
                testCase.verifyEmpty(localGrownFields(notGrows{k}), ...
                    sprintf('Grow matcher over-matched: %s', notGrows{k}));
            end

            testCase.verifyEqual(localInitFields('REL.out.pop_lognu = [];'), ...
                {'pop_lognu'});
            testCase.verifyEqual(localInitFields('REL.out.conv.data = {};'), ...
                {'conv'});
            % A subscripted write is an append or an overwrite of one column, not
            % the initialization that makes the field exist.
            testCase.verifyEmpty(localInitFields('REL.out.mu(:,3) = v;'));
            % An equality test is not an assignment.
            testCase.verifyEmpty(localInitFields('if REL.out.ndim == 2'));
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers. Kept as functions over a cellstr of source lines so the
% scanner can be driven by the fixtures above as well as by the real file.
% ---------------------------------------------------------------------------

function lines = localReadLines(file)
%LOCALREADLINES Source file to a cellstr of comment-stripped lines.
raw = regexp(fileread(file), '\r\n|\n|\r', 'split');
lines = cellfun(@localStripComment, raw, 'UniformOutput', false);
end

function code = localStripComment(line)
%LOCALSTRIPCOMMENT Drops a trailing MATLAB line comment.
%   Cuts at the first '%' preceded by an even number of single quotes, so a '%'
%   inside a quoted string survives. Same rule as TestBaseMatlabOnly, and safe
%   in the same direction: under-stripping can only leave prose in the scan,
%   which the matchers below ignore because prose carries no `REL.out.x(` form.
code = line;
idx = strfind(line, '%');
for k = 1:numel(idx)
    if mod(count(line(1:idx(k)-1), ''''), 2) == 0
        code = line(1:idx(k)-1);
        return
    end
end
end

function fields = localGrownFields(codeLine)
%LOCALGROWNFIELDS REL.out fields appended to on this line.
%   Matches `REL.out.F[.sub]*` followed by a '(' or '{' subscript containing
%   end+1. The [^=]* keeps the subscript from running past an assignment, so a
%   later `end+1` on the right-hand side cannot be mistaken for an append.
tok = regexp(codeLine, ...
    'REL\.out\.(\w+)(?:\.\w+)*\s*[({][^=]*end\s*\+\s*1', 'tokens');
fields = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
end

function fields = localInitFields(codeLine)
%LOCALINITFIELDS REL.out fields created by a plain assignment on this line.
%   Requires the '=' to follow the field path directly, with no subscript
%   between: `REL.out.b = {}` initializes, `REL.out.b(:,end+1) = {...}` does
%   not. The trailing [^=] rejects '==' so a comparison is not read as an
%   assignment.
tok = regexp(codeLine, 'REL\.out\.(\w+)(?:\.\w+)*\s*=[^=]', 'tokens');
fields = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
end

function chains = localFamilyChains(code, lo, hi, token)
%LOCALFAMILYCHAINS if/elseif/else chains selecting a family/parameterization.
%   A chain is opened by an `if` whose condition contains TOKEN. Its branches
%   and its closing `end` are found by matching the opening `if`'s indentation,
%   which is what keeps a nested block's `end` from closing the chain early.
%   Returns a struct array with .guard (branch labels, 'else' for the else arm)
%   and .range (one [firstLine lastLine] row per branch).
chains = struct('guard', {}, 'range', {});
k = lo;
while k <= hi
    tok = regexp(code{k}, '^(\s*)if\s+(.+?)\s*$', 'tokens', 'once');
    if ~isempty(tok) && contains(tok{2}, token)
        indent = numel(tok{1});
        guard = tok(2);
        starts = k;
        stop = [];
        j = k + 1;
        while j <= hi
            lead = regexp(code{j}, '^\s*', 'match', 'once');
            if numel(lead) == indent
                body = strtrim(code{j});
                elseifTok = regexp(body, '^elseif\s+(.+?)$', 'tokens', 'once');
                if strcmp(body, 'end')
                    stop = j;
                    break
                elseif strcmp(body, 'else')
                    guard{end+1} = 'else';   %#ok<AGROW>
                    starts(end+1) = j;       %#ok<AGROW>
                elseif ~isempty(elseifTok)
                    guard{end+1} = elseifTok{1}; %#ok<AGROW>
                    starts(end+1) = j;           %#ok<AGROW>
                end
            end
            j = j + 1;
        end
        if ~isempty(stop)
            lasts = [starts(2:end) - 1, stop - 1];
            chains(end+1) = struct('guard', {guard}, ...
                'range', [starts(:), lasts(:)]); %#ok<AGROW>
            k = stop;
        end
    end
    k = k + 1;
end
end

function [problems, stats] = localScanSource(code, token)
%LOCALSCANSOURCE Reports every REL.out grow with no matching branch init.
%   Splits the source into top-level `case` blocks of the analysis switch, then
%   within each block pairs grows with inits by family-guard text.

problems = {};
stats = struct('nCases', 0, 'nGrows', 0, 'nFamilyChains', 0);

isCase = ~cellfun(@isempty, regexp(code, '^    case\s', 'once'));
caseStart = find(isCase);
if isempty(caseStart)
    return
end
% The switch closes at the first column-0 `end` after the last case. Falling
% back to the file end keeps the fixtures (which carry no switch) scannable.
tail = find(~cellfun(@isempty, regexp(code, '^end\s*$', 'once')));
switchEnd = tail(find(tail > caseStart(end), 1, 'first'));
if isempty(switchEnd)
    switchEnd = numel(code) + 1;
end
bounds = [caseStart(:); switchEnd];
stats.nCases = numel(caseStart);

for c = 1:numel(caseStart)
    lo = bounds(c);
    hi = bounds(c+1) - 1;
    chains = localFamilyChains(code, lo, hi, token);
    stats.nFamilyChains = stats.nFamilyChains + numel(chains);

    % Tag every line with the family-guard branch it sits in ('' = unguarded,
    % i.e. reached whichever family branch was taken).
    label = repmat({''}, hi - lo + 1, 1);
    for ch = 1:numel(chains)
        for br = 1:size(chains(ch).range, 1)
            span = chains(ch).range(br,:);
            label(span(1)-lo+1 : span(2)-lo+1) = chains(ch).guard(br);
        end
    end

    [initByGuard, growByGuard] = localCollect(code, lo, hi, label);

    % Fields initialized in EVERY branch of some chain are reached
    % unconditionally, so they cover an unguarded grow just as a bare init does.
    coverAlways = localGet(initByGuard, '');
    for ch = 1:numel(chains)
        common = localGet(initByGuard, chains(ch).guard{1});
        for br = 2:numel(chains(ch).guard)
            common = intersect(common, localGet(initByGuard, chains(ch).guard{br}));
        end
        coverAlways = union(coverAlways, common);
    end

    guards = growByGuard.keys;
    for q = 1:numel(guards)
        g = guards{q};
        if isempty(g)
            cover = coverAlways;
        else
            cover = union(localGet(initByGuard, ''), localGet(initByGuard, g));
        end
        entries = growByGuard(g);
        seen = {};
        for e = 1:numel(entries)
            f = entries(e).field;
            if any(strcmp(f, cover)) || any(strcmp(f, seen))
                continue
            end
            seen{end+1} = f; %#ok<AGROW>
            problems{end+1} = struct('field', f, 'guard', g, ...
                'line', entries(e).line, 'case', lo, ...
                'text', strtrim(code{entries(e).line})); %#ok<AGROW>
        end
        stats.nGrows = stats.nGrows + numel(entries);
    end
end
problems = problems(:);
end

function [initByGuard, growByGuard] = localCollect(code, lo, hi, label)
%LOCALCOLLECT Init and grow sites in one case block, keyed by guard text.
initByGuard = containers.Map('KeyType', 'char', 'ValueType', 'any');
growByGuard = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = lo:hi
    line = code{k};
    if isempty(strtrim(line))
        continue
    end
    g = label{k - lo + 1};

    grown = localGrownFields(line);
    for t = 1:numel(grown)
        entry = struct('field', grown{t}, 'line', k);
        if growByGuard.isKey(g)
            growByGuard(g) = [growByGuard(g), entry];
        else
            growByGuard(g) = entry;
        end
    end

    created = localInitFields(line);
    if ~isempty(created)
        if initByGuard.isKey(g)
            initByGuard(g) = [initByGuard(g), created];
        else
            initByGuard(g) = created;
        end
    end
end
end

function v = localGet(map, key)
%LOCALGET Unique cellstr for KEY, or an empty cellstr when KEY is absent.
if map.isKey(key)
    v = unique(map(key));
else
    v = {};
end
end

function msg = localFormatReport(problems)
%LOCALFORMATREPORT Actionable failure text: which slot, which branch, the fix.
if isempty(problems)
    msg = '';
    return
end
parts = cell(numel(problems), 1);
for k = 1:numel(problems)
    p = problems{k};
    if isempty(p.guard)
        where = 'the case block (unguarded)';
    else
        where = sprintf('the branch guarded by `%s`', p.guard);
    end
    parts{k} = sprintf(['  psyrat_computevarcomp.m:%d  REL.out.%s is appended ' ...
        'to but never created in %s\n      %s'], p.line, p.field, where, p.text);
end
msg = sprintf(['%d REL.out slot(s) appended to without a matching branch ' ...
    'initialization.\nMATLAB cannot grow a field that does not exist: ' ...
    '`s.f(:,end+1) = v` throws MATLAB:nonExistentField rather than creating ' ...
    'f. Each of these therefore crashes at run time, AFTER the first fit ' ...
    'completes.\n\n%s\n\nHow to fix: add the slot to the init block branch ' ...
    'that carries the SAME guard (inside `if i == 1`, after `REL.out = []`). ' ...
    'Adding it to a sibling branch does not help - only the branch actually ' ...
    'taken runs.'], numel(problems), strjoin(parts, sprintf('\n')));
end

function lines = localFixture(corrected)
%LOCALFIXTURE Minimal reproduction of the defect this test exists to catch.
%   Two family-guarded chains inside one case block: an init chain and an
%   extraction chain, sharing branch guards. In the broken form the guarded
%   init branch omits pop_sdlog while the SIBLING `else` branch supplies it -
%   the arrangement that defeats a scanner unioning inits across the block.
if corrected
    guardedInit = '                    REL.out.pop_sdlog = [];';
else
    guardedInit = '                    REL.out.pop_lognu = [];';
end
lines = {
    '    case {11, 26}'
    '        for i = 1:n'
    '            if i == 1'
    '                REL.out = [];'
    '                if strcmp(REL.family,''new'')'
    guardedInit
    '                    REL.out.mu = [];'
    '                else'
    '                    REL.out.pop_sdlog = [];'
    '                    REL.out.mu = [];'
    '                end'
    '            end'
    '            if strcmp(REL.family,''new'')'
    '                REL.out.mu(:,end+1) = a;'
    '                REL.out.pop_sdlog(:,end+1) = b;'
    '            else'
    '                REL.out.mu(:,end+1) = a;'
    '                REL.out.pop_sdlog(:,end+1) = b;'
    '            end'
    '        end'
    '    case 12'
    };
end
