classdef LS27WordingRules
    %LS27WORDINGRULES Surfaces and content rules for the RESTRICTED disclosure
    %of the three gamma two-facet designs (analyses 14, 25 and 27; LS27 wording
    %residual) and for the dynamic-migration wording of
    %psyrat_reprocess_gamma_result (Q8 wording residual). Used by
    %TestRestrictedDisclosureWording.
    %
    %WHAT A RULE IS. Each rule checks what a sentence SAYS, never how it is
    %worded as a whole: no rule compares a paragraph with a stored copy. A rule
    %looks for a required fact (the analysis-27 dependability caveat, the unit
    %of a ratio, the number of captures behind a range) or a forbidden claim
    %(an unscoped universal, a certainty the evidence does not carry). Where a
    %rule needs a word list, such as the ways of saying that intervals do not
    %remedy the evidence, the list lives in one constant below. A correct
    %rewording that a rule rejects is fixed by widening that list, never by
    %weakening what the rule requires.
    %
    %SCOPE OF EACH CHECK. Most facts attach to a sentence or a clause, not to a
    %whole surface. The clearest case is the analysis-27 dependability (D)
    %caveat: the export header already says "dependability" in its clause about
    %analyses 14 and 25, so a check that the header merely CONTAINS
    %"dependab" passes on the text that omits the caveat. The rule therefore
    %asks for one clause that names analysis 27, both coefficients (G and D)
    %and a downward bias, and clauses are split at sentence ends and at
    %semicolons.
    %
    %EVERY RULE RETURNS VIOLATIONS, not a Boolean: a cell row of messages that
    %name the surface, the rule and the offending text, so one run lists every
    %place that still needs a change.
    %
    %THE SURFACES. The one text source (psyrat_gamma_twofacet_restriction: its
    %export header, viewer note, saved-figure note, Analyze notice and help
    %text), the CHANGELOG entry and design-maturity row, the Restricted text of
    %manual Section 5.4, the availability matrix's footnote [t] and section 5
    %paragraph, and, for the scope rules only, the Restricted sentences of
    %manual chapters 9 and 15. Each surface carries flags that say which
    %conditional rules apply to it (a short note is not asked for
    %dataset-specific offsets). Chapter 16's status cards were read and need no
    %rule; chapter 14 (scripted use) belongs to another increment.
    %
    %THE FACTS BEHIND THE RULES (row 27 contract, 2026-09-22, and its review).
    %The ranges O 2.4-12.1, R 2.2-7.2 and E 1.1-1.9 are posterior mean over
    %planted value, on the variance scale, over the five scope captures; the
    %retained analysis-27 location-scale control reads O 1.42, R 1.66, E 1.02
    %and is not a pass. Analysis 27's offsets are medians, across participants,
    %of the signed error of each participant's posterior-mean coefficient
    %(posterior mean minus planted value): D -.1304 in the converged
    %prospective fit A, -.097 in the log-nu capture, -.023 in the retained
    %control; B's +.0766 (printed .08) comes from a fit that did not converge.
    %The "within about .03 and .06" observation belongs to analyses 14 and 25
    %only. The 2e-12 agreement between production and an independent
    %enumeration is implementation consistency, not recovery, identifiability
    %or sampling validation. No universal maximum is claimed.
    %
    %SEEDED CONTROLS. A forbid-only rule is green when its subject is absent,
    %so every rule is also run on seeded bad and good sentences
    %(TestRestrictedDisclosureWording): the bad seed must be caught and its
    %correction must pass.

    properties (Constant)
        %an unscoped statement about every in-repo check or capture
        UNIVERSALS = { ...
            '\<(every|all|any)\s+in-repo\s+(checks?|captures?|tests?|fits?)\>', ...
            ['\<no\s+in-repo\s+(checks?|captures?|tests?|fits?)\>' ...
            '(?!\s+(of|for|on)\s+(these|the three|the restricted|analyses 14|analysis 14))'], ...
            '\<every check the toolbox.s own validation'};
        %the definite article turns a set of designs into all of them
        DEFINITE = '\<the (gamma )?two-facet (designs|arms)\>(?! of analys)';
        %a stated upper bound on an offset or ratio
        MAXIMUM = '\<(up to|at most|no more than|as much as|as large as|a maximum of)\s+[-+]?\d*\.\d';
        %the words that say a ratio is on the variance scale
        VARSCALE = ['variance[- ]scale|variance[- ](component )?ratios?|' ...
            'ratios? of (the )?(posterior[- ]mean )?variances'];
        %a ratio statement ("2.4 to 12.1 times")
        RATIO = '\d\s+times\>';
        %an affirmative pass or recovery claim (checked for a negation before it)
        AFFIRM = '\<(pass|passes|passed|recovered|recovers|at level|within tolerance|meets)\>';
        NEGATION = '\<(not|no|never|nor|without)\>';
        %ways of saying that intervals do not repair the evidence
        REMEDY = ['\<(not|no|cannot|nor|never)\>.{0,80}\<(remed|repair|rescu|compensat|' ...
            'fix|cure|resolv|substitut|offset|remov|overcom|correct)'];
        %the interval-miss statement that must stay
        MISS = ['\<intervals?\>.{0,60}(exclud|miss|above the truth|outside|not cover)|' ...
            '(exclud|miss)\w*.{0,30}\<truth\>'];
        %the 2e-12 agreement framed as consistency, with recovery negated
        CONSISTENCY = 'consisten|implement';
        NOTRECOVERY = '\<(not|no|neither|nor)\>.{0,120}(recover|identif|sampl|validat)';
        %causal or "not a formula error" certainty
        CERTAINTY = {'not a formula error', 'not at fault', ...
            '\<(this|it|the reason) is a (summary-and-design )?property\>', ...
            'is a summary-and-design property'};
        HEDGE = ['consistent with|likel|attribut|suggest|appear|plausib|probab|' ...
            '\<may\>|\<might\>|explained by|interpret'];
        %Q8: "no number changes" and its variants, which do not say WHAT is unchanged
        AMBIGUOUS = '\<(no numeric changes?|nothing numeric changes|no numbers? changes?)\>';
        %Q8: the raw draws are unchanged
        RAWUNCHANGED = ['\<raw\>( posterior)? draws?[^.;]{0,60}(unchanged|not changed|' ...
            'not modified|do not change|does not change|untouched|stay the same)|' ...
            '(unchanged|not changed)[^.;]{0,30}\<raw\>'];
    end

    methods (Static)

        % ------------------------------------------------------------------
        % Surfaces
        % ------------------------------------------------------------------

        function S = surfaces(root)
            %Every RESTRICTED-disclosure surface: the text source, then the
            %shipped documents under ROOT (the repository root).
            S = [LS27WordingRules.sourceSurfaces(), LS27WordingRules.documentSurfaces(root)];
        end

        function S = sourceSurfaces()
            %The one text source. The Analyze notice is read twice: the
            %all-designs form and the analysis-27 form (the notice for 14 or
            %25 may legitimately drop the analysis-27 sentence).
            R = psyrat_gamma_twofacet_restriction();
            R27 = psyrat_gamma_twofacet_restriction(27);
            S = [ ...
                LS27WordingRules.surface('export header (R.header)', R.header, ...
                'three', true, 'coefficients', true, 'long', true, 'offsets', true, ...
                'remedy', true, 'intervalMiss', true), ...
                LS27WordingRules.surface('viewer note (R.screen)', R.screen, ...
                'coefficients', true), ...
                LS27WordingRules.surface('saved-figure note (R.figure)', ...
                LS27WordingRules.fieldOr(R, 'figure'), 'coefficients', true), ...
                LS27WordingRules.surface('Analyze notice, all designs (R.dialog)', ...
                R.dialog, 'coefficients', true), ...
                LS27WordingRules.surface('Analyze notice, analysis 27 (R.dialog)', ...
                R27.dialog, 'coefficients', true), ...
                LS27WordingRules.surface('help text of psyrat_gamma_twofacet_restriction', ...
                help('psyrat_gamma_twofacet_restriction'), ...
                'coefficients', true, 'long', true, 'offsets', true, 'intervalMiss', true)];
        end

        function S = documentSurfaces(root)
            %The shipped documents that restate the finding. Each block is
            %found by a content anchor, not a line number, so an edit above it
            %does not move it; a missing anchor is reported as a violation.
            cl = fileread(fullfile(root, 'CHANGELOG.md'));
            ma = fileread(fullfile(root, 'documentation', 'model_availability_matrix.md'));
            m5 = fileread(fullfile(root, 'documentation', 'manual', '05_choosing-an-analysis.md'));
            m9 = fileread(fullfile(root, 'documentation', 'manual', '09_viewing-results.md'));
            m15 = fileread(fullfile(root, 'documentation', 'manual', '15_advanced-designs.md'));
            S = [ ...
                LS27WordingRules.surface('CHANGELOG 0.1.1-beta Restricted entry', ...
                LS27WordingRules.block(cl, '^- ', {'Restricted', '12\.1'}), ...
                'coefficients', true, 'long', true, 'intervalMiss', true), ...
                LS27WordingRules.surface('CHANGELOG design-maturity Restricted row', ...
                LS27WordingRules.lineStarting(cl, '| **Restricted**'), ...
                'coefficients', true, 'long', true, 'remedy', true, 'intervalMiss', true), ...
                LS27WordingRules.surface('manual Section 5.4 Restricted text', ...
                LS27WordingRules.manual54(m5), ...
                'coefficients', true, 'long', true, 'offsets', true, 'remedy', true, ...
                'intervalMiss', true), ...
                LS27WordingRules.surface('availability matrix footnote [t]', ...
                LS27WordingRules.lineStarting(ma, '- **[t]**')), ...
                LS27WordingRules.surface('availability matrix section 5 RESTRICTED paragraph', ...
                LS27WordingRules.block(ma, '^\*\*RESTRICTED', {'RESTRICTED'}), ...
                'coefficients', true, 'long', true, 'intervalMiss', true), ...
                LS27WordingRules.surface('manual chapter 9 viewer-note paragraph (scope rules)', ...
                LS27WordingRules.block(m9, '.', {'Restricted two-facet designs'})), ...
                LS27WordingRules.surface('manual chapter 15 Restricted paragraphs (scope rules)', ...
                LS27WordingRules.paragraphsWith(m15, 'Restricted'))];
        end

        function t = fieldOr(R, f)
            %R.(f) when the text source has that field, else '' (reported as a
            %missing anchor by every rule).
            t = '';
            if isfield(R, f); t = R.(f); end
        end

        function s = surface(name, text, varargin)
            %One surface: its name, its normalized text and the flags that
            %switch on the conditional rules. found is false when the anchor
            %did not match (text is then empty).
            s = struct('name', name, 'text', '', 'found', false, 'three', false, ...
                'coefficients', false, 'long', false, 'offsets', false, ...
                'remedy', false, 'intervalMiss', false);
            s.text = LS27WordingRules.normalize(text);
            s.found = ~isempty(s.text);
            for k = 1:2:numel(varargin)
                s.(varargin{k}) = varargin{k + 1};
            end
        end

        % ------------------------------------------------------------------
        % The five corrections, one entry point each
        % ------------------------------------------------------------------

        function v = scopeWording(s)
            %(1) Name the three designs; no unscoped universal about every
            %in-repo check or capture (analysis 5 is fitted at five occasions
            %in the repository); no "the gamma two-facet designs".
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v); return; end
            t = s.text;
            if s.three && ~LS27WordingRules.has(t, '\<three gamma two-facet designs\>')
                v{end + 1} = sprintf('%s: does not say "three gamma two-facet designs (14, 25 and 27)"', s.name);
            end
            v = [v, LS27WordingRules.matches(s, LS27WordingRules.DEFINITE, ...
                'reads as every gamma two-facet design; name the three')];
            for k = 1:numel(LS27WordingRules.UNIVERSALS)
                v = [v, LS27WordingRules.matches(s, LS27WordingRules.UNIVERSALS{k}, ...
                    'universal over in-repo checks; limit it to the reviewed captures of these designs')]; %#ok<AGROW>
            end
        end

        function v = analysis27Caveat(s)
            %(2) Analysis 27's D caveat beside its G caveat, in one clause;
            %the "within about .03 and .06" observation sits in a clause that
            %names analyses 14 and 25 and not 27 (a comma cannot carry it over
            %to 27); long surfaces also give the dataset-specific D offsets; a
            %surface that prints analysis 27's D offsets names the statistic
            %(medians across participants of the posterior-mean coefficient's
            %signed error); B's positive D offset, in any rounding, only beside
            %its non-convergence.
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v) || ~s.coefficients; return; end
            c = LS27WordingRules.clauses(s.text);
            sen = LS27WordingRules.sentences(s.text);
            ok = false;
            for k = 1:numel(c)
                if LS27WordingRules.states27Bias(c{k}); ok = true; break; end
            end
            if ~ok
                v{end + 1} = sprintf(['%s: no clause names analysis 27 with both G and D ' ...
                    'and a downward bias (the D caveat is missing)'], s.name);
            end
            for k = 1:numel(c)
                isWithin = LS27WordingRules.has(c{k}, '\<within\>') && ...
                    (LS27WordingRules.gTerm(c{k}) || LS27WordingRules.dTerm(c{k}));
                if isWithin && LS27WordingRules.mentions(c{k}, 27)
                    v{end + 1} = sprintf(['%s: analysis 27 shares the clause of the "within" ' ...
                        'observation, which belongs to 14 and 25 only: "%s"'], ...
                        s.name, LS27WordingRules.clip(c{k})); %#ok<AGROW>
                end
                if isWithin && LS27WordingRules.hasCase(c{k}, '(^|[^\d])\.03(?!\d)') && ...
                        ~(LS27WordingRules.mentions(c{k}, 14) && LS27WordingRules.mentions(c{k}, 25))
                    v{end + 1} = sprintf(['%s: the "within about .03 and .06" observation does ' ...
                        'not name analyses 14 and 25: "%s"'], s.name, LS27WordingRules.clip(c{k})); %#ok<AGROW>
                end
            end
            dOffset = cellfun(@(x) LS27WordingRules.dTerm(x) && ...
                LS27WordingRules.hasCase(x, '(^|[^\d])-?\.13\d*(?!\d)'), sen);
            named = cellfun(@(x) LS27WordingRules.has(x, '\<median') && ...
                LS27WordingRules.has(x, 'participant') && ...
                LS27WordingRules.has(x, 'posterior[- ]mean'), sen);
            if any(dOffset) && ~any(named)
                v{end + 1} = sprintf(['%s: prints analysis 27''s D offsets without naming the ' ...
                    'statistic (the median across participants of the posterior-mean ' ...
                    'coefficient minus its planted value)'], s.name);
            end
            if s.offsets
                if ~LS27WordingRules.hasCase(s.text, '(^|[^\d])-?\.13\d*(?!\d)')
                    v{end + 1} = sprintf(['%s: lacks the converged prospective fit''s D offset ' ...
                        '(-.13; median signed error -.1304)'], s.name);
                end
                if ~LS27WordingRules.hasCase(s.text, '(^|[^\d])-?\.02\d?(?!\d)')
                    v{end + 1} = sprintf('%s: lacks the retained location-scale control''s D offset (-.023)', s.name);
                end
            end
            for k = 1:numel(sen)
                if LS27WordingRules.dTerm(sen{k}) && ...
                        LS27WordingRules.hasCase(sen{k}, '(^|[^\d])\+?\.0(766|77|8)(?!\d)') && ...
                        ~LS27WordingRules.has(sen{k}, '(non-?|not |did not |without )converg')
                    v{end + 1} = sprintf(['%s: cites B''s D offset (+.0766, printed .08) without ' ...
                        'its non-convergence: "%s"'], s.name, LS27WordingRules.clip(sen{k})); %#ok<AGROW>
                end
            end
        end

        function v = ratioWording(s)
            %(3) Ratios are on the variance scale; the ranges belong to the
            %five scope captures; long surfaces give the retained control's
            %1.42, 1.66 and 1.02; no sentence makes 1.02 a pass.
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v); return; end
            t = s.text;
            if LS27WordingRules.has(t, LS27WordingRules.RATIO) && ...
                    ~LS27WordingRules.has(t, LS27WordingRules.VARSCALE)
                v{end + 1} = sprintf(['%s: prints component ratios without saying they are ' ...
                    'on the variance scale'], s.name);
            end
            sen = LS27WordingRules.sentences(t);
            hasRanges = false;
            for k = 1:numel(sen)
                if LS27WordingRules.hasCase(sen{k}, '12\.1(?!\d)')
                    hasRanges = true;
                    if ~LS27WordingRules.has(sen{k}, '\<five\>')
                        v{end + 1} = sprintf(['%s: the 2.4-12.1 / 2.2-7.2 / 1.1-1.9 ranges are not ' ...
                            'attributed to the five scope captures: "%s"'], s.name, ...
                            LS27WordingRules.clip(sen{k})); %#ok<AGROW>
                    end
                end
            end
            if s.long && hasRanges
                vals = {'1\.42', '1\.66', '1\.02'};
                for k = 1:numel(vals)
                    if ~LS27WordingRules.hasCase(t, ['(^|[^\d])' vals{k} '(?!\d)'])
                        v{end + 1} = sprintf('%s: lacks the retained control value %s', ...
                            s.name, strrep(vals{k}, '\', '')); %#ok<AGROW>
                    end
                end
            end
            for k = 1:numel(sen)
                if LS27WordingRules.hasCase(sen{k}, '(^|[^\d])1\.02(?!\d)')
                    v = [v, LS27WordingRules.affirmations(s.name, sen{k})]; %#ok<AGROW>
                end
            end
        end

        function v = intervalWording(s)
            %(4) Keep the interval-miss statement; where intervals are
            %advised, say that reading components against them does not
            %remedy the weak evidence.
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v); return; end
            if s.intervalMiss && ~LS27WordingRules.has(s.text, LS27WordingRules.MISS)
                v{end + 1} = sprintf('%s: dropped the statement that some intervals miss the truth', s.name);
            end
            if s.remedy
                sen = LS27WordingRules.sentences(s.text);
                ok = false;
                for k = 1:numel(sen)
                    if LS27WordingRules.has(sen{k}, '\<intervals?\>') && ...
                            LS27WordingRules.has(sen{k}, LS27WordingRules.REMEDY)
                        ok = true; break;
                    end
                end
                if ~ok
                    v{end + 1} = sprintf(['%s: does not say that reading the components against ' ...
                        'their intervals does not remedy the weak evidence'], s.name);
                end
            end
        end

        function v = consistencyWording(s)
            %(5) The 2e-12 agreement is implementation consistency, not
            %recovery, identifiability or sampling validation; no "not a
            %formula error" and no unhedged causal claim.
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v); return; end
            sen = LS27WordingRules.sentences(s.text);
            for k = 1:numel(sen)
                if LS27WordingRules.hasCase(sen{k}, '2e-12') && ...
                        ~(LS27WordingRules.has(sen{k}, LS27WordingRules.CONSISTENCY) && ...
                        LS27WordingRules.has(sen{k}, LS27WordingRules.NOTRECOVERY))
                    v{end + 1} = sprintf(['%s: states the 2e-12 agreement without calling it ' ...
                        'implementation consistency, not recovery: "%s"'], s.name, ...
                        LS27WordingRules.clip(sen{k})); %#ok<AGROW>
                end
                if LS27WordingRules.has(sen{k}, 'exponent') && ...
                        LS27WordingRules.has(sen{k}, '\<because\>') && ...
                        ~LS27WordingRules.has(sen{k}, LS27WordingRules.HEDGE)
                    v{end + 1} = sprintf('%s: states the mechanism as established cause: "%s"', ...
                        s.name, LS27WordingRules.clip(sen{k})); %#ok<AGROW>
                end
            end
            for k = 1:numel(LS27WordingRules.CERTAINTY)
                v = [v, LS27WordingRules.matches(s, LS27WordingRules.CERTAINTY{k}, ...
                    'certainty the 2e-12 agreement does not carry')]; %#ok<AGROW>
            end
        end

        function v = maximumWording(s)
            %No universal maximum ("up to .13" is not a bound).
            v = LS27WordingRules.anchorMissing(s);
            if ~isempty(v); return; end
            v = LS27WordingRules.matches(s, LS27WordingRules.MAXIMUM, 'states a universal maximum');
        end

        % ------------------------------------------------------------------
        % Q8: raw draws versus derived values
        % ------------------------------------------------------------------

        function v = q8Wording(name, units, needRaw, needCoef27, needSurface)
            %UNITS is a cell of strings, each checked on its own (record
            %strings, or the sentences of a help text or paragraph), so a
            %requirement cannot be met by words scattered across unrelated
            %strings. needRaw: some unit says the raw draws are unchanged.
            %needCoef27: some sentence says analysis 27's per-participant
            %coefficients change. needSurface: some sentence says a surface
            %can change.
            v = {};
            units = cellfun(@LS27WordingRules.normalize, units, 'UniformOutput', false);
            units = units(~cellfun(@isempty, units));
            if isempty(units)
                v{end + 1} = sprintf('%s: no text found', name);
                return;
            end
            sen = {};
            for k = 1:numel(units)
                m = regexpi(units{k}, LS27WordingRules.AMBIGUOUS, 'match');
                for j = 1:numel(m)
                    v{end + 1} = sprintf(['%s: "%s" does not say which values are unchanged ' ...
                        '(the raw draws; derived values can change): "%s"'], name, m{j}, ...
                        LS27WordingRules.clip(units{k})); %#ok<AGROW>
                end
                sen = [sen, LS27WordingRules.sentences(units{k})]; %#ok<AGROW>
            end
            if needRaw && ~any(cellfun(@(u) LS27WordingRules.has(u, LS27WordingRules.RAWUNCHANGED), sen))
                v{end + 1} = sprintf('%s: does not say that the raw draws are unchanged', name);
            end
            if needCoef27 && ~any(cellfun(@(u) LS27WordingRules.mentions(u, 27) && ...
                    LS27WordingRules.has(u, 'coefficient') && LS27WordingRules.has(u, '\<chang'), sen))
                v{end + 1} = sprintf(['%s: does not say that analysis 27''s per-participant ' ...
                    'coefficients change at the next summary'], name);
            end
            if needSurface && ~any(cellfun(@(u) LS27WordingRules.has(u, 'surface') && ...
                    LS27WordingRules.has(u, '\<(can|may|will|could)\s+(also\s+)?change'), sen))
                v{end + 1} = sprintf('%s: does not say that an explicit surface target can change a surface', name);
            end
        end

        function units = recordUnits(rec)
            %The user-visible strings of a migration record: method, notes,
            %limitations, invalidation and each stratum's transform.
            units = {};
            if isfield(rec, 'method'); units{end + 1} = rec.method; end
            for f = {'notes', 'limitations', 'invalidation'}
                if isfield(rec, f{1}) && iscell(rec.(f{1}))
                    units = [units, rec.(f{1})(:)']; %#ok<AGROW>
                end
            end
            if isfield(rec, 'transforms') && isstruct(rec.transforms) && isfield(rec.transforms, 'transform')
                units = [units, {rec.transforms.transform}];
            end
            units = units(cellfun(@(x) ischar(x) || isstring(x), units));
            units = cellfun(@char, units, 'UniformOutput', false);
        end

        function t = appendixDynamicParagraph(root)
            %Appendix C's paragraph on the dynamic designs (anchored on its
            %bold lead-in), or '' when the anchor is gone.
            ac = fileread(fullfile(root, 'documentation', 'manual', 'appendix-c_reprocessing.md'));
            t = LS27WordingRules.block(ac, '^\*\*Dynamic', {'27'});
        end

        % ------------------------------------------------------------------
        % Text utilities
        % ------------------------------------------------------------------

        function t = normalize(t)
            %One line of plain text: joined cells, no markdown emphasis or
            %code ticks, no help-browser tags, ASCII minus and spaces.
            if iscell(t); t = strjoin(cellfun(@char, t(:)', 'UniformOutput', false), ' '); end
            if isstring(t); t = strjoin(t(:)', ' '); end
            t = char(t);
            t = regexprep(t, '</?strong>', '');
            t = regexprep(t, '[`*]', '');
            t = strrep(t, char(8722), '-');
            t = strrep(t, char(160), ' ');
            t = strtrim(regexprep(t, '\s+', ' '));
        end

        function s = sentences(t)
            %Split after a sentence end (. ! ?) followed by spaces and a
            %capital, an opening parenthesis or an opening bracket. Decimals
            %(".03") and "e.g." do not split. Boundaries are located by match
            %indices, with no lookbehind and no escaped bracket inside a
            %character class: the executed RED run showed the earlier
            %lookbehind 'split' pattern returning every text unsplit, which
            %widened each sentence-scoped rule to the whole surface.
            [st, en] = regexp(t, '(\.|!|\?) +([A-Z]|\(|\[)', 'start', 'end');
            s = cell(1, numel(st) + 1);
            from = 1;
            for k = 1:numel(st)
                s{k} = strtrim(t(from:st(k)));
                from = en(k);   % the capital or bracket that opens the next sentence
            end
            s{end} = strtrim(t(from:end));
            s = s(~cellfun(@isempty, s));
        end

        function c = clauses(t)
            %Sentences, split again at semicolons.
            c = {};
            s = LS27WordingRules.sentences(t);
            for k = 1:numel(s)
                c = [c, regexp(s{k}, '\s*;\s*', 'split')]; %#ok<AGROW>
            end
        end

        function tf = has(t, pat)
            tf = ~isempty(regexpi(t, pat, 'once'));
        end

        function tf = hasCase(t, pat)
            tf = ~isempty(regexp(t, pat, 'once'));
        end

        function tf = mentions(t, n)
            %N as a standalone number, not part of a decimal or a longer number.
            tf = LS27WordingRules.hasCase(t, sprintf('(^|[^\\d.])%d(?!\\d)', n));
        end

        function tf = gTerm(t)
            tf = LS27WordingRules.hasCase(t, '\<G\>') || LS27WordingRules.has(t, 'generaliz');
        end

        function tf = dTerm(t)
            tf = LS27WordingRules.hasCase(t, '\<D\>') || LS27WordingRules.has(t, 'dependab');
        end

        function tf = biasTerm(t)
            tf = LS27WordingRules.has(t, 'bias|downward|\<low\>|\<lower\>|below|underestimat|negative') ...
                || LS27WordingRules.hasCase(t, '-\.\d');
        end

        function tf = states27Bias(c)
            %A clause that names analysis 27, both coefficients and a bias.
            tf = LS27WordingRules.mentions(c, 27) && LS27WordingRules.gTerm(c) && ...
                LS27WordingRules.dTerm(c) && LS27WordingRules.biasTerm(c);
        end

        function t = clip(t)
            if numel(t) > 160; t = [t(1:157) '...']; end
        end

    end

    methods (Static, Access = private)

        function v = anchorMissing(s)
            v = {};
            if ~s.found
                v{1} = sprintf(['%s: its content anchor was not found, so the rules could ' ...
                    'not run (update LS27WordingRules if the block moved)'], s.name);
            end
        end

        function v = matches(s, pat, why)
            %One violation per match of PAT in the surface text.
            v = {};
            [st, en] = regexpi(s.text, pat, 'start', 'end');
            for k = 1:numel(st)
                a = max(1, st(k) - 50); b = min(numel(s.text), en(k) + 50);
                v{end + 1} = sprintf('%s: %s: "...%s..."', s.name, why, s.text(a:b)); %#ok<AGROW>
            end
        end

        function v = affirmations(name, sentence)
            %An affirmative pass or recovery word with no negation in the 25
            %characters before it.
            v = {};
            st = regexpi(sentence, LS27WordingRules.AFFIRM, 'start');
            for k = 1:numel(st)
                pre = sentence(max(1, st(k) - 25):max(1, st(k) - 1));
                if ~LS27WordingRules.has(pre, LS27WordingRules.NEGATION)
                    v{end + 1} = sprintf('%s: the retained control''s E 1.02 reads as a pass: "%s"', ...
                        name, LS27WordingRules.clip(sentence)); %#ok<AGROW>
                end
            end
        end

        function t = block(txt, startPat, mustHave)
            %The first blank-line-delimited block whose first line matches
            %STARTPAT and that contains every pattern in MUSTHAVE.
            t = '';
            blocks = regexp(txt, '\r?\n[ \t]*\r?\n', 'split');
            for k = 1:numel(blocks)
                b = strtrim(blocks{k});
                if isempty(regexp(b, startPat, 'once')); continue; end
                if all(cellfun(@(p) ~isempty(regexp(b, p, 'once')), mustHave))
                    t = b;
                    return;
                end
            end
        end

        function t = paragraphsWith(txt, pat)
            %Every blank-line-delimited block that contains PAT (case-sensitive),
            %joined; '' when none does.
            blocks = regexp(txt, '\r?\n[ \t]*\r?\n', 'split');
            keep = blocks(~cellfun(@isempty, regexp(blocks, pat, 'once')));
            t = strjoin(cellfun(@strtrim, keep, 'UniformOutput', false), ' ');
        end

        function t = lineStarting(txt, prefix)
            %The first line that starts with PREFIX, with its indented
            %continuation lines.
            t = '';
            lines = regexp(txt, '\r?\n', 'split');
            i = find(startsWith(lines, prefix), 1);
            if isempty(i); return; end
            t = lines{i};
            j = i + 1;
            while j <= numel(lines) && ~isempty(regexp(lines{j}, '^[ \t]+\S', 'once'))
                t = [t ' ' strtrim(lines{j})]; %#ok<AGROW>
                j = j + 1;
            end
        end

        function t = manual54(m5)
            %Section 5.4's Restricted text: the table row and the paragraphs
            %that carry the capitalized label (the lower-case "restricted arm"
            %preamble is left out).
            t = '';
            sec = regexp(m5, '\n## 5\.4 .*?(?=\n## |$)', 'match', 'once');
            if isempty(sec); return; end
            keep = {};
            blocks = regexp(sec, '\r?\n[ \t]*\r?\n', 'split');
            for k = 1:numel(blocks)
                b = strtrim(blocks{k});
                if startsWith(b, '|')
                    rows = regexp(b, '\r?\n', 'split');
                    keep = [keep, rows(contains(rows, 'Restricted'))]; %#ok<AGROW>
                elseif contains(b, 'Restricted') && ~startsWith(b, '#')
                    keep{end + 1} = b; %#ok<AGROW>
                end
            end
            t = strjoin(keep, ' ');
        end

    end
end
