classdef TestSingViewerCoefficientControl < PsyRATTestBase
    % RC-32 regression: on every single-session viewer screen, the coefficient
    % control the user SEES must write the preference the summary route
    % actually CONSUMES.
    %
    % psyrat_startview_sing sets sserrvar for BOTH 'ic_sserrvar' and
    % 'ic_diff_sserrvar' but isdiff for 'ic_diff' only, so 'ic_diff_sserrvar'
    % used to match the sserrvar popup branch and get a control labeled
    % "G-Theory Coefficient" writing prefs.view.gcoeff. That route never reads
    % gcoeff: psyrat_relfigures passes prefs.view.diffgcoeff and
    % psyrat_relsummary_sing_diff_sserr takes only diffgcoeff, assigning
    % relsummary.gcoeff from it. The visible control was inert.
    %
    % Confirmed live on 2026-07-29 before the fix: with the popup set to
    % Generalizability, every rendered output still read "dependability" --
    % window titles ("Dependability Analyses", "Subject-Level Dependability
    % Estimates") and table headers ("Group SEM (dependability)") alike.
    %
    % Viewer code is usually called headlessly untestable, but `matlab -batch`
    % builds real figures, so the constructed controls can be read back. These
    % tests do that; they do not replace a click-through of the callbacks.

    methods (Test)

        function testCoefficientControlWritesThePreferenceTheRouteReads(testCase)
            % The load-bearing assertion, and deliberately NOT a label check:
            % gcoeff and diffgcoeff are given DIFFERENT values, so the popup's
            % initial Value identifies which preference it is bound to, whatever
            % it is called on screen.
            %
            %   gcoeff = 1 (Dependability), diffgcoeff = 2 (Generalizability)
            %
            % so a popup showing Value 2 is bound to diffgcoeff and one showing
            % Value 1 is bound to gcoeff.
            cases = {
                'ic',                1, 'population one-facet: relsummary takes gcoeff'
                'ic_sserrvar',       1, 'subject-level one-facet: relsummary takes gcoeff'
                'ic_diff',           2, 'difference: relsummary.gcoeff comes FROM diffgcoeff'
                'ic_diff_sserrvar',  2, 'subject-level difference: likewise (RC-32)'
                };

            for k = 1:size(cases, 1)
                analysis = cases{k, 1};
                expected = cases{k, 2};
                why      = cases{k, 3};

                value = localOpenViewerPopup(testCase, analysis);
                testCase.verifyEqual(value, expected, sprintf( ...
                    ['%s: coefficient popup should initialize from the ' ...
                    'preference its summary route reads (%s). Value %d ' ...
                    'means it is bound to the other one.'], ...
                    analysis, why, value));
            end
        end

        function testDifferenceScreensAreLabelledAsDifferenceScreens(testCase)
            % The user-facing half: a difference design must not present the
            % non-difference control's label.
            expectLabel = {
                'ic',               'G-Theory Coefficient:'
                'ic_sserrvar',      'G-Theory Coefficient:'
                'ic_diff',          'Difference-Score Coefficient:'
                'ic_diff_sserrvar', 'Difference-Score Coefficient:'
                };

            for k = 1:size(expectLabel, 1)
                analysis = expectLabel{k, 1};
                want     = expectLabel{k, 2};

                [~, label] = localOpenViewerPopup(testCase, analysis);
                testCase.verifyEqual(label, want, sprintf( ...
                    '%s should show "%s", got "%s".', analysis, want, label));
            end
        end

        function testOnlyPlainOneFacetOffersCriterionScoring(testCase)
            % RC-05, guarded BEHAVIORALLY rather than by matching source text.
            %
            % An earlier version of this test asserted that two literal strings
            % were still present in psyrat_startview_sing.m, under the name
            % ...WideningIsdiffWouldBreakTheCriterionGate. Both the name and the
            % method were wrong. The gate is a conjunction of negations
            % (showCriterionButton = ~sserrvar && ~isdiff && hasEstimatedMu), so
            % WIDENING isdiff could only ever suppress the button further, never
            % expose it -- the hazard RC-05 actually names is REMOVING a term.
            % And a source-text assertion passes or fails on whitespace, which
            % is not the invariant anyone cares about.
            %
            % What matters is the outcome: criterion scoring is offered on the
            % plain one-facet screen and on no other, because only there is
            % mu estimated rather than substituted with mean(bp).
            expectButton = {
                'ic',               true
                'ic_sserrvar',      false
                'ic_diff',          false
                'ic_diff_sserrvar', false
                };

            for k = 1:size(expectButton, 1)
                analysis = expectButton{k, 1};
                want     = expectButton{k, 2};

                has = localScreenHasCriterionButton(testCase, analysis);
                testCase.verifyEqual(has, want, sprintf( ...
                    ['%s: criterion-scoring button present=%d, expected %d. ' ...
                    'Only the plain one-facet screen may offer it (RC-05).'], ...
                    analysis, has, want));
            end
        end

    end
end


function has = localScreenHasCriterionButton(testCase, analysis)
%True when the viewer for this analysis builds a criterion-scoring button.
%Read inside the helper for the same reason the popup helper returns values:
%the onCleanup below deletes the handles when this returns.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force'));

psyrat_startview('psyrat_prefs', localViewPrefs(), ...
    'psyrat_data', localViewData(analysis));

f = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
testCase.assertNotEmpty(f, sprintf('%s opened no viewer figure.', analysis));

has = false;
btns = findall(f(1), 'Style', 'pushbutton');
for i = 1:numel(btns)
    s = btns(i).String;
    if iscell(s); s = strjoin(s, ' '); end
    if contains(string(s), 'Criterion', 'IgnoreCase', true)
        has = true;
    end
end
end


function [value, label] = localOpenViewerPopup(testCase, analysis)
%Open the single-session viewer for one analysis label and return its
%coefficient popup's VALUE plus the static text immediately to its left.
%
%Values, not handles: the onCleanup below closes the figure when this helper
%returns, so a returned handle would already be deleted at the call site.
close('all', 'force');
cleanup = onCleanup(@() close('all', 'force'));

psyrat_data = localViewData(analysis);

prefs = localViewPrefs();

psyrat_startview('psyrat_prefs', prefs, 'psyrat_data', psyrat_data);

f = findall(groot, 'Type', 'figure', 'Tag', 'psyrat_gui');
testCase.assertNotEmpty(f, sprintf('%s opened no viewer figure.', analysis));
f = f(1);

pops = findall(f, 'Style', 'popupmenu');
testCase.assertNumElements(pops, 1, sprintf( ...
    '%s should build exactly one coefficient popup, found %d.', ...
    analysis, numel(pops)));
pop = pops(1);
value = pop.Value;

%the label is the static text on the same row, nearest to the popup's left
label = '';
bestdx = inf;
txts = findall(f, 'Style', 'text');
pp = pop.Position;
for j = 1:numel(txts)
    tp = txts(j).Position;
    if abs(tp(2) - pp(2)) < 14 && tp(1) < pp(1)
        dx = pp(1) - tp(1);
        if dx < bestdx
            bestdx = dx;
            label = txts(j).String;
        end
    end
end
if iscell(label); label = strjoin(label, ' '); end
label = strtrim(char(label));
end


function prefs = localViewPrefs()
%Start from the shipped defaults rather than a hand-rolled struct, so this test
%does not have to track every field the viewer reads (it already broke once on
%a missing view.plotdep).
%
%gcoeff and diffgcoeff are then set to DIFFER on purpose -- that is what makes
%the popup's binding observable from its initial Value alone.
prefs = psyrat_defaults;
prefs.view.gcoeff = 1;       % Dependability
prefs.view.diffgcoeff = 2;   % Generalizability
end


function psyrat_data = localViewData(analysis)
%Summarize a real fixture through the production relsummary path, then stamp
%the analysis label the viewer routes on.
switch analysis
    case 'ic'
        pd0 = PsyRATTestDataFactory.makeICSummaryData();
        psyrat_data = pd0;
    case 'ic_sserrvar'
        pd0 = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();
        psyrat_data = psyrat_relsummary('psyrat_data', pd0, ...
            'analysis', 'sing_sserr', 'depcutoff', 0.01, 'meascutoff', 2, ...
            'depcentmeas', 1, 'gcoeff', 1, 'CI', 0.95);
    case 'ic_diff'
        pd0 = PsyRATTestDataFactory.makeDoDRelData();
        vd = psyrat_dod_build_virtual_rel('psyrat_data', pd0, ...
            'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);
        psyrat_data = psyrat_relsummary('psyrat_data', vd, ...
            'analysis', 'sing_diff', 'depcutoff', 0.5, 'meascutoff', 2, ...
            'depcentmeas', 1, 'diffgcoeff', 1, 'CI', 0.95);
    case 'ic_diff_sserrvar'
        pd0 = PsyRATTestDataFactory.makeICDiffSSErrRelData();
        psyrat_data = psyrat_relsummary('psyrat_data', pd0, ...
            'analysis', 'sing_diff_sserr', 'depcutoff', 0.5, ...
            'meascutoff', 2, 'depcentmeas', 1, 'diffgcoeff', 1, 'CI', 0.95);
end

psyrat_data.rel.analysis = analysis;
if ~isfield(psyrat_data, 'ver'); psyrat_data.ver = '0-test'; end
if ~isfield(psyrat_data.rel, 'filename'); psyrat_data.rel.filename = 'test'; end
%the viewer prints the measurement header; the fixtures carry no proc block
if ~isfield(psyrat_data, 'proc') || ~isfield(psyrat_data.proc, 'measheader')
    psyrat_data.proc.measheader = 'meas';
end
end
