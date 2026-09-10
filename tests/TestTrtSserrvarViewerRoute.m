classdef TestTrtSserrvarViewerRoute < PsyRATTestBase
    % RC-24 regression: the subject-level test-retest design (trt_sserrvar,
    % analysis 25) must survive the GUI view route.
    %
    % psyrat_startview_trt hands psyrat_relfigures the label 'trt' for every
    % test-retest file except trt_diff, and relfigures re-detects the
    % subject-level design from rel.analysis. That contract was already in place
    % for the FIGURES (the subject-level plot gate, the psyrat_variancet table,
    % and the duplicate-suppression check) but not for the relsummary CALL: the
    % 'trt' branch dereferences REL.out.sig_id, a field analysis 25 never writes,
    % so View Results died with "Unrecognized field name sig_id".
    %
    % These tests are the headless half. The on-screen half needs a click-through.

    methods (Test)

        function testRelsummaryTrtBranchStillNeedsSigId(testCase)
            % Pin the ROOT CAUSE, so the regression is anchored to why the fix
            % is needed and not just to its symptom. Sending a trt_sserrvar REL
            % through the plain 'trt' branch must still fail -- the fix routes
            % around this branch, it does not change it.
            pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
            testCase.verifyError(@() psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt', 'gcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.7, 'meascutoff', 2, 'relcentmeas', 1, ...
                'noccmode', 1), 'MATLAB:nonExistentField');
        end

        function testRelfiguresRoutesTrtSserrvarWithoutError(testCase)
            % The RC-24 crash itself: the GUI entry point with the label the trt
            % viewer actually sends. Every output pref is ON, so if any of the
            % five unavailable group-level outputs were still attempted, this
            % would throw (each was verified to error, not to return a wrong
            % number). Producing figures without error is the whole assertion.
            pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
            prefs = localAllOnPrefs();
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            close('all','force');
            nBefore = numel(findall(groot,'Type','figure'));
            psyrat_relfigures('psyrat_data', pd, 'psyrat_prefs', prefs, ...
                'analysis', 'trt');
            nAfter = numel(findall(groot,'Type','figure'));

            testCase.verifyGreaterThan(nAfter, nBefore, ...
                'psyrat_relfigures(trt) on a trt_sserrvar file produced nothing.');
        end

        function testSubjectLevelPlotsCarryTheTwoFacetReferenceLine(testCase)
            % The route must reach the subject-level caterpillar plots, and the
            % reference line drawn there must be the stored two-facet value
            % (RC-06's gro_rel/gro_icc), not the one-facet fallback. This is the
            % headless stand-in for the click-through RC-24 blocked.
            pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
            [summ, ~] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', 3, ...
                'depcutoff', 0.7, 'meascutoff', 2, 'depcentmeas', 1, ...
                'noccmode', 1);
            T = summ.relsummary.group(1).event(1).ssrel_table;
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            spec = {'dep', 'gro_rel'; 'icc', 'gro_icc'};
            for k = 1:size(spec,1)
                close('all','force');
                psyrat_ssrelplot('psyrat_data', summ, 'stat', spec{k,1});
                h = findobj(gcf, 'Type', 'line', 'LineWidth', 1.5);
                testCase.verifyNotEmpty(h, sprintf( ...
                    'no reference line on the %s panel', spec{k,1}));
                y = get(h(1), 'YData');
                testCase.verifyEqual(y(1), T.(spec{k,2})(1), 'AbsTol', 1e-12, ...
                    sprintf('%s reference line is not the stored %s', ...
                    spec{k,1}, spec{k,2}));
            end
        end

        function testSummaryRecordsWhichCoefficientProducedTheNumbers(testCase)
            % RC-26. The 'trt_sserr' branch computes its per-participant
            % coefficient FROM reltype but did not store it, so the exported
            % header could not say whether a number was CE, CS or CES -- the gap
            % RC-12's provenance line exists to close, left open on this design
            % because that line gates on the presence of relsummary.reltype.
            %
            % Asserted at BOTH ends: the field on the summary, and the header
            % line the CSV exports actually write (psyrat_savesscoeffs and the
            % variance-table saves all route through psyrat_fprintf_provenance,
            % which is psyrat_provenance_lines). Two reltype values are run, so
            % the line has to track the selection rather than being a constant.
            spec = { ...
                2, 'trt',    'Coefficient of Stability (CS)'; ...
                3, 'ic_trt', 'Coefficient of Equivalence and Stability (CES)'};

            for k = 1:size(spec,1)
                reltype = spec{k,1};
                pd = PsyRATTestDataFactory.makeTRTSSErrRelData();
                [summ, ~] = psyrat_relsummary('psyrat_data', pd, ...
                    'analysis', 'trt_sserr', 'gcoeff', 1, 'reltype', reltype, ...
                    'depcutoff', 0.7, 'meascutoff', 2, 'depcentmeas', 1, ...
                    'noccmode', 1);

                testCase.assertTrue(isfield(summ.relsummary,'reltype'), ...
                    ['the trt_sserr summary must record reltype -- without it ' ...
                    'psyrat_trt_coeflabel returns '''' and every export is blind ' ...
                    'to which of CE/CS/CES produced its numbers']);
                testCase.verifyEqual(summ.relsummary.reltype, reltype, ...
                    'the trt_sserr summary must record the reltype it computed from');
                testCase.verifyEqual(summ.relsummary.reltype_name, spec{k,2}, ...
                    'reltype_name must match the sibling branches'' mapping');

                prov = psyrat_provenance_lines(summ);
                hit = prov(startsWith(prov, 'Coefficient:'));
                testCase.verifyNumElements(hit, 1, ...
                    'exported headers must carry exactly one Coefficient: line');
                testCase.verifySubstring(hit{1}, spec{k,3}, ...
                    'the Coefficient: line must name the selected estimand');
                % The error type and n'_o clauses ride along on the same line;
                % without them the label is not a complete estimand.
                testCase.verifySubstring(hit{1}, 'Dependability (absolute error)');
                testCase.verifySubstring(hit{1}, 'occasion');
            end
        end

        function testViewerDisablesUnavailableOutputsForTrtSserrvar(testCase)
            % The five group-level outputs cannot be produced for this design, so
            % the viewer must not offer them. Plain 'trt' must be unaffected.
            prefs = localStartviewPrefs();
            rel = struct('analysis','trt','filename','synthetic', ...
                'time',{{'o1','o2'}},'groups','none','events','none', ...
                'nchains',2,'niter',1);
            pd = struct('rel',rel,'proc',struct('measheader','meas'),'ver','0-test');
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            close('all','force');
            psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
            testCase.verifyEqual(localNumDisabledCheckboxes(), 0, ...
                'plain trt should not disable any output checkbox');

            pd.rel.analysis = 'trt_sserrvar';
            close('all','force');
            psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
            testCase.verifyEqual(localNumDisabledCheckboxes(), 5, ...
                ['trt_sserrvar should disable exactly the five group-level ' ...
                'outputs (trials curve, ICC plot, between-person SD plot, ' ...
                'cutoff table, overall table)']);
        end

        function testViewerOffersCriterionButtonOnlyWhereItCanRun(testCase)
            % G38: the Criterion Score Outputs button was offered
            % unconditionally, but the criterion route needs the group-level
            % components psyrat_cutscore_trt consumes. trt_sserrvar estimates a
            % per-participant residual instead of them, and the trt_diff
            % D-study stores no mu, so on both designs the button was a live
            % crash path (MATLAB:nonExistentField out of a uicontrol callback,
            % the same shape RC-24 fixed on the figures route). The viewer must
            % suppress the button for those two and keep it for plain trt --
            % the gate psyrat_startview_sing already applies on its side.
            prefs = localStartviewPrefs();
            rel = struct('analysis','trt','filename','synthetic', ...
                'time',{{'o1','o2'}},'groups','none','events','none', ...
                'nchains',2,'niter',1);
            pd = struct('rel',rel,'proc',struct('measheader','meas'),'ver','0-test');
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            close('all','force');
            psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
            testCase.verifyTrue(localHasCriterionButton(), ...
                'plain trt must keep offering Criterion Score Outputs');
            nWithButton = localNumPushButtons();

            for a = {'trt_sserrvar','trt_diff'}
                pd.rel.analysis = a{1};
                close('all','force');
                psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);
                testCase.verifyFalse(localHasCriterionButton(), sprintf( ...
                    '%s must not offer the criterion route it cannot run', a{1}));
                testCase.verifyEqual(localNumPushButtons(), nWithButton - 1, ...
                    sprintf(['%s: exactly the criterion button must be gone; ' ...
                    'the rest of the row stays'], a{1}));
            end
        end

        function testViewerControlsAllFitInsideTheWindow(testCase)
            %RC-25. This screen accumulates fixed-pixel row offsets that ran off
            %the bottom of its own figure. The bottom button row sat at
            %y = -23.75 on 'trt' and 'trt_sserrvar' -- 26 of its 50 px inside --
            %and at y = -68.75 on 'trt_diff', where NONE of it was inside and the
            %user could not see or press Back to Home, Open Another, Generate,
            %Output Prefs, Criterion Score or Close PsyRAT at all. The entry
            %filed it as cosmetic and font-size dependent; it is neither. Every
            %Position on this screen is an expression in
            %figheight/rowspace/lcol/rcol, and fsize only sets 'fontsize'.
            %
            %'trt_diff' is the case that matters and is easy to miss: it adds one
            %more row than the other two, so a check that only covers 'trt'
            %passes while the worst case stays broken. All three labels are
            %driven here for that reason.
            %
            %The window-fits assertion is the other half of the ruling. The
            %origin is hardcoded, so raising figheight without lowering the
            %origin would move the window's TOP off a 1728x1117 display -- which
            %the pre-fix 850-tall window at y = 310 already did (top 1160). 1117
            %is asserted rather than the 1440 of the current machine precisely so
            %a future height increase cannot silently reintroduce that.
            prefs = localStartviewPrefs();
            rel = struct('analysis','trt','filename','synthetic', ...
                'time',{{'o1','o2'}},'groups','none','events','none', ...
                'nchains',2,'niter',1);
            pd = struct('rel',rel,'proc',struct('measheader','meas'),'ver','0-test');
            cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

            for a = {'trt','trt_sserrvar','trt_diff'}
                pd.rel.analysis = a{1};
                close('all','force');
                psyrat_startview_trt('psyrat_prefs',prefs,'psyrat_data',pd);

                fig = findobj('Tag','psyrat_gui');
                testCase.assertNotEmpty(fig, ...
                    sprintf('%s did not build a viewer figure', a{1}));
                fig = fig(1);
                u = findall(fig,'Type','uicontrol');
                testCase.assertNotEmpty(u, ...
                    sprintf('%s viewer figure holds no uicontrols', a{1}));

                ys = arrayfun(@(h) h.Position(2), u);
                hs = arrayfun(@(h) h.Position(4), u);

                testCase.verifyGreaterThanOrEqual(min(ys), 0, sprintf( ...
                    ['%s: every control must start at or above the figure ' ...
                    'bottom edge (lowest y = %.2f)'], a{1}, min(ys)));
                testCase.verifyLessThanOrEqual(max(ys + hs), fig.Position(4), ...
                    sprintf(['%s: every control must end at or below the figure ' ...
                    'top edge (highest y+h = %.2f, figheight = %.2f)'], ...
                    a{1}, max(ys + hs), fig.Position(4)));
                testCase.verifyLessThanOrEqual( ...
                    fig.Position(2) + fig.Position(4), 1117, sprintf( ...
                    ['%s: the whole window must fit a 1728x1117 display ' ...
                    '(origin %.0f + height %.0f)'], ...
                    a{1}, fig.Position(2), fig.Position(4)));
            end
        end

    end
end

function n = localNumDisabledCheckboxes()
%Count checkbox uicontrols currently disabled across open figures.
h = findall(groot, 'Style', 'checkbox', 'Enable', 'off');
n = numel(h);
end

function tf = localHasCriterionButton()
%True when any open pushbutton carries the Criterion Score Outputs label (a
%2-row cellstr on screen, so match on its first row).
h = findall(groot, 'Style', 'pushbutton');
tf = false;
for i = 1:numel(h)
    s = cellstr(h(i).String);
    if any(contains(s, 'Criterion Score'))
        tf = true;
        return;
    end
end
end

function n = localNumPushButtons()
%Count pushbutton uicontrols across open figures.
n = numel(findall(groot, 'Style', 'pushbutton'));
end

function v = localViewPrefsStruct()
%Every output turned ON, so a missing gate surfaces as an error.
v = struct();
v.gcoeff = 1; v.reltype = 3; v.diffgcoeff = 1;
v.noccmode = 1; v.nocc = [];
v.relvalue = 0.5; v.relcentmeas = 1; v.meascutoff = 2;
v.plotrel = 1; v.plotrelline = 2;
v.ploticc = 1; v.inctrltable = 1; v.overalltable = 1;
v.ntrials = 20; v.showstddevt = 1; v.showstddevf = 1;
v.plotssrel = 1; v.plotssicc = 1; v.tablessrel = 1;
end

function prefs = localAllOnPrefs()
prefs = struct('view', localViewPrefsStruct(), 'ver', '0-test');
end

function prefs = localStartviewPrefs()
prefs = struct('guis', struct('fsize', 10), ...
    'view', localViewPrefsStruct(), 'ver', '0-test');
end
