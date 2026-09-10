classdef TestModalRejectionDialogs < PsyRATTestBase
    %Source pins for B21 stage 1 -- the SEVEN input-rejection errordlg sites
    %must be modal and carry DISTINCT, purpose-named titles.
    %
    %WHY. A single-argument errordlg is non-modal, untitled ('Error Dialog'),
    %and opens at a fixed position that the parent figure covers at default
    %placement: live click-throughs (2026-08-14 ntrls; re-verified 2026-08-16)
    %showed the user sees NOTHING on a rejected entry, and repeated rejections
    %stack identical dialogs at one spot so working OK clicks look ignored.
    %MATLAB's 'modal' CreateMode both raises the dialog and REPLACES an
    %existing message box with the same Title (verified empirically on R2025b:
    %same-title modal replaces; it also collapses previously stacked
    %non-modal same-title boxes), which dissolves the stacking half of B22.
    %Titles are per PURPOSE: same-purpose replacement is wanted (a stale
    %rejection should be superseded), while distinct purposes must not
    %clobber each other's pending messages -- hence three titles across the
    %seven sites, not one.
    %
    %The four post-computation NOTIFICATION sites in psyrat_criterionfigures /
    %psyrat_relfigures are deliberately NOT in this rollout (separate ruling;
    %see B21 in PSYRAT_AUDIT_FINDINGS.md).
    %
    %WAVE 2 (2026-08-17): the titled two-argument errordlg sites classified as
    %INPUT/CONFIG REJECTIONS on interactive paths -- 67 sites across nine
    %files -- are modal with their existing purpose titles. The twelve
    %interactive operation-FAILURE reporters (catch-ME sites in
    %psyrat_startimportscore / psyrat_startview_dodiff / psyrat_formula_guide,
    %plus startproc's end-of-batch 'All Measurements Failed') joined them
    %under their own ruling later the same day: modal with their existing
    %titles, and the end-of-batch site additionally REORDERED because it
    %fired before its own psyrat_startproc_gui relaunch. The notification /
    %advisory class then took its own ruling the same day: the four D-study
    %notices KEEP raised-last non-modal and gain purpose titles, the
    %relsummary threshold notice stays as-is (reachable only on nogooddata
    %routes, where every caller returns before drawing output; the modal
    %question there is its own ruling), and the iteration advisory MOVES to
    %fire after the Preferences screen builds. All still non-modal, pinned
    %below so an unruled modal conversion trips a test. Eight
    %relaunch-then-dialog sites (psyrat_startview, the six
    %installer/updater cancels, and the end-of-batch report) are REORDERED to
    %create the dialog LAST -- the B20 creation-order lesson: modality binds
    %only the windows that exist when the dialog is created.
    %
    %WHY SOURCE SCANS. The behavioral halves are covered elsewhere: the
    %reloverallt dialog's census in TestNonGuiComplexLiteralValidation keys on
    %the new title (and asserts modality), and the viewer flows are exercised
    %by TestViewerComplexLiteralValidation through exclusion logic that is
    %title-agnostic. These pins hold the seven call sites' argument form.
    %
    %B44 (2026-08-22) ADDS A DIFFERENT KIND OF PIN. Everything above is about
    %the CreateMode argument of an errordlg/warndlg. psyrat_start's TWO
    %update prompts are not dialog functions at all -- they are plain figures
    %-- so their modality is the figure's WindowStyle property. They are the
    %toolbox's only figure-level modals, are pinned separately below, and are
    %NOT part of the eleven CreateMode censuses. A third prompt in the same
    %file, psyrat_ask2install, stays non-modal on purpose: psyrat_start
    %returns immediately after raising it, without building the home screen,
    %so a modal there would have nothing to bind.
    %
    %B45 (2026-08-22) ADDS THE FIRST MESSAGE-BODY WORDING PIN. Not modality
    %at all. Dialog TITLES are already pinned as exact strings above; what is
    %new is pinning a sentence inside a message. The two figure-level prompts
    %share that sentence and their Tag, so it did not say which prompt it
    %belonged to -- and the destructive one carried the sibling's noun. They
    %are NOT otherwise identical (different opening sentence, only one has the
    %WARNING clause, text controls 400x40 vs 500x75), so the pin is per
    %FUNCTION, not a file-level count. It is NOT part of the eleven
    %CreateMode censuses and does not use ModalArg.
    %
    %B35 (2026-08-21) -- WHAT THESE PINS SCAN, AND WHY IT IS NOT THE RAW FILE.
    %Every count below used to be a plain substring search over the raw file
    %text, which moved for reasons unrelated to any dialog's modality, in
    %both directions and both silently:
    %  (a) UNDER-COUNT, the dangerous one. A site whose modal argument sits
    %      after a line continuation contains no ",'modal')" substring. The
    %      first draft of the B33 dialog was written that way: this census
    %      read 38 while a 39th modal dialog existed. A dialog could be
    %      turned NON-MODAL by reformatting alone and the pin would not
    %      notice -- which is the exact defect class the pin exists to catch.
    %  (b) OVER-COUNT. A COMMENT quoting the needle was counted. The comment
    %      written to document this very pin pushed the count to 40.
    %The fix is two halves, and BOTH are needed -- neither closes B35 alone:
    %  * readCode (not readSource) strips comments and joins continuations,
    %    which answers (b) and makes the source text safe to quote in prose.
    %  * the counts that must not silently drift are REGEXES tolerant of the
    %    whitespace a joined continuation leaves behind, which answers (a).
    %    A substring over normalized text does NOT do this: joining a
    %    continuation leaves the comma and the quote separated by blanks,
    %    which the fixed needle ",'modal')" still fails to match.
    %DIVISION OF LABOUR. Exactly FIVE pins still read readSource, and they
    %are the relaunch-order pins plus the '%recommend' placement. Two of
    %them REQUIRE raw text -- one matches a COMMENT, one a literal
    %continuation -- and read zero against normalized code; the other three
    %assert statement ORDER, where normalizing would additionally let an
    %intervening comment satisfy them. All five assert ==1 or ==3, so a
    %reformat there fails LOUDLY rather than silently, and that residual is
    %accepted rather than papered over. Everything else reads readCode.
    %Do not move a pin to readSource for convenience: a raw pin can be
    %satisfied by a COMMENT that merely describes the code it claims to
    %hold, which is how the uiwait pin was written until review caught it.

    properties (Constant, Access = private)
        %the modal CreateMode argument as it appears at a dialog call site,
        %tolerant of the run of blanks that normalizeSource leaves wherever
        %a line continuation was joined away. Kept in one place because
        %eleven file-level censuses, across nine call sites, share it.
        ModalArg = ',[ \t]*''modal''[ \t]*\)'
    end

    methods (Test)

        function testRecalcTrialCountSiteIsModalTitled(testCase)
            code = testCase.readCode(fullfile('subroutines', 'tables', ...
                'psyrat_trt_reloverallt.m'));
            testCase.verifyEqual( ...
                testCase.countOf(code, '''Invalid trial count'',''modal'');'), 1, ...
                ['the recalculation trial-count rejection must be modal with ' ...
                'its purpose title (B21 stage 1)']);
        end

        function testTrtViewerSitesAreModalTitled(testCase)
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startview_trt.m'));
            testCase.verifyEqual( ...
                testCase.countOf(code, 'errordlg(errorstr,''Invalid reliability cutoff'',''modal'');'), 2, ...
                'both trt cutoff rejections must be modal with the cutoff title');
            testCase.verifyEqual( ...
                testCase.countOf(code, '''Invalid plot trial count'',''modal'');'), 1, ...
                'the trt plot trial-count rejection must be modal with its title');
            %tolerant regex, not a substring: this is a ==0 pin, so a bare
            %errordlg reintroduced across a continuation would EVADE a
            %substring search and the suite would stay green (B35, under-count)
            testCase.verifyEqual(testCase.countRe(code, ...
                'errordlg[ \t]*\([ \t]*errorstr[ \t]*\)'), 0, ...
                'no bare single-argument errordlg(errorstr) may remain in the trt viewer');
        end

        function testSingViewerSitesAreModalTitled(testCase)
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startview_sing.m'));
            testCase.verifyEqual( ...
                testCase.countOf(code, 'errordlg(errorstr,''Invalid reliability cutoff'',''modal'');'), 2, ...
                'both sing cutoff rejections must be modal with the cutoff title');
            testCase.verifyEqual( ...
                testCase.countOf(code, '''Invalid plot trial count'',''modal'');'), 1, ...
                'the sing plot trial-count rejection must be modal with its title');
            testCase.verifyEqual(testCase.countRe(code, ...
                'errordlg[ \t]*\([ \t]*errorstr[ \t]*\)'), 0, ...
                'no bare single-argument errordlg(errorstr) may remain in the sing viewer');
        end

        function testWaveTwoStartprocRejectionsAreModal(testCase)
            %src is the RAW file, read by the two ORDERING pins that
            %deliberately match a comment and a literal continuation; code
            %is comments-stripped and continuations joined, read by
            %everything else (B35)
            src = testCase.readSource(fullfile('subroutines', 'guis', ...
                'psyrat_startproc.m'));
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startproc.m'));
            %36 rejection errordlg + the end-of-batch failure report + the
            %B20 modal warndlg + the B28 DoD single-measurement advisory + the
            %G8 preflight-warnings advisory (pre-beta batch, 2026-09-05: the
            %few-participants / low-median-trials / incomplete-cell warnings
            %the preflight computes, raised modally before the fit so a long
            %run does not start on a design the user never saw flagged; the
            %same precondition class as B28) = 40.
            %The 36th rejection is the B33 GUI follow-up: difference-of-
            %differences combined with residual covariance, rejected inside
            %psyrat_exec BEFORE the uiputfile save dialog so the engine's
            %varargin:familyconcurrent (analysis 9) and varargin:dynrel
            %(analyses 28/29) stop surfacing through the per-measurement catch
            %only after the user has already chosen an output location.
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 40, ...
                ['the 36 startproc input-rejection dialogs plus the ' ...
                'end-of-batch failure report must be modal (37 errordlg + ' ...
                'the B20 warndlg + the B28 advisory + the G8 preflight ' ...
                'advisory = 40); a different count ' ...
                'means a site was missed, dropped, or an allowlisted site ' ...
                'was converted without its ruling']);
            %B28 (ruled 2026-08-18, after the live click-through CONFIRMED the
            %defect): the DoD single-measurement advisory is uiwait + MODAL. It
            %reports a silent NARROWING OF ANALYSIS SCOPE -- the user selected N
            %measurements and gets 1 -- which is the operation/precondition
            %class ruling 1 already made modal, not the post-output
            %notification class ruling 2 kept non-modal. Observed defect: as a
            %non-modal warndlg fired BEFORE the setup GUI closes, it survived
            %the close, survived the psyrat_startproc_gui relaunch buried
            %underneath it, and ESCAPED the end-of-batch modal's binding
            %(created first, so modality never bound it) -- a handle census
            %found it still open with vis=on after the modal was dismissed.
            %
            %uiwait IS PINNED SEPARATELY BECAUSE IT IS THE HALF THAT FIXES B28.
            %An earlier draft of this ruling used 'modal' alone and claimed it
            %removed the orphaning "by construction"; that was refuted before
            %commit. warndlg returns immediately whatever its CreateMode -- only
            %uiwait blocks -- and psyrat_startproc has no uiwait/drawnow/pause
            %elsewhere, so without it execution runs through close(psyrat_gui)
            %into the measurement loop with no yield to the event queue, the
            %user never gets to dismiss the box, and it can still be open at the
            %relaunch. A pin on the modal argument alone would PASS on that
            %broken form, so it cannot be the only assertion here.
            testCase.verifyEqual(testCase.countOf(code, ...
                '''DoD mode: single measurement only'',''modal''));'), 1, ...
                ['the DoD single-measurement advisory must be modal: it ' ...
                'discloses that the requested analysis scope was narrowed ' ...
                'to one measurement, and non-modal it was provably ' ...
                'unseeable on the all-fail path (B28)']);
            %reads CODE, not raw: on raw text this pin is satisfied by a
            %COMMENT describing the wrapper, so the real uiwait( could be
            %deleted and the pin would stay green. It is not an ordering
            %pin and reads the same value either way.
            testCase.verifyEqual(testCase.countRe(code, ...
                'uiwait\(warndlg\(\{''Difference-of-differences mode'), 1, ...
                ['the advisory must be wrapped in uiwait: modality alone ' ...
                'does not block, so without it the run starts before the ' ...
                'user can acknowledge that their analysis scope was ' ...
                'narrowed, and the dialog can still reach the relaunch']);
            %the iteration advisory (2026-08-17 notification ruling): stays
            %non-modal and MOVED to fire after the Preferences screen is
            %fully built, so the screen no longer buries it by creation order
            testCase.verifyEqual( ...
                testCase.countOf(code, '''Recommended Iteration Increase'');'), 1, ...
                'the iteration advisory must stay non-modal (notification ruling)');
            %RAW on purpose: this pattern matches a COMMENT, so it reads
            %zero against comment-stripped text
            testCase.verifyEqual(testCase.countRe(src, ...
                'Tag = ''psyrat_gui'';\s+%recommend'), 1, ...
                ['the advisory must fire AFTER the Preferences screen ' ...
                'builds (it previously fired before the figure existed, ' ...
                'so the screen buried it by creation order)']);
            %the failure-reporter ruling (2026-08-17): modal, and reordered
            %to relaunch FIRST so the new screen neither covers nor escapes
            %the dialog (B20 creation order)
            testCase.verifyEqual( ...
                testCase.countOf(code, '''All Measurements Failed'',''modal'');'), 1, ...
                'the end-of-batch failure report must be modal with its title');
            %RAW on purpose: this pattern matches a literal line
            %continuation, which normalizeSource joins away
            testCase.verifyEqual(testCase.countRe(src, ...
                ['psyrat_startproc_gui\([^\n]*\);\s+errordlg\(\{\.\.\.' ...
                '\s+sprintf\(''All %d selected']), 1, ...
                ['the end-of-batch failure report must relaunch FIRST and ' ...
                'create the modal dialog LAST (B20 creation order)']);
        end

        function testWaveTwoImportScoreGuardsAreModal(testCase)
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startimportscore.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 14, ...
                ['all 14 import/score dialogs must be modal: the eight ' ...
                'precondition guards (wave 2) plus the six catch-ME ' ...
                'operation-failure reporters (the 2026-08-17 ruling)']);
            testCase.verifyEqual(testCase.countOf(code, 'errordlg({ME.message}'), 6, ...
                ['the six catch-ME operation-failure reporters must remain ' ...
                '(modal per the failure-reporter ruling)']);
        end

        function testWaveTwoEngineEventCountRejectionsAreModal(testCase)
            code = testCase.readCode(fullfile('subroutines', 'estimation', ...
                'psyrat_computevarcomp.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 11, ...
                ['the eleven showgui-gated engine event-count rejections ' ...
                'must be modal (ten ''Two events required'' + one ' ...
                '''Difference-score event count'')']);
        end

        function testWaveTwoDodViewerRejectionsAreModal(testCase)
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startview_dodiff.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 8, ...
                ['all eight DoD-viewer dialogs must be modal: five input ' ...
                'rejections (wave 2) plus the three catch-ME builders ' ...
                '(the 2026-08-17 failure-reporter ruling)']);
            testCase.verifyEqual(testCase.countOf(code, 'errordlg({ME.message}'), 3, ...
                'the three catch-ME reporters must remain (modal per the ruling)');
        end

        function testWaveTwoFormulaSites(testCase)
            code = testCase.readCode(fullfile('subroutines', 'filefunctions', ...
                'psyrat_formula_prompt.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 1, ...
                'the formula input rejection must be modal');
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_formula_guide.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 2, ...
                ['both formula-guide operation-failure reporters (clipboard ' ...
                '/ formula run) must be modal per the 2026-08-17 ruling']);
        end

        function testWaveTwoRelaunchSitesReorderedAndModal(testCase)
            %the seven cancel-then-relaunch sites: the relaunch must precede
            %the dialog (B20: modality binds only windows existing at the
            %dialog's creation, so a dialog fired first is covered by the
            %relaunched screen AND escapes it). The ORDER pins read the raw
            %file; the modality censuses read normalized code (B35).
            src = testCase.readSource(fullfile('subroutines', 'guis', ...
                'psyrat_startview.m'));
            code = testCase.readCode(fullfile('subroutines', 'guis', ...
                'psyrat_startview.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 1, ...
                'the results-file cancel rejection must be modal');
            testCase.verifyEqual(testCase.countRe(src, ...
                'psyrat_start;\s+errordlg\('), 1, ...
                'psyrat_startview must relaunch FIRST and create the dialog LAST');

            src = testCase.readSource(fullfile('subroutines', ...
                'installation', 'psyrat_updatepsyrat.m'));
            code = testCase.readCode(fullfile('subroutines', ...
                'installation', 'psyrat_updatepsyrat.m'));
            %four since the G87 ruling (2026-09-09): the three uiputfile cancel
            %rejections plus the 'Update unavailable' dialog raised when the
            %releases API reports nothing downloadable
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 4, ...
                ['all three cancel rejections and the update-unavailable ' ...
                'dialog in psyrat_updatepsyrat.m must be modal']);
            testCase.verifyEqual(testCase.countRe(src, ...
                'psyrat_start;\s+errordlg\(''Location not selected'',''File Error'',''modal''\);'), 3, ...
                ['psyrat_updatepsyrat.m must relaunch FIRST and create each ' ...
                'modal dialog LAST (B20 creation order)']);
            testCase.verifyEqual(testCase.countRe(src, 'psyrat_start;\s+errordlg\('), 4, ...
                ['every dialog in psyrat_updatepsyrat.m, the update-unavailable ' ...
                'one included, must follow its relaunch (B20 creation order)']);

            %psyrat_installdependents carries its three cancel rejections
            %PLUS, since G56, the three empty-wrkdir update guards (one per
            %OS arm) -- six modal dialogs, each preceded by its relaunch
            src = testCase.readSource(fullfile('subroutines', ...
                'installation', 'psyrat_installdependents.m'));
            code = testCase.readCode(fullfile('subroutines', ...
                'installation', 'psyrat_installdependents.m'));
            testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 6, ...
                ['psyrat_installdependents.m: three cancel rejections plus ' ...
                'the three G56 update guards, all modal']);
            testCase.verifyEqual(testCase.countRe(src, ...
                'psyrat_start;\s+errordlg\(''Location not selected'',''File Error'',''modal''\);'), 3, ...
                ['psyrat_installdependents.m must relaunch FIRST on each ' ...
                'cancel rejection (B20 creation order)']);
            testCase.verifyEqual(testCase.countRe(src, ...
                'psyrat_start;\s+errordlg\(dlgtext,''Update location not found'',''modal''\);'), 3, ...
                ['psyrat_installdependents.m must relaunch FIRST on each ' ...
                'G56 update guard too (B20 creation order)']);
        end

        function testNotificationSitesFollowTheirRuling(testCase)
            %the 2026-08-17 notification ruling: the four post-output D-study
            %notices KEEP their deliberate raised-last non-modal design (the
            %design comment in psyrat_criterionfigures documents why: raising
            %them last keeps the later-drawn output from burying them, and
            %modal would lock the user out of the very output the notice
            %describes) and gain distinct purpose titles. Converting any of
            %them to modal is a new ruling and must update these pins with
            %its record.
            for f = {'psyrat_criterionfigures.m', 'psyrat_relfigures.m'}
                code = testCase.readCode(fullfile('subroutines', 'plotting', f{1}));
                testCase.verifyEqual(testCase.countOf(code, ...
                    'errordlg(errorstr,''Cutoff not calculable'');'), 1, ...
                    sprintf('%s: the trlcutoff notice must be titled, non-modal', f{1}));
                testCase.verifyEqual(testCase.countOf(code, ...
                    'errordlg(errorstr,''Extrapolation beyond data'');'), 1, ...
                    sprintf('%s: the trlmax notice must be titled, non-modal', f{1}));
                %tolerant regex, not a substring: a ==0 pin whose violation
                %could otherwise hide behind a line continuation (B35, under-count)
                testCase.verifyEqual(testCase.countRe(code, testCase.ModalArg), 0, ...
                    sprintf('%s: the post-output notices are non-modal BY DESIGN', f{1}));
            end
            %the relsummary threshold notice is MODAL as of the 2026-08-18
            %ruling. It is the terminal event of a run that drew NOTHING, so
            %it belongs to the same class as the twelve operation-failure
            %reporters, not to the post-output notification class above --
            %the reason those stay non-modal (a modal would lock the user out
            %of the very output the notice describes) has no analogue here,
            %because there is no output. Traced and load-bearing: it can fire
            %only on the nogooddata routes (since G45 removed the trt_diff
            %helper's site, 1 of its 11 call sites is reachable -- the
            %sing_diff_sserr helper's; the ten guarded ones are RC-40/G39
            %expected-dead), and on
            %those routes every caller returns before drawing, so nothing is
            %created after this dialog and the B20 burial mechanism cannot
            %apply. Modality also answers the B21 risk that a non-modal box
            %does not raise above an ALREADY-open viewer on a re-run.
            code = testCase.readCode(fullfile('subroutines', 'calculation', ...
                'psyrat_relsummary.m'));
            testCase.verifyEqual(testCase.countOf(code, ...
                'errordlg(dlg, ''Data do not meet reliability threshold'',''modal'');'), 1, ...
                'the threshold notice must be titled, modal, single-site');
            testCase.verifyEqual(testCase.countRe(code, ...
                ['errordlg[ \t]*\([ \t]*dlg[ \t]*,[ \t]*''Data do not meet ' ...
                'reliability threshold''[ \t]*\)']), 0, ...
                'no non-modal spelling of the threshold notice may remain');
        end

        function testStartupUpdatePromptsAreModalFigures(testCase)
            %B44 (ruled 2026-08-22, after the click-through OBSERVED the
            %defect): the update-dependents prompt confirms a DESTRUCTIVE
            %action -- Yes runs updatedeps -> psyrat_installdependents ->
            %rmolddeps, which deletes the old dependency directories -- and
            %non-modal it could be buried while still open and live. Measured
            %on screen, not inferred: a figure census taken while the screen
            %showed NO dialog returned it as vis=on style=normal behind the
            %MATLAB desktop window, with the main menu clickable underneath a
            %pending delete confirmation.
            %
            %WHY THERE IS NO uiwait PIN HERE, unlike B28 -- and NOT for the
            %reason this comment used to give. It said the call is the LAST
            %statement in psyrat_start, so nothing runs after it to race.
            %REFUTED before the B44 commit, and the correction landed in the
            %production comment but not here until now: it is the last
            %statement only INSIDE psyrat_start, which returns to its callers,
            %and TWELVE call sites continue immediately with a modal dialog
            %of their own -- a bare "psyrat_start;" followed by a modal
            %errordlg: the 'File Error' cancel rejections in
            %psyrat_startview.m:73, psyrat_startproc.m:83,
            %psyrat_updatepsyrat (three cancel arms plus the G87 'Update
            %unavailable' refusal; lines not restated) and psyrat_installdependents
            %(three uiputfile-cancel arms; the old 112/150/183 line citation
            %had drifted and per file convention is not restated), plus the
            %three G56 'Update location not found' guards in
            %psyrat_installdependents. The B44 comment said "two", naming
            %only the GUI pair; corrected under B46, re-counted eight to
            %eleven under G56. The count matters because the LOAD-BEARING
            %clause below is about those dialogs staying modal, and the
            %wave-2 census pins them.
            %
            %uiwait is still wrong, for a different reason: it would hold
            %psyrat_start open until this prompt is answered, so the caller's
            %own "No data file was selected" message would be withheld behind
            %an unrelated update prompt. B28 needed uiwait because
            %psyrat_startproc raced on into its measurement loop with nothing
            %to stop it; here the follow-on dialogs are themselves modal and
            %MATLAB gives focus to the most recently created modal, so the
            %ordering resolves itself. LOAD-BEARING: that resolution holds
            %only while those follow-on dialogs STAY modal.
            %
            %Reads readCode, so a comment that merely QUOTES the property
            %cannot satisfy this pin (B35) -- which matters here, because the
            %ruling is documented at length in a comment at the call site.
            code = testCase.readCode('psyrat_start.m');
            testCase.verifyEqual(testCase.countRe(code, ...
                '''WindowStyle''[ \t]*,[ \t]*''modal'''), 2, ...
                ['BOTH psyrat_start update prompts must be MODAL figures: ' ...
                'the dependents one confirms deletion of the dependency ' ...
                'directories, and non-modal it was observed open, live and ' ...
                'buried behind another window (B44)']);
            %The count is 2: psyrat_ask2updatedeps AND its sibling
            %psyrat_ask2updatetoolbox, converted together on 2026-08-22. The
            %sibling's action is NOT destructive, so its justification is
            %consistency rather than blast radius -- the prompts are
            %CONFUSABLE AT THE SENTENCE LEVEL, which is B45's adjudicated
            %hedge and not a claim that the whole message is
            %indistinguishable: the opening sentence and the WARNING clause
            %differ, and always did. READ FROM SOURCE, not
            %measured on screen: until B45 they shared their second sentence;
            %they still share the sentence "Would you like to do so now?"
            %verbatim (it CLOSES the toolbox prompt but sits mid-string in
            %the dependents one, which ends on the WARNING) and a 200 px
            %figure height; and they still share the Tag psyrat_gui_update,
            %which is invisible to a user but is why the findobj lookups
            %need care. One prompt dismissable-by-burial while the other is
            %not is worse than either state alone.
            %
            %They are NOT near-identical -- that was the B45-refuted
            %overstatement, and of its three conjuncts ("wording, layout and
            %Tag") only Tag is unqualifiedly true. Differing: the opening
            %sentence, the second sentence since B45, the WARNING clause
            %(dependents only), the figure width (400 toolbox vs 500
            %dependents) and the text control (400x40 vs 500x75, same
            %order). The button-position FORMULAS are character-identical in
            %both, but figwidth differs, so the rendered buttons are not the
            %same size or place -- do not cite them as user-visible evidence.
            %B50 removed the overstatement from psyrat_start.m; this copy
            %survived HERE and is corrected under B52, which also found a
            %second present-tense copy still standing in this repo's own
            %findings ledger, fixed 2026-08-23 on the same maintainer ruling.
            %
            %A COUNT ALONE WOULD NOT HOLD THIS. Two modal figures in one file
            %could be two conversions of the same prompt, so the count is
            %paired with a per-function check below.
            blocks = strsplit(code, sprintf('\nfunction '));
            for fn = {'psyrat_ask2updatetoolbox', 'psyrat_ask2updatedeps'}
                hit = blocks(startsWith(strtrim(blocks), fn{1}));
                testCase.assertNumElements(hit, 1, ...
                    sprintf('expected exactly one %s definition', fn{1}));
                testCase.verifyEqual(numel(regexp(hit{1}, ...
                    '''WindowStyle''[ \t]*,[ \t]*''modal''', 'start')), 1, ...
                    sprintf(['%s must build its prompt as a MODAL figure ' ...
                    '(B44); the file-level count above cannot tell which ' ...
                    'of the two prompts was converted'], fn{1}));
            end
            %AND the property must not be undone after construction. Both
            %pins above read the CREATION call, so a later
            %set(h,'WindowStyle','normal') would leave them green while the
            %prompt displayed non-modally -- the same shape as the B28 trap
            %where a pin on modality alone passed on the broken form.
            testCase.verifyEqual(testCase.countRe(code, ...
                'set\([^)\n]*''WindowStyle'''), 0, ...
                ['nothing may reassign WindowStyle after construction: the ' ...
                'creation-site pins cannot see it, and the prompt would ' ...
                'display non-modally with the suite still green (B44)']);
            testCase.verifyEqual(testCase.countRe(code, ...
                '''WindowStyle''[ \t]*,[ \t]*''normal'''), 0, ...
                'no prompt in psyrat_start may declare WindowStyle normal');
        end

        function testUpdatePromptsNameTheArtifactTheyActOn(testCase)
            %B45 (ruled 2026-08-22, found on screen during the B42
            %click-through): psyrat_ask2updatedeps told the user "It is
            %recommended that you update the toolbox" while its Yes deletes
            %the DEPENDENCY directories. The sentence had been copied verbatim
            %from psyrat_ask2updatetoolbox, where it is correct.
            %
            %WHY THIS IS PINNED AT ALL for a one-word wording defect. The
            %defect ARRIVED by verbatim copy from the sibling, so the failure
            %mode is re-copying, and that is what a pin catches. Note the
            %prompts are NOT near-identical overall -- different opening
            %sentence, only the dependents one carries the WARNING clause,
            %text controls 400x40 vs 500x75 -- but they DO share this sentence
            %and their Tag, so the sentence alone did not say which prompt it
            %belonged to.
            %
            %PER FUNCTION, not per file, for the same reason the modality pin
            %above is: both prompts live in psyrat_start.m and the sibling
            %legitimately says "the toolbox", so a file-level count of either
            %phrase cannot tell which prompt it came from.
            code = testCase.readCode('psyrat_start.m');
            blocks = strsplit(code, sprintf('\nfunction '));

            deps = blocks(startsWith(strtrim(blocks), 'psyrat_ask2updatedeps'));
            testCase.assertNumElements(deps, 1, ...
                'expected exactly one psyrat_ask2updatedeps definition');
            toolbox = blocks(startsWith(strtrim(blocks), 'psyrat_ask2updatetoolbox'));
            testCase.assertNumElements(toolbox, 1, ...
                'expected exactly one psyrat_ask2updatetoolbox definition');

            %THE PIN. Tolerant spacing follows the ModalArg house style, but be
            %honest about what it buys HERE: this phrase sits inside a single
            %quoted literal, so a continuation split would insert quote
            %characters that [ \t]+ cannot bridge anyway. WHAT THIS PIN
            %ACTUALLY COVERS is the realistic regression -- a verbatim re-copy
            %of the sibling's sentence. It is case-sensitive (regexp, not
            %regexpi), and it cannot see the noun being got wrong in different
            %words. It is a re-copy guard, not a proof of correct wording.
            testCase.verifyEqual(testCase.countRe(deps{1}, ...
                'update[ \t]+the[ \t]+toolbox'), 0, ...
                ['the dependents prompt must not tell the user to update ' ...
                'the TOOLBOX (B45): its Yes deletes the dependency ' ...
                'directories, and this sentence was copied verbatim from ' ...
                'the sibling prompt, whose Tag it still shares']);

            %POSITIVE CONTROLS, and they are not decoration. The assertion
            %above does fire on the regression it exists for -- a re-copied
            %sentence trips it. What it CANNOT do is notice that it read
            %nothing: a readCode returning '' or a block split matching
            %nothing satisfies a ==0 without reading a character of the
            %prompt. These two are what make that case fail.
            testCase.verifyEqual(testCase.countRe(deps{1}, ...
                'update[ \t]+the[ \t]+dependents'), 1, ...
                ['the dependents prompt must name the dependents, and this ' ...
                'also proves the block split and readCode above actually ' ...
                'reached the prompt text']);
            testCase.verifyEqual(testCase.countRe(toolbox{1}, ...
                'update[ \t]+the[ \t]+toolbox'), 1, ...
                ['the sibling toolbox prompt is CORRECT as-is and must keep ' ...
                'saying "the toolbox" -- without this, the B45 pin could be ' ...
                'satisfied by deleting the wording from both prompts']);
        end

        function testStartupTeardownIsScopedToPsyratFigures(testCase)
            %B48 (ruled 2026-08-22): psyrat_start's four startup callbacks
            %ran a bare `close all`, which closes EVERY figure in the session
            %whose HandleVisibility is on, whoever owns it. All four fire
            %during startup prompts -- precisely when a user is likely to
            %have EEGLAB/ERPLAB figures open from earlier work -- so
            %accepting an update, and even DECLINING the install prompt via
            %giveup, silently destroyed their unsaved figures.
            %
            %WHY readCode AND NOT readSource. The production fix leaves the
            %words "close all" in the helper's own docstring, which explains
            %both what it replaced and the one class of window it no longer
            %reaches. A raw substring scan would read those as live calls and
            %this pin would fail on correct code. Comment stripping is
            %load-bearing here, not incidental. (An earlier version of this
            %note also credited the B44 comment at psyrat_ask2updatetoolbox --
            %wrong: the same diff rewrote that line and it no longer contains
            %the phrase. Both surviving occurrences are in the docstring.)
            code = testCase.readCode('psyrat_start.m');

            testCase.verifyEqual(testCase.countRe(code, 'close[ \t]+all'), 0, ...
                ['psyrat_start must not run a bare `close all`: it destroys ' ...
                'the user''s own unsaved figures, not just PsyRAT''s (B48)']);

            %POSITIVE CONTROLS. The ==0 above is satisfied by a readCode that
            %returned '' -- or by deleting the teardown altogether, which
            %would strand a MODAL prompt across the whole installer run. The
            %B45 lesson applies unchanged: a negative assertion with no
            %positive beside it cannot distinguish "fixed" from "not read".
            testCase.verifyEqual(testCase.countRe(code, ...
                'psyrat_closepsyratfigs;'), 4, ...
                ['all four startup callbacks must still tear the PsyRAT GUI ' ...
                'down, and this also proves readCode reached them']);
            testCase.verifyEqual(testCase.countRe(code, ...
                'findall\(0,''Type'',''figure'',''-regexp'',''Tag'',''\^psyrat_'''), 1, ...
                ['the teardown must be SCOPED by psyrat_ Tag prefix -- ' ...
                'without this the call count above is satisfied by a helper ' ...
                'whose body is itself `close all`']);

            %AND THE HANDLES IT FINDS MUST ACTUALLY BE CLOSED. Every assertion
            %above reads only the SELECTOR, so deleting the close() leaves a
            %helper that finds the windows and returns -- with this class
            %still green and the MODAL psyrat_gui_update left up for the whole
            %installer run. Verified by mutation: with the close removed, all
            %five assertions above still passed. A selector pin is not a
            %teardown pin.
            testCase.verifyEqual(testCase.countRe(code, ...
                'close\(psyratfigs\)'), 1, ...
                ['the teardown must CLOSE the handles it finds (B48); the ' ...
                'selector pin alone is satisfied by a no-op helper']);

            %PER FUNCTION, not per file. A file-level count of 4 is satisfied
            %by four calls inside one callback and none in the other three,
            %which is exactly the shape B42 had (an if/elseif chain that
            %covered one dependent while the deletion covered three).
            blocks = strsplit(code, sprintf('\nfunction '));
            for name = {'updatepsyrat', 'updatedeps', 'giveup', 'installdeps'}
                blk = blocks(startsWith(strtrim(blocks), [name{1} '(']));
                testCase.assertNumElements(blk, 1, ...
                    sprintf('expected exactly one %s definition', name{1}));
                testCase.verifyEqual(testCase.countRe(blk{1}, ...
                    'psyrat_closepsyratfigs;'), 1, ...
                    sprintf(['%s must scope its teardown to PsyRAT''s own ' ...
                    'windows (B48)'], name{1}));
            end
        end

        function testTeardownReachableWindowsCarryNamespacedTags(testCase)
            %THE OTHER HALF OF B48, and the half that is easy to lose. Once
            %the teardown matches on a 'psyrat_' Tag PREFIX, an untagged
            %PsyRAT window is no longer swept up -- it is stranded.
            %
            %NAMED "REACHABLE", NOT "EVERY", ON PURPOSE. This does not assert
            %that every shipped PsyRAT figure is prefixed, because one class
            %is not: MATLAB tags msgbox/warndlg/errordlg windows
            %'Msgbox_<title>', which the teardown deliberately does not match
            %(that prefix is MATLAB's, not PsyRAT's, so matching it would
            %close EEGLAB's dialogs too -- see psyrat_closepsyratfigs). A test
            %called "every" would be a false universal, which is exactly the
            %class of defect the suites cannot see.
            %
            %THE TRACE-PLOT WINDOW USED TO BE A SECOND EXCEPTION and is not
            %any more. It was excused on the grounds that it cannot outlive
            %estimation; that was REFUTED -- dismissing the rerun prompt with
            %the window X makes guidata throw in psyrat_computevarcompwarp
            %before its own close/delete of that figure is reached. It is now
            %tagged 'psyrat_tplots' and is covered like everything else.
            %
            %Four were
            %found untagged when the fix was scoped, and one of them is
            %reachable at the teardown sites by an ordinary route:
            %psyrat_startimportscore's localBackToStart closes only
            %psyrat_import_gui before calling psyrat_start, and
            %psyrat_import_trialbrowser.m contains no close, delete or
            %CloseRequestFcn anywhere, so nothing in the toolbox would ever
            %have closed the browser again.
            %
            %NOT A STYLE PIN. Each assertion below is the only thing standing
            %between a specific window and being orphaned for the session.
            expected = {
                'subroutines/guis/psyrat_import_trialbrowser.m', 'psyrat_gui_trialbrowser'
                'subroutines/guis/psyrat_startview_dynrel.m',    'psyrat_gui_dynrelview'
                'subroutines/guis/psyrat_startview_splits.m',    'psyrat_gui_splitsview'
                };

            for i = 1:size(expected, 1)
                code = testCase.readCode(expected{i,1});
                testCase.verifyGreaterThanOrEqual(testCase.countOf(code, ...
                    ['''' expected{i,2} '''']), 1, ...
                    sprintf(['%s must tag its figure %s, or psyrat_start''s ' ...
                    'prefix-scoped teardown leaves it stranded (B48)'], ...
                    expected{i,1}, expected{i,2}));
            end

            %THE DYNREL AND SPLITS VIEWERS MUST NOT REUSE 'psyrat_gui'. Both
            %files look that Tag up to find and close the SETUP screen before
            %building their own window (psyrat_startview_dynrel.m:106,
            %psyrat_startview_splits.m:62). Tagging the viewer with the same
            %name would make those lookups ambiguous, so the distinct names
            %above are required rather than cosmetic. Asserted as a count of
            %the lookup, which must stay at exactly one per file.
            for f = {'subroutines/guis/psyrat_startview_dynrel.m', ...
                    'subroutines/guis/psyrat_startview_splits.m'}
                code = testCase.readCode(f{1});
                testCase.verifyEqual(testCase.countOf(code, ...
                    'findobj(''Tag'',''psyrat_gui'')'), 1, ...
                    sprintf(['%s must keep exactly one setup-screen lookup; ' ...
                    'its own viewer window carries a distinct Tag'], f{1}));
            end

            %B49, the same defect class one file over. The Command Line Tools
            %wait prompt was tagged 'f' and closed with findobj('Tag','f'),
            %which searches the whole graphics root and is NOT restricted to
            %figures -- any object in the user's session tagged 'f',
            %including a uicontrol inside somebody else's figure, matched and
            %was closed. It was masked in practice because installdeps ran
            %`close all` first.
            inst = testCase.readCode( ...
                'subroutines/installation/psyrat_installdependents.m');
            testCase.verifyEqual(testCase.countRe(inst, ...
                '''Tag''[ \t]*,[ \t]*''f'''), 0, ...
                ['the installer must not look figures up by the ' ...
                'one-character Tag ''f'' (B49)']);
            testCase.verifyEqual(testCase.countRe(inst, ...
                '\.Tag[ \t]*=[ \t]*''f'''), 0, ...
                'the installer must not TAG a figure ''f'' either (B49)');

            %POSITIVE CONTROL for the two ==0 assertions above: the namespaced
            %name must appear exactly twice -- the assignment and the lookup.
            %Without this, deleting the prompt's teardown entirely would pass.
            testCase.verifyEqual(testCase.countOf(inst, ...
                '''psyrat_gui_cltwait'''), 2, ...
                ['the Command Line Tools prompt must be tagged AND looked ' ...
                'up by its namespaced Tag, and this proves the two ==0 ' ...
                'assertions above read a real file']);
        end

        function testNormalizerCountsGenuineSitesWhateverTheirLayout(testCase)
            %POSITIVE CONTROLS for normalizeSource (B35). Every "must not be
            %counted" case in the sibling test needs one of these standing
            %beside it: a normalizer that simply returned '' would satisfy
            %every negative assertion on its own and prove nothing.
            %
            %EVERY FIXTURE HERE WAS MUTATION-CHECKED. The implementation was
            %broken five ways in turn -- blank-skipping lookback, no
            %transpose rule, no closing-quote exclusion, first-percent-wins,
            %no continuation joining -- and each mutant is killed by at
            %least one assertion below. Two earlier drafts of this method
            %were VACUOUS and were rewritten or deleted; if a fixture here
            %stops discriminating, delete it rather than leave decoration.
            testCase.verifyEqual(testCase.countModalSites( ...
                'errordlg(m,''T'',''modal'');'), 1, ...
                'the plain one-line call form must still be counted');

            %THE UNDER-COUNT DEFECT ITSELF. This is the layout the first B33
            %draft used, which the raw substring scan read as zero.
            testCase.verifyEqual(testCase.countModalSites( ...
                sprintf('errordlg(m,''T'', ...\n    ''modal'');')), 1, ...
                ['a modal argument reached across a line continuation must ' ...
                'be counted; missing it is how a dialog could be made ' ...
                'NON-modal without any census moving (B35, under-count)']);

            %models psyrat_formula_prompt.m:302 verbatim in shape. A "first
            %percent on the line wins" stripper deletes this real dialog
            %site and drives that file's pin from 1 to 0.
            testCase.verifyEqual(testCase.countModalSites( ...
                'errordlg(sprintf(''Input %d is not numeric.'',i),''T'',''modal'');'), 1, ...
                ['a percent sign INSIDE a quoted literal is not a comment; ' ...
                'treating it as one silently deletes a real dialog site']);

            testCase.verifyEqual(testCase.countModalSites( ...
                's = ''a...b''; errordlg(m,''T'',''modal'');'), 1, ...
                ['an ellipsis inside a quoted literal is not a line ' ...
                'continuation and must not swallow the rest of the line']);

            %BRACKET CONCATENATION -- the shape that broke an earlier draft.
            %[ident 'text'] is concatenation, NOT a transpose followed by a
            %literal. Reading it as a transpose flips the parity of every
            %quote after it, so the line is cut at the percent INSIDE the
            %literal and the dialog vanishes from the census. Shapes like
            %this are live in the scanned files; without this fixture the
            %defect ships green. The percent is load-bearing.
            testCase.verifyEqual(testCase.countModalSites( ...
                'errordlg([msg '' is 50% done''],''T'',''modal'');'), 1, ...
                ['bracket concatenation must not be read as a transpose: ' ...
                'that silently deletes the rest of the line']);

            %the same defect in its second mechanism: a spuriously opened
            %literal also MASKS the trailing ..., so the continuation is
            %never joined and the modal argument is never reached
            testCase.verifyEqual(testCase.countModalSites( ...
                sprintf('errordlg([a '' | '' b],ttl, ...\n    ''modal'');')), 1, ...
                ['concatenation followed by a continuation must still join: ' ...
                'a mis-read quote masks the ... and suppresses the join']);

            %pins the transpose rule itself. Under "no transpose rule" the
            %quote after a opens a literal, parity flips, and the percent
            %below is read as a comment start -> 0.
            testCase.verifyEqual(testCase.countModalSites( ...
                'y = a''; errordlg(''50% done'',''T'',''modal'');'), 1, ...
                'a transpose must not be mistaken for the start of a literal');

            %pins the closing-quote exclusion specifically: with a'' the
            %second quote follows a TRANSPOSE (unmarked), not a literal
            %close (marked), so dropping the mask clause misreads it
            testCase.verifyEqual(testCase.countModalSites( ...
                'y = a''''; errordlg(''50% done'',''T'',''modal'');'), 1, ...
                'a repeated transpose must not be mistaken for a literal');

            testCase.verifyEqual(testCase.countModalSites( ...
                'x = [''a'' ''b%c'']; errordlg(m,''T'',''modal'');'), 1, ...
                ['two adjacent literals on one line must not leave the scan ' ...
                'believing it is still inside a literal']);

            testCase.verifyEqual(testCase.countModalSites( ...
                'errordlg(''it''''s 50% bad'',''T'',''modal'');'), 1, ...
                'an escaped quote inside a literal must not close it');
        end

        function testNormalizerIgnoresCommentsAndBlockDelimiters(testCase)
            %NEGATIVE CONTROLS for normalizeSource (B35). Read these paired
            %with the positive controls above -- alone they are satisfiable
            %by a normalizer that discards everything.
            needle = [',' '''' 'modal' '''' ')'];

            testCase.verifyEqual(testCase.countModalSites( ...
                sprintf('%%prose mentioning %s in passing\nx = 1;', needle)), 0, ...
                ['a full-line comment quoting the modal argument must not ' ...
                'be counted; one written to document this very pin pushed ' ...
                'the startproc census to 40 and failed the suite ' ...
                '(B35, over-count)']);

            testCase.verifyEqual(testCase.countModalSites( ...
                sprintf('x = 1; %%was errordlg(m,''T'',''modal'');')), 0, ...
                'a trailing comment quoting the modal argument must not be counted');

            testCase.verifyEqual(testCase.countModalSites( ...
                sprintf('y = a''; %%note %s here', needle)), 0, ...
                ['a trailing comment must still be stripped on a line that ' ...
                'also contains a transpose']);

            %the over-count face of the bracket-concatenation defect: a
            %mis-read quote leaves the trailing comment UNSTRIPPED, so a
            %comment merely describing a dialog inflates a census
            testCase.verifyEqual(testCase.countModalSites( ...
                'labels{1} = [g{1} '' - '' e{1}]; %was errordlg(m,''T'',''modal'');'), 0, ...
                ['a trailing comment on a bracket-concatenation line must ' ...
                'still be stripped']);

            %a block comment would be stripped one delimiter at a time by a
            %per-line scan, leaving the commented-out CODE in the census.
            %There are none in the scanned files today; this makes the day
            %one appears loud instead of silent.
            caught = '';
            try
                testCase.normalizeSource(sprintf('%%{\nerrordlg(m,''T'',''modal'');\n%%}'));
            catch ME
                caught = ME.identifier;
            end
            testCase.verifyEqual(caught, 'TestModalRejectionDialogs:blockComment', ...
                ['a block comment must raise rather than be silently ' ...
                'treated as two ordinary comment lines around live code']);
        end

    end

    methods (Access = private)

        function src = readSource(~, relpath)
            %absolute path from THIS file's location, never from pwd or the
            %MATLAB path (worktree copies have shadowed the toolbox before)
            here = fileparts(mfilename('fullpath'));
            src = fileread(fullfile(fileparts(here), relpath));
        end

        function code = readCode(testCase, relpath)
            %the same file as readSource, reduced to CODE: comments removed
            %and line continuations joined. Every count in this class reads
            %this; only the ordering pins read readSource (B35).
            code = testCase.normalizeSource(testCase.readSource(relpath));
        end

        function n = countOf(~, src, needle)
            n = numel(strfind(src, needle));
        end

        function n = countRe(~, src, pat)
            n = numel(regexp(src, pat, 'start'));
        end

        function n = countModalSites(testCase, srctext)
            %count modal dialog sites in a literal source STRING. Used only
            %by the normalizer fixtures, and deliberately routed through the
            %same normalizeSource + ModalArg pair the censuses use, so a
            %fixture cannot pass against machinery the pins do not share.
            n = testCase.countRe(testCase.normalizeSource(srctext), ...
                testCase.ModalArg);
        end

        function code = normalizeSource(testCase, src)
            %Return SRC with comments removed and ... continuations joined,
            %leaving quoted literals verbatim.
            %
            %WHY PER LINE. A MATLAB quoted literal cannot span a newline,
            %so literal state always resets at a line break and each
            %physical line can be scanned independently; the majority that
            %contain no quote at all skip the scan entirely. Note the limit
            %of that guarantee: it makes the LITERAL state exact, while the
            %transpose-vs-literal call in literalMask remains a left-
            %adjacency rule rather than a parse. It is correct on every
            %construct in the scanned tree and pinned by fixture, but it is
            %a rule, not a theorem -- see literalMask.
            src = strrep(src, char(13), '');
            lines = strsplit(src, newline, 'CollapseDelimiters', false);
            out = cell(1, numel(lines));
            nout = 0;
            pending = '';
            for k = 1:numel(lines)
                trimmed = strtrim(lines{k});
                if strcmp(trimmed, '%{') || strcmp(trimmed, '%}')
                    error('TestModalRejectionDialogs:blockComment', ...
                        ['line %d opens or closes a block comment. A ' ...
                        'per-line scan strips the delimiter and keeps the ' ...
                        'block''s CODE, which would inflate every census ' ...
                        'here. Teach normalizeSource about block comments ' ...
                        'before introducing one.'], k);
                end
                [linecode, continued] = testCase.codePart(lines{k});
                if continued
                    pending = [pending linecode ' ']; %#ok<AGROW>
                else
                    nout = nout + 1;
                    out{nout} = [pending linecode];
                    pending = '';
                end
            end
            if ~isempty(pending)
                nout = nout + 1;
                out{nout} = pending;
            end
            code = strjoin(out(1:nout), newline);
        end

        function [linecode, continued] = codePart(testCase, txt)
            %Resolve ONE physical line into its code part, dropping any
            %comment, and report whether it ends in a ... continuation.
            %
            %A percent sign inside a quoted literal is NOT a comment. That
            %is not hypothetical: psyrat_formula_prompt.m:302 passes a
            %sprintf format string to errordlg, and a naive "first percent
            %wins" split would delete that real dialog site and drive its
            %pin from 1 to 0.
            continued = false;
            if ~any(txt == '''') && ~any(txt == '"')
                %fast path -- no literal can be open on this line, so any
                %percent starts a comment and any ellipsis continues. Most
                %lines of the scanned files land here.
                cut = find(txt == '%', 1);
                cont = strfind(txt, '...');
                if ~isempty(cont) && (isempty(cut) || cont(1) < cut)
                    linecode = txt(1:cont(1)-1);
                    continued = true;
                    return
                end
                if isempty(cut)
                    linecode = txt;
                else
                    linecode = txt(1:cut-1);
                end
                return
            end
            inliteral = testCase.literalMask(txt);
            cut = find(txt == '%' & ~inliteral, 1);
            if isempty(cut)
                linecode = txt;
            else
                linecode = txt(1:cut-1);
            end
            for i = 1:numel(linecode)-2
                if strcmp(linecode(i:i+2), '...') && ~any(inliteral(i:i+2))
                    linecode = linecode(1:i-1);
                    continued = true;
                    return
                end
            end
        end

        function mask = literalMask(~, txt)
            %MASK(i) is true where TXT(i) lies inside a quoted literal, the
            %delimiting quotes included.
            %
            %TRANSPOSE vs LITERAL is the one genuine ambiguity in MATLAB's
            %lexer that a scanner has to resolve by rule. Outside a literal,
            %a single quote is a TRANSPOSE when the character IMMEDIATELY
            %before it could end a value -- an identifier character, digit,
            %closing bracket, dot, or an earlier transpose -- and opens a
            %LITERAL otherwise. Adjacency is the whole rule: MATLAB reads
            %a' as a transpose but [a 'b'] as concatenation with a string,
            %and only the intervening blank distinguishes them. Excluding a
            %preceding literal-CLOSING quote is what additionally keeps
            %a'' correct: that closing quote is marked in MASK, an earlier
            %transpose is not.
            %
            %DO NOT "TIDY" THIS BY SKIPPING BLANKS BACKWARDS. An earlier
            %draft did, and adversarial review killed it before it landed.
            %Skipping blanks erases exactly the distinction above, so
            %[msg ' is 50% done'] read as a transpose, the parity of every
            %quote after it flipped, and the line was cut at a percent sign
            %sitting INSIDE a literal -- deleting real code. Shapes like
            %that are live in the scanned files, and the ==0 pins would
            %have gone on passing with a modal dialog present.
            %
            %FAILURE DIRECTIONS. Calling a real opening quote a transpose
            %is the dangerous one: it silently LOWERS a census, and it also
            %masks a following ... so the continuation is never joined --
            %two silent under-counts, not a loud over-count. It is pinned
            %by the bracket-concatenation fixtures. The opposite direction
            %leaves a comment unstripped and fails loudly against a
            %hardcoded count.
            %
            %The escape branch below is DEFENSIVE, not pinned, and a
            %fixture claiming otherwise would be decoration: because a
            %doubled quote is adjacent, re-pairing it as two literals
            %yields the SAME set of in-literal characters, so dropping the
            %branch cannot change any count. It is kept so MASK's
            %boundaries are faithful for any future consumer.
            n = numel(txt);
            mask = false(1, n);
            valueend = ['abcdefghijklmnopqrstuvwxyz' ...
                'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_)]}.'];
            quote = char(0);
            i = 1;
            while i <= n
                c = txt(i);
                if quote ~= char(0)
                    mask(i) = true;
                    if c == quote
                        if i < n && txt(i+1) == quote
                            %doubled quote inside a literal is an escape,
                            %not a close
                            mask(i+1) = true;
                            i = i + 2;
                            continue
                        end
                        quote = char(0);
                    end
                    i = i + 1;
                    continue
                end
                if c == ''''
                    %IMMEDIATELY preceding character -- do NOT skip blanks.
                    %Blank-skipping here was a real defect, caught by
                    %adversarial review before this landed: it made
                    %[msg ' text'] read as a transpose, flipped the parity
                    %of every quote after it, and silently deleted code.
                    j = i - 1;
                    istranspose = j >= 1 && (any(txt(j) == valueend) || ...
                        (txt(j) == '''' && ~mask(j)));
                    if istranspose
                        i = i + 1;
                        continue
                    end
                    quote = '''';
                    mask(i) = true;
                    i = i + 1;
                    continue
                end
                if c == '"'
                    quote = '"';
                    mask(i) = true;
                    i = i + 1;
                    continue
                end
                i = i + 1;
            end
        end

    end
end
