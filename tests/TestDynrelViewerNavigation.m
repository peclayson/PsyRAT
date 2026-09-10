classdef TestDynrelViewerNavigation < PsyRATTestBase
    %Behavioral pins for finding B18 -- the dynamic-reliability viewer must
    %carry the two standard navigation buttons every sibling viewer has
    %('Back to Home', 'Open Another .psyrat File'), wired through the same
    %local-callback idiom (psyrat_svb / psyrat_view_loadnewfile). Without
    %them the setup GUI is already closed on entry, so Done drops the user
    %to a bare MATLAB prompt with every output figure still open and no way
    %back or to another stored result.
    %
    %The buttons' presence, wiring, and on-canvas geometry are pinned
    %headlessly (the viewer builds headlessly from the factory fixture; no
    %CmdStan). CLICKING them relaunches psyrat_start / psyrat_startview, so
    %the click itself is deferred to the live-GUI session (B20 precedent);
    %these pins hold the structure until then.

    methods (Test)

        function testNavigationButtonsExistAndAreWired(testCase)
            pd = PsyRATTestDataFactory.makeDynrelDiffSummaryData();
            name = 'Dynamic Reliability: View Results';
            close(findall(0, 'Type', 'figure', 'Name', name));
            psyrat_startview_dynrel('psyrat_prefs', psyrat_defaults(), ...
                'psyrat_data', pd);
            fig = findall(0, 'Type', 'figure', 'Name', name);
            testCase.assertNotEmpty(fig, 'Dynrel viewer figure not created');
            fig = fig(1);
            cleaner = onCleanup(@() close(fig, 'force'));

            back = testCase.buttonByString(fig, 'Back to Home');
            testCase.assertNotEmpty(back, ...
                'the dynrel viewer must carry a Back to Home button (B18)');
            testCase.verifyEqual(func2str(back.Callback{1}), 'psyrat_svb', ...
                'Back to Home must route through the sibling psyrat_svb idiom');
            testCase.verifyEqual(back.Callback{2}, fig, ...
                'psyrat_svb must receive the viewer figure to close');

            openbtn = testCase.buttonByString(fig, 'Open Another');
            testCase.assertNotEmpty(openbtn, ...
                'the dynrel viewer must carry an Open Another .psyrat File button (B18)');
            testCase.verifyEqual(func2str(openbtn.Callback{1}), ...
                'psyrat_view_loadnewfile', ...
                'Open Another must route through psyrat_view_loadnewfile');
            testCase.verifyEqual(openbtn.Callback{2}, fig, ...
                'psyrat_view_loadnewfile must receive the viewer figure to close');

            %geometry sanity: both buttons fully on the canvas. The real
            %look is the live click-through's job (headless geometry passes
            %at equality but misses visual clipping); this catches only an
            %off-canvas row.
            figpos = fig.Position;
            for h = [back openbtn]
                pos = h.Position;
                testCase.verifyGreaterThanOrEqual(pos(2), 0, ...
                    'navigation button sits below the figure bottom');
                testCase.verifyLessThanOrEqual(pos(2) + pos(4), figpos(4), ...
                    'navigation button clipped by the figure top');
            end
        end

    end

    methods (Access = private)

        function h = buttonByString(~, fig, fragment)
            %find a pushbutton whose String (joined, for multi-line cell
            %labels) contains the fragment; empty when absent
            h = [];
            btns = findall(fig, 'Style', 'pushbutton');
            for b = btns(:)'
                s = b.String;
                if iscell(s); s = strjoin(s, ' '); end
                if contains(string(s), fragment)
                    h = b;
                    return;
                end
            end
        end

    end
end
