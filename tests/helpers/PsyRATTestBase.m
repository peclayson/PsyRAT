classdef (Abstract) PsyRATTestBase < matlab.unittest.TestCase
    % Shared setup/teardown for deterministic PsyRAT unit tests.

    properties (Access = private)
        OriginalPath
        OriginalRng
        OriginalFigureVisible
    end

    methods (TestMethodSetup)
        function setupEnvironment(testCase)
            root = testCase.projectRoot();

            testCase.OriginalPath = path;
            testCase.OriginalRng = rng;
            testCase.OriginalFigureVisible = get(groot, 'DefaultFigureVisible');

            rng(42, 'twister');
            set(groot, 'DefaultFigureVisible', 'off');

            addpath(root);
            addpath(genpath(fullfile(root, 'subroutines')));
            addpath(genpath(fullfile(root, 'tests', 'helpers')));
        end
    end

    methods (TestMethodTeardown)
        function teardownEnvironment(testCase)
            close all force;
            path(testCase.OriginalPath);
            rng(testCase.OriginalRng);
            set(groot, 'DefaultFigureVisible', testCase.OriginalFigureVisible);
        end
    end

    methods
        function root = projectRoot(~)
            thisFile = mfilename('fullpath');
            root = fileparts(fileparts(fileparts(thisFile)));
        end
    end
end
