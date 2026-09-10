classdef TestScoringSpecBuilder < PsyRATTestBase
    % Tests for the shared scoring-method registry (psyrat_scoring_methods) and
    % the pure GUI spec builder (psyrat_scoring_spec_from_inputs). The builder
    % holds the string-to-spec parsing that the import/scoring GUI used to do
    % inline, so it can be unit-tested without a figure. The GUI callbacks
    % themselves (uigetfile/inputdlg/listdlg) are not headlessly testable and
    % are exercised by these helpers plus manual review.

    methods (Test)

        function testMethodsRegistryIncludesNewMethods(testCase)
            m = psyrat_scoring_methods();
            testCase.verifyTrue(iscellstr(m)); %#ok<ISCLSTR>
            testCase.verifyEqual(numel(m),7);
            testCase.verifyEqual(numel(unique(m)),7); % no duplicates
            testCase.verifyTrue(ismember('fractional_area_latency',m));
            testCase.verifyTrue(ismember('integrated_area_amplitude',m));
        end

        function testBuilderAbsolutePeakSingleChannel(testCase)
            in = struct('name','P3','window_start_ms','250','window_end_ms','500', ...
                'channel','1','polarity','positive');
            spec = psyrat_scoring_spec_from_inputs('absolute_peak_amplitude',in);
            testCase.verifyEqual(spec.method,'absolute_peak_amplitude');
            testCase.verifyEqual(spec.channel_selector,1);
            testCase.verifyEqual(spec.polarity,'positive');
            testCase.verifyEqual(spec.window_start_ms,250);
            testCase.verifyEqual(spec.window_end_ms,500);
        end

        function testBuilderChannelRoiLabelsParsed(testCase)
            in = struct('name','ROI','window_start_ms','0','window_end_ms','400', ...
                'channel','Fz,Cz,Pz','polarity','positive');
            spec = psyrat_scoring_spec_from_inputs('time_window_mean_amplitude',in);
            testCase.verifyEqual(spec.channel_selector,{'Fz','Cz','Pz'});
        end

        function testBuilderNumericRoiParsed(testCase)
            in = struct('name','ROI','window_start_ms','0','window_end_ms','400', ...
                'channel','1,2,3','polarity','positive');
            spec = psyrat_scoring_spec_from_inputs('time_window_mean_amplitude',in);
            testCase.verifyEqual(spec.channel_selector,[1 2 3]);
        end

        function testBuilderFractionalAreaLatencyFractionParsed(testCase)
            in = struct('name','FAL','window_start_ms','0','window_end_ms','400', ...
                'channel','1','polarity','positive','area_fraction','0.25');
            spec = psyrat_scoring_spec_from_inputs('fractional_area_latency',in);
            testCase.verifyEqual(spec.method_params.area_fraction,0.25);
        end

        function testBuilderFractionalAreaLatencyBlankFractionDefaults(testCase)
            % Blank area_fraction -> spec_create applies the 0.5 default.
            in = struct('name','FAL','window_start_ms','0','window_end_ms','400', ...
                'channel','1','polarity','positive','area_fraction','');
            spec = psyrat_scoring_spec_from_inputs('fractional_area_latency',in);
            testCase.verifyEqual(spec.method_params.area_fraction,0.5);
        end

        function testBuilderIntegratedAreaAmplitude(testCase)
            in = struct('name','Area','window_start_ms','0','window_end_ms','400', ...
                'channel','1','polarity','negative');
            spec = psyrat_scoring_spec_from_inputs('integrated_area_amplitude',in);
            testCase.verifyEqual(spec.method,'integrated_area_amplitude');
            testCase.verifyEqual(spec.polarity,'negative');
        end

        function testBuilderLocalPeakNeighborhood(testCase)
            in = struct('name','LP','window_start_ms','0','window_end_ms','400', ...
                'channel','1','polarity','positive','neighborhood','5');
            spec = psyrat_scoring_spec_from_inputs('local_peak_amplitude',in);
            testCase.verifyEqual(spec.method_params.neighborhood_samples,5);

            % Blank neighborhood falls back to the conventional default of 3.
            in.neighborhood = '';
            spec2 = psyrat_scoring_spec_from_inputs('local_peak_amplitude',in);
            testCase.verifyEqual(spec2.method_params.neighborhood_samples,3);
        end

        function testBuilderCentroidWeightMode(testCase)
            in = struct('name','Cen','window_start_ms','0','window_end_ms','400', ...
                'channel','1','polarity','positive','weight_mode','positive_only');
            spec = psyrat_scoring_spec_from_inputs('centroid_latency',in);
            testCase.verifyEqual(spec.method_params.weight_mode,'positive_only');
        end

        function testBuilderEmptyChannelErrors(testCase)
            in = struct('name','x','window_start_ms','0','window_end_ms','400', ...
                'channel','','polarity','positive');
            testCase.verifyError(@() psyrat_scoring_spec_from_inputs('absolute_peak_amplitude',in), ...
                'psyrat_scoring:channels');
        end

    end
end
