classdef TestProjectionReportContract < PsyRATTestBase
    % Public report/export contracts, independent of the retention repair.
    methods (Test)
        function case24ObservedReportCarriesDesign(t)
            pd=PsyRATTestDataFactory.makeSplitsTRTSummaryData();
            pd.rel.data(1:4,:)=[]; % one missing occasion
            r=psyrat_report(pd);
            t.assertTrue(isfield(r,'splits_design'));
            t.verifyEqual(r.splits_design.mode,'observed_available_occasions');
            t.verifyTrue(r.splits_design.missing_occasions);
            t.verifySubstring(r.design_note,'available occasion');
            t.verifyFalse(contains(r.gbias_note,'downward-biased'));
        end
        function explicitOverridesReachReportsAndSavedFiles(t)
            pd=PsyRATTestDataFactory.makeSplitsTRTSummaryData();
            observed=psyrat_report(pd);
            for ext={'.csv','.xlsx'}
                folder=tempname; mkdir(folder); t.addTeardown(@()rmdir(folder,'s'));
                r=psyrat_report(pd,'obs',2,'nocc',1,'outdir',folder,'format',ext{1});
                t.assertTrue(isfield(r,'splits_design'));
                t.verifyEqual(r.splits_design.mode,'balanced_projection');
                t.verifyEqual(r.splits_design.projection_nsplit,2);
                t.verifyEqual(r.splits_design.projection_nocc,1);
                t.verifyNotEqual(r.tables.coefficients,observed.tables.coefficients);
                for k=1:numel(r.files)
                    txt=localWords(r.files{k});
                    t.verifySubstring(txt,'balanced_projection');
                    t.verifySubstring(txt,'projection');
                end
            end
        end
        function numericProjectionNoteAndExports(t)
            pd=localProjected();
            for ext={'.csv','.xlsx'}
                folder=tempname; mkdir(folder); t.addTeardown(@()rmdir(folder,'s'));
                r=psyrat_report(pd,'outdir',folder,'format',ext{1});
                t.verifySubstring(r.cutoff_note,'extrapolation beyond the data');
                t.verifySubstring(r.cutoff_note,'projected coefficient');
                t.verifyFalse(contains(r.cutoff_note,'print -1'));
                for k=1:numel(r.files)
                    txt=localWords(r.files{k});
                    t.verifySubstring(txt,'extrapolation beyond the data');
                end
            end
        end
        function emptyRetainedReportDisclosesReferenceProvenance(t)
            pd=localProjected();
            pd.relsummary.group.goodids={};
            pd.relsummary.group.badids={'s1';'s2';'s3'};
            pd.relsummary.group.goodn=0; pd.relsummary.group.badn=3;
            pd.relsummary.group.event.goodn=0;
            pd.relsummary.group.event.summary_population='full_sample_reference';
            r=psyrat_report(pd);
            t.assertTrue(isfield(r,'retention_note'));
            t.verifySubstring(r.retention_note,'full-sample reference');
            t.verifySubstring(r.retention_note,'no participants retained');
        end
        function cleanControl(t)
            pd=PsyRATTestDataFactory.makeICSummaryData();
            r=psyrat_report(pd);
            t.verifyEmpty(r.cutoff_note);
        end
        function oldProjectionSentinelStillExplained(t)
            pd=localProjected();
            pd.relsummary.group.event=rmfield(pd.relsummary.group.event,'cutoff_extrapolated');
            pd.relsummary.group.event.rel.m=-1;
            r=psyrat_report(pd);
            t.verifySubstring(r.cutoff_note,'print -1');
            t.verifySubstring(r.cutoff_note,'extrapolation beyond the data');
        end
    end
end

function pd=localProjected()
pd=PsyRATTestDataFactory.makeICSummaryData();
pd.relsummary.group.event.trlcutoff=81;
pd.relsummary.group.event.rel=struct('m',.9,'ll',.9,'ul',.9);
pd.relsummary.group.event.cutoff_extrapolated=true;
pd.relsummary.group.event.observed_trlinfo=struct('min',6,'max',6,'mean',6,'med',6,'std',0);
end

function txt=localWords(file)
if endsWith(file,'.csv'); txt=fileread(file);
else
    cells=readcell(file);
    txt=strjoin(string(cells(cellfun(@ischar,cells))),newline);
end
end
