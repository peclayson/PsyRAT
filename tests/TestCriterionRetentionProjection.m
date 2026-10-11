classdef TestCriterionRetentionProjection < PsyRATTestBase
    % G84/G88 consumer regression: real summaries and actual saved tables.
    % Fixed P=1, residual=9, mu=cut=0 gives criterion D(n)=n/(n+9).
    methods (Test)
        function supportedRetainedNumericControl(t)
            for route={'sing','trt'}
                pd=localSummary(route{1},.2);
                fig=localTable(pd,route{1}); tab=findobj(fig,'Type','uitable');
                t.verifyEqual(tab.Data{1,2},3);
                t.verifyEqual(tab.Data{1,4},5);
                t.verifyEqual(cell2mat(tab.Data(1,6:8)),repmat(5/14,1,3),'AbsTol',1e-12);
                close all force;
            end
        end
        function emptyProjectedRetentionDashesOrdinaryTable(t)
            for route={'sing','trt'}
                pd=localSummary(route{1},.9);
                t.verifyEmpty(pd.relsummary.group.goodids);
                t.verifyTrue(pd.relsummary.group.event.cutoff_extrapolated);
                fig=localTable(pd,route{1}); tab=findobj(fig,'Type','uitable');
                localCheckEmpty(t,tab);
                close all force;
            end
        end
        function emptyFullSampleReferenceExportsRemainAuditable(t)
            localInstallSaveShadow(t);
            for route={'sing','trt'}
                % The D-study search never finds this threshold. Its raw
                % overall fallback still uses all three persons (mean n=5).
                pd=localSummary(route{1},.99999);
                t.verifyEqual(pd.relsummary.group.event.summary_population,'full_sample_reference');
                t.verifyEqual(pd.relsummary.group.event.trlcutoff,-1);
                fig=localTable(pd,route{1}); tab=findobj(fig,'Type','uitable');
                localCheckEmpty(t,tab);
                for ext={'.csv','.xlsx'}
                    file=[tempname ext{1}]; t.addTeardown(@()localDelete(file));
                    localSave(fig,file);
                    if strcmp(ext{1},'.csv'); raw=readcell(file,'Delimiter',',','NumHeaderLines',0);
                    else; raw=readcell(file); end
                    words=string(raw(cellfun(@ischar,raw)));
                    t.verifyTrue(any(contains(lower(words),'no participants retained')));
                    t.verifyTrue(any(contains(lower(words),'full-sample reference')));
                    t.verifyTrue(any(contains(lower(words),'raw criterion reference')));
                    t.verifyTrue(any(contains(words,'0.357142'))); % preserved 5/14 reference
                    r=find(strcmp(raw(:,1),'Label'),1); t.assertNotEmpty(r);
                    hdr=string(raw(r,:)); body=raw(r+1,:);
                    t.verifyEqual(body{hdr=="n_Included"},0);
                    t.assertTrue(any(hdr=="n_Excluded"));
                    t.verifyEqual(body{hdr=="n_Excluded"},3);
                    for col=[4 6 7 8]; t.verifyEqual(body{col},'---'); end
                    t.verifyEqual(body{3},0); % actual criterion cutpoint
                end
                close all force;
            end
        end
        function extrapolationWarningDistinguishesDStudyCutoff(t)
            for route={'sing','trt'}
                pd=localRaw(route{1}); prefs=localPrefs();
                prefs.view.depvalue=.9; prefs.view.relvalue=.9;
                prefs.view.tablecriterion=0; prefs.view.plotcriterion=0;
                prefs.view.noccmode=2;
                psyrat_criterionfigures('psyrat_prefs',prefs,'psyrat_data',pd,'analysis',route{1});
                f=findall(0,'Type','figure','Name','Extrapolation beyond data');
                t.assertNumElements(f,1);
                textObjects=findall(f,'Type','text'); words=strings(0,1);
                for k=1:numel(textObjects)
                    lines=string(textObjects(k).String); words=[words;lines(:)]; %#ok<AGROW>
                end
                msg=strjoin(words,newline);
                t.verifyTrue(contains(msg,'D-study'));
                t.verifyTrue(contains(lower(msg),'criterion'));
                t.verifyTrue(contains(lower(msg),'projection'));
                if strcmp(route{1},'trt')
                    t.verifyTrue(contains(msg,'retained participants'));
                end
                close all force;
            end
        end
        function emptyPopulationPlotExplicitlyNamesModelProjection(t)
            for route={'sing','trt'}
                pd=localSummary(route{1},.9); prefs=localPrefs();
                prefs.view.tablecriterion=0; prefs.view.plotcriterion=1;
                psyrat_criterionfigures('psyrat_prefs',prefs,'psyrat_data',pd,'analysis',route{1});
                f=findobj(0,'Type','figure','Tag','psyrat_output');
                t.assertNumElements(f,1);
                t.verifyTrue(contains(lower(string(f.Name)),'model projection'));
                line=findobj(f,'Type','line');
                t.verifyTrue(any(arrayfun(@(h) any(isfinite(h.YData)),line)));
                close all force;
            end
        end
        function unsupportedFamilyControl(t)
            pd=localRaw('trt'); pd.rel.analysis='trt_splits';
            t.verifyError(@()psyrat_criterionfigures('psyrat_prefs',localPrefs(),...
                'psyrat_data',pd,'analysis','trt'),'criterion:unsupportedDesign');
        end
        function legacyRecordedCountDoesNotAssertFullRoster(t)
            localInstallSaveShadow(t);
            for route={'sing','trt'}
                pd=localSummary(route{1},.99999);
                pd.relsummary.group.event=rmfield(pd.relsummary.group.event,'summary_population');
                fig=localTable(pd,route{1});
                file=[tempname '.csv']; t.addTeardown(@()localDelete(file));
                localSave(fig,file); text=fileread(file);
                t.verifyTrue(contains(text,'legacy_recorded_count_reference'));
                t.verifyTrue(contains(text,'original summary population unavailable'));
                t.verifyFalse(contains(text,'full-sample reference'));
                t.verifyTrue(contains(text,'0.357142'));
                t.verifyTrue(contains(text,'Mu is the fitted model mean'));
                close all force;
            end
        end
        function fittedMuIsQualifiedWithRetainedRows(t)
            pd=localSummary('sing',.2); fig=localTable(pd,'sing');
            controls=findall(fig,'Style','text'); words=strings(0,1);
            for k=1:numel(controls)
                lines=string(controls(k).String); words=[words;lines(:)]; %#ok<AGROW>
            end
            t.verifyTrue(any(contains(words,'Mu is the fitted model mean')));
        end
    end
end

function localCheckEmpty(t,tab)
t.verifyEqual(tab.Data{1,2},0);
t.verifyEqual(tab.Data{1,3},0);
for col=[4 6 7 8]; t.verifyEqual(tab.Data{1,col},'---'); end
names=string(tab.ColumnName);
t.assertTrue(any(names=="n Excluded"));
t.verifyEqual(tab.Data{1,names=="n Excluded"},3);
end

function pd=localRaw(route)
pd=PsyRATTestDataFactory.makeRawSingRelData();
pd.rel.out.sig_u(:)=1; pd.rel.out.sig_e(:)=3;
pd.rel.out.sig_trl=zeros(4,1); pd.rel.out.mu=zeros(4,1);
if strcmp(route,'trt')
    pd.rel.analysis='trt'; pd.rel.time={'T1';'T2'};
    a=pd.rel.data; a.time=repmat({'T1'},height(a),1);
    b=pd.rel.data; b.time=repmat({'T2'},height(b),1);
    pd.rel.data=[a;b];
    for field={'sig_id','sig_occ','sig_trlxid','sig_occxid','sig_trlxocc','sig_err'}
        pd.rel.out.(field{1})=zeros(4,1);
    end
    pd.rel.out.sig_id(:)=1; pd.rel.out.sig_err(:)=3;
end
pd.ver='test'; pd.proc=struct('interactive',false);
pd.rel.filename='criterion_retention_projection'; pd.rel.nchains=2; pd.rel.niter=4;
end

function pd=localSummary(route,target)
pd=localRaw(route);
args={'psyrat_data',pd,'analysis',route,'meascutoff',2,'gcoeff',1,'CI',.95,'interactive',false};
if strcmp(route,'sing'); args=[args,{'depcutoff',target,'depcentmeas',1}];
else; args=[args,{'relcutoff',target,'relcentmeas',1,'reltype',1}]; end
pd=psyrat_relsummary(args{:});
end

function prefs=localPrefs()
prefs=psyrat_defaults(); prefs.view.criterioncutoff=0;
prefs.view.plotcriterion=0; prefs.view.tablecriterion=1; prefs.view.ntrials=6;
prefs.view.depcentmeas=1; prefs.view.relcentmeas=1; prefs.view.reltype=1;
end

function fig=localTable(pd,route)
psyrat_criterionfigures('psyrat_prefs',localPrefs(),'psyrat_data',pd,'analysis',route);
fig=findobj(0,'Type','figure','Tag','psyrat_output');
end

function localInstallSaveShadow(t)
folder=tempname; mkdir(folder); t.addTeardown(@()rmdir(folder,'s'));
fid=fopen(fullfile(folder,'uiputfile.m'),'w');
fprintf(fid,'%s\n','function [filename,pathname,filterindex]=uiputfile(varargin)',...
    'global PSYRAT_CRITERION_SAVE','[p,n,e]=fileparts(PSYRAT_CRITERION_SAVE);',...
    'filename=[n e]; pathname=[p filesep]; filterindex=1;','end');
fclose(fid); addpath(folder,'-begin'); rehash;
end

function localSave(fig,file)
global PSYRAT_CRITERION_SAVE
PSYRAT_CRITERION_SAVE=file;
button=findobj(fig,'Style','pushbutton','String','Save Table');
cb=button.Callback; feval(cb{1},button,[],cb{2:end});
end

function localDelete(file)
if isfile(file); delete(file); end
end
