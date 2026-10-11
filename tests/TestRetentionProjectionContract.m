classdef TestRetentionProjectionContract < PsyRATTestBase
    % G84/G88: expectations are explicit set identities and D(n)=n/(n+9),
    % independently specified before production repair. No estimator runs.
    methods (Test)
        function completeRosterControl(tc)
            [out,err] = localSummary(localGroups(),.5,'sing');
            tc.verifyEqual(err.nogooddata,0);
            for g=1:2
                tc.verifyEqual(sort(out.relsummary.group(g).goodids),{'s1';'s2';'s3'});
                tc.verifyEmpty(out.relsummary.group(g).badids);
            end
        end

        function missingFirstAndLastPreserveObservedPartition(tc)
            for absent=1:2
                pd=localGroups(); t=pd.rel.data{absent};
                pd.rel.data{absent}=t(~(strcmp(t.group,'GroupA') & strcmp(t.id,'s3')),:);
                [out,err]=localSummary(pd,.5,'sing');
                tc.verifyEqual(err.nogooddata,0);
                g=out.relsummary.group(1);
                tc.verifyEqual(sort(g.goodids),{'s1';'s2'});
                tc.verifyEqual(g.badids,{'s3'});
                tc.verifyEqual(sort(out.relsummary.group(2).goodids),{'s1';'s2';'s3'});
                tbl=psyrat_depoverallt('psyrat_data',out,'gui',0);
                tc.verifyEqual(cell2mat(tbl.n_Included(1:2)),[2;2]);
                tc.verifyEqual(cell2mat(tbl.n_Excluded(1:2)),[1;1]);
            end
        end

        function eventReorderingCannotChangeExcludedIDs(tc)
            pd=localGroups(); t=pd.rel.data{1};
            pd.rel.data{1}=t(~(strcmp(t.group,'GroupA') & strcmp(t.id,'s3')),:);
            a=localSummary(pd,.5,'sing');
            pd.rel.events=pd.rel.events([2 1]); pd.rel.data=pd.rel.data([2 1]);
            b=localSummary(pd,.5,'sing');
            for g=1:2
                tc.verifyEqual(sort(a.relsummary.group(g).goodids),sort(b.relsummary.group(g).goodids));
                tc.verifyEqual(sort(a.relsummary.group(g).badids),sort(b.relsummary.group(g).badids));
            end
        end

        function emptyGroupDashesRetainedSummariesAndPreservesReference(tc)
            pd=localGroups(); pd.rel.out.sig_e(:,2)=30;
            [out,err]=localSummary(pd,.5,'sing');
            tc.verifyEqual(err.nogooddata,0);
            tc.verifyEmpty(out.relsummary.group(1).goodids);
            tc.verifyEqual(sort(out.relsummary.group(1).badids),{'s1';'s2';'s3'});
            tbl=psyrat_depoverallt('psyrat_data',out,'gui',0);
            for row=1:2
                for col=4:9; tc.verifyEqual(tbl{row,col},{'---'}); end
            end
            tc.verifyEqual(cell2mat(tbl.n_Included(1:2)),[0;0]);
            tc.verifyEqual(cell2mat(tbl.n_Excluded(1:2)),[3;3]);
            tc.verifyNotEqual(tbl{3,4},{'---'});
            tc.verifyGreaterThan(out.relsummary.group(1).event(1).dep.m,0);
            tc.assertTrue(isfield(out.relsummary.group(1).event,'summary_population'));
            tc.verifyEqual(out.relsummary.group(1).event(1).summary_population,'full_sample_reference');
        end

        function singleProjectionRetainsValueAndMarksOriginalSupport(tc)
            localProjectionCheck(tc,'sing',.9,true);
        end

        function singleSupportedCutoffControl(tc)
            localProjectionCheck(tc,'sing',.2,false);
        end

        function trtProjectionRetainsValueAndMarksOriginalSupport(tc)
            localProjectionCheck(tc,'trt',.9,true);
        end

        function trtSupportedCutoffControl(tc)
            localProjectionCheck(tc,'trt',.2,false);
        end

        function actualIDExportsPreserveGroupIdentity(tc)
            pd=localGroups(); t=pd.rel.data{1};
            pd.rel.data{1}=t(~(strcmp(t.group,'GroupA') & strcmp(t.id,'s3')),:);
            out=localSummary(pd,.5,'sing');
            localInstallSaveShadow(tc);
            psyrat_depoverallt('psyrat_data',out,'gui',1);
            fig=findobj(0,'Type','figure','Tag','psyrat_output');
            for ext={'.csv','.xlsx'}
                file=[tempname ext{1}]; tc.addTeardown(@()localDelete(file));
                localFire(tc,fig(1),file,'Save IDs');
                raw=localReadExport(file); r=find(strcmp(raw(:,1),'Good IDs'),1);
                tc.assertNotEmpty(r);
                tc.assertGreaterThanOrEqual(size(raw,2),4,'Export must retain group identity beside each ID');
                tc.verifyEqual(raw(r,1:4),{'Good IDs','Bad IDs','Good Group','Bad Group'});
                body=raw(r+1:end,:);
                pairs=strcat(string(body(:,3)),"/",string(body(:,1)));
                tc.verifyEqual(sort(pairs),sort(["GroupA/s1";"GroupA/s2";"GroupB/s1";"GroupB/s2";"GroupB/s3"]));
                tc.verifyEqual(body(1,[2 4]),{'s3','GroupA'});
            end
        end

        function actualProjectionAndEmptyTableExports(tc)
            localInstallSaveShadow(tc);
            for family={'sing','trt'}
                pd=localSimple(family{1}); out=localSummary(pd,.9,family{1});
                if strcmp(family{1},'sing'); builders={@psyrat_depcutofft,@psyrat_depoverallt};
                else; builders={@psyrat_trt_relcutofft,@psyrat_trt_reloverallt}; end
                for k=1:2
                    close all force; builders{k}('psyrat_data',out,'gui',1);
                    fig=findobj(0,'Type','figure','Tag','psyrat_output');
                    for ext={'.csv','.xlsx'}
                        file=[tempname ext{1}]; tc.addTeardown(@()localDelete(file));
                        localFire(tc,fig(1),file,'Save Table');
                        raw=localReadExport(file); words=string(raw(cellfun(@ischar,raw)));
                        if k==1
                            tc.verifyTrue(any(contains(words,'projection beyond observed data')));
                        else
                            tc.verifyGreaterThanOrEqual(sum(words=="---"),6);
                        end
                    end
                end
            end
        end

        function literalNoneIDAndUngroupedExportControl(tc)
            localInstallSaveShadow(tc);
            for family={'sing','trt'}
                pd=localSimple(family{1});
                pd.rel.data.id(strcmp(pd.rel.data.id,'s3'))={'none'};
                out=localSummary(pd,.2,family{1});
                tc.verifyEqual(sort(out.relsummary.group(1).goodids),{'none';'s1';'s2'});
                close all force;
                if strcmp(family{1},'sing'); psyrat_depoverallt('psyrat_data',out,'gui',1);
                else; psyrat_trt_reloverallt('psyrat_data',out,'gui',1); end
                fig=findobj(0,'Type','figure','Tag','psyrat_output');
                for ext={'.csv','.xlsx'}
                    file=[tempname ext{1}]; tc.addTeardown(@()localDelete(file));
                    localFire(tc,fig(1),file,'Save IDs'); raw=localReadExport(file);
                    r=find(strcmp(raw(:,1),'Good IDs'),1); tc.assertNotEmpty(r);
                    tc.assertGreaterThanOrEqual(size(raw,2),4,'Export must retain group identity beside each ID');
                    tc.verifyEqual(raw(r,1:4),{'Good IDs','Bad IDs','Good Group','Bad Group'});
                    body=raw(r+1:end,:);
                    tc.verifyEqual(sort(string(body(:,1))),["none";"s1";"s2"]);
                    tc.verifyEqual(string(body(:,3)),repmat("none",3,1));
                end
            end
        end
    end
end

function localProjectionCheck(tc,family,target,extrapolated)
pd=localSimple(family); [out,err]=localSummary(pd,target,family);
ev=out.relsummary.group(1).event(1);
tc.verifyEqual(err.trlmax,double(extrapolated));
if strcmp(family,'sing'); coeff=ev.rel; builder=@psyrat_depcutofft;
else; coeff=ev.relcutoff; builder=@psyrat_trt_relcutofft; end
tc.verifyEqual(coeff.m,ev.trlcutoff/(ev.trlcutoff+9),'AbsTol',1e-12);
tc.verifyEqual(coeff.ll,coeff.m,'AbsTol',1e-12);
tc.verifyEqual(coeff.ul,coeff.m,'AbsTol',1e-12);
tc.assertTrue(isfield(ev,'observed_trlinfo'));
tc.verifyEqual(ev.observed_trlinfo.max,6);
tc.verifyEqual(ev.cutoff_extrapolated,extrapolated);
tbl=builder('psyrat_data',out,'gui',0);
tc.verifyEqual(contains(string(tbl{1,3}),'projection beyond observed data'),extrapolated);
if extrapolated
    tc.verifyGreaterThan(ev.trlcutoff,6); tc.verifyEqual(ev.goodn,0);
else
    tc.verifyEqual(ev.trlcutoff,3); tc.verifyEqual(ev.goodn,3);
end
end

function pd=localGroups()
pd=PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
pd.rel.out.sig_u(:)=2; pd.rel.out.sig_e(:)=1; pd.rel.out.sig_trl(:)=.2;
pd=localHeaders(pd);
end

function pd=localSimple(family)
pd=PsyRATTestDataFactory.makeRawSingRelData();
pd.rel.out.sig_u(:)=1; pd.rel.out.sig_e(:)=3;
pd.rel.out.sig_trl=zeros(4,1);
if strcmp(family,'trt')
    pd.rel.analysis='trt'; pd.rel.time={'T1';'T2'};
    a=pd.rel.data; a.time=repmat({'T1'},height(a),1);
    b=pd.rel.data; b.time=repmat({'T2'},height(b),1);
    pd.rel.data=[a;b];
    for f={'sig_id','sig_occ','sig_trlxid','sig_occxid','sig_trlxocc','sig_err'}
        pd.rel.out.(f{1})=zeros(4,1);
    end
    pd.rel.out.sig_id(:)=1; pd.rel.out.sig_err(:)=3;
end
pd=localHeaders(pd);
end

function pd=localHeaders(pd)
pd.ver='test'; pd.proc=struct('interactive',false);
pd.rel.filename='retention_projection_fixture'; pd.rel.nchains=2; pd.rel.niter=4;
end

function [out,err]=localSummary(pd,target,family)
args={'psyrat_data',pd,'analysis',family,'meascutoff',2,'gcoeff',1,'CI',.95,'interactive',false};
if strcmp(family,'trt'); args=[args,{'reltype',1,'relcutoff',target,'relcentmeas',1}];
else; args=[args,{'depcutoff',target,'depcentmeas',1}]; end
[out,err]=psyrat_relsummary(args{:});
end

function localInstallSaveShadow(tc)
folder=tempname; mkdir(folder); tc.addTeardown(@()rmdir(folder,'s'));
fid=fopen(fullfile(folder,'uiputfile.m'),'w');
fprintf(fid,'%s\n','function [filename,pathname,filterindex]=uiputfile(varargin)',...
    'global PSYRAT_FAKE_SAVE','[p,n,e]=fileparts(PSYRAT_FAKE_SAVE);',...
    'filename=[n e]; pathname=[p filesep]; filterindex=1;','end');
fclose(fid); addpath(folder,'-begin'); rehash;
end

function localFire(tc,fig,file,label)
global PSYRAT_FAKE_SAVE
PSYRAT_FAKE_SAVE=file; found=false;
for b=findall(fig,'Style','pushbutton')'
    if contains(string(b.String),label)
        cb=b.Callback; feval(cb{1},b,[],cb{2:end}); found=true; break;
    end
end
tc.assertTrue(found);
end

function localDelete(file)
if isfile(file); delete(file); end
end

function raw=localReadExport(file)
%All-good IDs leave one CSV column empty throughout the body. Do not allow
%delimiter inference to reinterpret the comma-separated header as one field.
if endsWith(file,'.csv'); raw=readcell(file,'Delimiter',',','NumHeaderLines',0);
else; raw=readcell(file); end
end
