classdef TestRetentionGroupEdges < PsyRATTestBase
    %G84/G88 group edges. Expected sets/counts are stated directly, and the
    %planted model has D(n)=n/(n+9): no production formula is an oracle.
    %The data factory supplies only input shape and six observed rows/cell.
    methods (Test)
        function completeTRTRosterControl(tc)
            pd=localFixture('trt',true);
            [out,err]=localSummary(pd,'trt',.2);
            tc.verifyEqual(err.nogooddata,0);
            tc.verifyEqual(err.trlmax,0);
            for g=1:2
                group=out.relsummary.group(g);
                tc.verifyEqual(sort(group.observedids),{'s1';'s2';'s3'});
                tc.verifyEqual(sort(group.goodids),{'s1';'s2';'s3'});
                tc.verifyEmpty(group.badids);
                for e=1:2
                    tc.verifyEqual(group.event(e).trlcutoff,3);
                    tc.verifyEqual(group.event(e).relcutoff.m,1/4,'AbsTol',1e-12);
                    tc.verifyEqual(group.event(e).observed_trlinfo.max,6);
                    tc.verifyEqual(group.event(e).goodn,3);
                end
            end
            localVerifyNamedMapping(tc,out);
        end

        function trtMissingFirstAndLastUseObservedGroupUnion(tc)
            for missingEvent={'cor','err'}
                pd=localMissing(missingEvent{1});
                [out,err]=localSummary(pd,'trt',.2);
                localVerifyMissingPartition(tc,out,err);
                localVerifyNamedMapping(tc,out);
                g=out.relsummary.group(1);
                missing=find(strcmp({g.event.name},missingEvent{1}));
                present=find(~strcmp({g.event.name},missingEvent{1}));
                tc.verifyEqual(sort(g.event(missing).eventgoodids),{'s1';'s2'});
                tc.verifyEqual(sort(g.event(present).eventgoodids),{'s1';'s2';'s3'});
                %Missing is not itself an observed below-cutoff event row.
                tc.verifyEmpty(g.event(missing).eventbadids);
            end
        end

        function trtReorderingPreservesPartitionAndNamedComponentMapping(tc)
            for missingEvent={'cor','err'}
                pd=localMissing(missingEvent{1});
                [a,ea]=localSummary(pd,'trt',.2);
                pd.rel.events=flipud(pd.rel.events);
                pd.rel.data=pd.rel.data(end:-1:1,:);
                %The fitted labels retain Group_;_Event identity. Distinct
                %planted mu values below prove correct named reattachment.
                [b,eb]=localSummary(pd,'trt',.2);
                localVerifyMissingPartition(tc,a,ea);
                localVerifyMissingPartition(tc,b,eb);
                localVerifyNamedMapping(tc,b);
                for g=1:2
                    tc.verifyEqual(sort(a.relsummary.group(g).goodids),sort(b.relsummary.group(g).goodids));
                    tc.verifyEqual(sort(a.relsummary.group(g).badids),sort(b.relsummary.group(g).badids));
                    tc.verifyEqual(sort(a.relsummary.group(g).observedids),sort(b.relsummary.group(g).observedids));
                end
            end
        end

        function singGroupSupportUsesSixNotPooledTwelve(tc)
            localVerifyGroupSupport(tc,'sing',.45,8,true);
        end

        function trtGroupSupportUsesSixNotPooledTwelve(tc)
            localVerifyGroupSupport(tc,'trt',.45,8,true);
        end

        function singGroupSupportedCutoffControl(tc)
            localVerifyGroupSupport(tc,'sing',.35,5,false);
        end

        function trtGroupSupportedCutoffControl(tc)
            localVerifyGroupSupport(tc,'trt',.35,5,false);
        end
    end
end

function pd=localMissing(eventName)
pd=localFixture('trt',true);
t=pd.rel.data;
drop=strcmp(t.group,'GroupA') & strcmp(t.event,eventName) & strcmp(t.id,'s3');
pd.rel.data=t(~drop,:);
end

function localVerifyMissingPartition(tc,out,err)
tc.verifyEqual(err.nogooddata,0);
tc.verifyEqual(err.trlmax,0);
g=out.relsummary.group(1);
tc.verifyEqual(sort(g.observedids),{'s1';'s2';'s3'});
tc.verifyEqual(sort(g.goodids),{'s1';'s2'});
tc.verifyEqual(g.badids,{'s3'});
tc.verifyEmpty(intersect(g.goodids,g.badids));
tc.verifyEqual(sort(union(g.goodids,g.badids)),{'s1';'s2';'s3'});
tc.verifyEqual(sort(out.relsummary.group(2).goodids),{'s1';'s2';'s3'});
tc.verifyEmpty(out.relsummary.group(2).badids);
tbl=psyrat_trt_reloverallt('psyrat_data',out,'gui',0);
tc.verifyEqual(cell2mat(tbl.n_Included),[2;2;3;3]);
tc.verifyEqual(cell2mat(tbl.n_Excluded),[1;1;0;0]);
for g=1:2
    for e=1:2
        tc.verifyEqual(out.relsummary.group(g).event(e).goodn,g+1);
    end
end
end

function localVerifyNamedMapping(tc,out)
eventNames={'cor','err'};
for g=1:2
    for e=1:2
        name=out.relsummary.group(g).event(e).name;
        originalEvent=find(strcmp(eventNames,name));
        tc.assertNotEmpty(originalEvent);
        tc.verifyEqual(out.relsummary.data.g(g).e(e).mu.raw,...
            repmat(10*g+originalEvent,4,1));
    end
end
end

function localVerifyGroupSupport(tc,family,target,expectedCutoff,extrapolated)
pd=localFixture(family,false);
%An independent input census exposes the dangerous pooling: each source
%group has six rows per ID (and occasion), while reused IDs pool to twelve.
keys={'id'};
if strcmp(family,'trt'); keys={'id','time'}; end
pooled=varfun(@length,pd.rel.data,'GroupingVariables',keys);
tc.assertEqual(max(pooled.GroupCount),12);
for name={'GroupA','GroupB'}
    rows=pd.rel.data(strcmp(pd.rel.data.group,name{1}),:);
    perGroup=varfun(@length,rows,'GroupingVariables',keys);
    tc.assertEqual(perGroup.GroupCount,repmat(6,height(perGroup),1));
end
tc.assertLessThan((expectedCutoff-1)/(expectedCutoff+8),target);
tc.assertGreaterThan(expectedCutoff/(expectedCutoff+9),target);
[out,err]=localSummary(pd,family,target);
tc.verifyEqual(err.nogooddata,0);
tc.verifyEqual(err.trlmax,double(extrapolated));
for g=1:2
    group=out.relsummary.group(g); ev=group.event(1);
    tc.verifyEqual(sort(group.observedids),{'s1';'s2';'s3'});
    tc.verifyEqual(ev.observed_trlinfo.max,6);
    tc.verifyEqual(ev.observed_trlinfo.mean,6);
    tc.verifyEqual(ev.trlcutoff,expectedCutoff);
    tc.verifyEqual(ev.cutoff_extrapolated,extrapolated);
    if strcmp(family,'trt'); coefficient=ev.relcutoff;
    else; coefficient=ev.rel; end
    expected=expectedCutoff/(expectedCutoff+9);
    tc.verifyEqual([coefficient.ll coefficient.m coefficient.ul],repmat(expected,1,3),'AbsTol',1e-12);
    if extrapolated
        tc.verifyEmpty(group.goodids);
        tc.verifyEqual(sort(group.badids),{'s1';'s2';'s3'});
        tc.verifyEqual(ev.goodn,0);
    else
        tc.verifyEqual(sort(group.goodids),{'s1';'s2';'s3'});
        tc.verifyEmpty(group.badids);
        tc.verifyEqual(ev.goodn,3);
    end
end
if strcmp(family,'trt')
    tbl=psyrat_trt_relcutofft('psyrat_data',out,'gui',0);
else
    tbl=psyrat_depcutofft('psyrat_data',out,'gui',0);
end
tc.verifyEqual(contains(string(tbl{:,3}),'projection beyond observed data'),repmat(extrapolated,2,1));
end

function pd=localFixture(family,withEvents)
pd=PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
if withEvents
    pd.rel.data=vertcat(pd.rel.data{:}); nevents=2;
else
    pd.rel.data=pd.rel.data{1}; pd.rel.data.event=[];
    pd.rel.events='none'; nevents=1;
end
ngroups=2; ncolumns=ngroups*nevents;
out=struct('labels',{cell(ncolumns,1)},'mu',zeros(4,ncolumns));
for g=1:ngroups
    for e=1:nevents
        col=(g-1)*nevents+e;
        out.mu(:,col)=10*g+e;
        if withEvents
            %The TRT case-4 contract is Group_;_Event, unlike plain IC.
            out.labels{col}=[pd.rel.groups{g} '_;_' pd.rel.events{e}];
        else
            out.labels{col}=pd.rel.groups{g};
        end
    end
end
if strcmp(family,'trt')
    pd.rel.analysis='trt'; pd.rel.time={'T1';'T2'};
    a=pd.rel.data; a.time=repmat({'T1'},height(a),1);
    b=pd.rel.data; b.time=repmat({'T2'},height(b),1);
    pd.rel.data=[a;b];
    for f={'sig_id','sig_occ','sig_trl','sig_trlxid','sig_occxid','sig_trlxocc','sig_err'}
        out.(f{1})=zeros(4,ncolumns);
    end
    out.sig_id(:)=1; out.sig_err(:)=3;
else
    out.sig_u=ones(4,ncolumns); out.sig_e=3*ones(4,ncolumns);
    out.sig_trl=zeros(4,ncolumns);
end
pd.rel.out=out; pd.ver='test'; pd.proc=struct('interactive',false);
pd.rel.filename='retention_group_edge_fixture'; pd.rel.nchains=2; pd.rel.niter=4;
end

function [out,err]=localSummary(pd,family,target)
args={'psyrat_data',pd,'analysis',family,'meascutoff',2,'gcoeff',1,'CI',.95,'interactive',false};
if strcmp(family,'trt'); args=[args,{'reltype',1,'relcutoff',target,'relcentmeas',1}];
else; args=[args,{'depcutoff',target,'depcentmeas',1}]; end
[out,err]=psyrat_relsummary(args{:});
end
