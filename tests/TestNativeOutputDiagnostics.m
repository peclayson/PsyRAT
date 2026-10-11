classdef TestNativeOutputDiagnostics < PsyRATTestBase
    methods (Test)
        function preservesDrawsAndChainCoordinates(t)
            [fit,x]=localFit(t,'different_scale');
            fit.draws.ind_sd=repmat(x(:),1,2); fit.draws.Intercept_sigma=zeros(numel(x),1);
            fit.draws.b=zeros(numel(x),0);
            before=fit.draws;
            out=psyrat_diagnose_native_outputs(fit,'sserr',localData(),size(x,1),size(x,2));
            t.verifyEqual(out.draws,before);
            row=find(strcmp(out.conv(:,1),'returned.ind_sd[1]'));
            d=psyrat_rank_diagnostics(x);
            t.verifyEqual(out.conv{row,3},d.rhat,'AbsTol',1e-12);
            t.verifyGreaterThan(out.conv{row,3},1.1);
            row=find(strcmp(out.conv(:,1),'returned.residual_variance_at_dimension_zero[1]'));
            d=psyrat_rank_diagnostics(exp(2*x));
            t.verifyEqual(out.conv{row,6},d.mcse_mean,'RelTol',1e-12);
            t.verifyTrue(any(contains(out.meta.diagnostics.absent_fields,'b: zero')));
        end
        function enabledAndDisabledIdenticalValues(t)
            [fit,x]=localFit(t,'iid','diff_sserr'); n=numel(x);
            fit.draws.rescor=zeros(n,1); fit.draws.wp_cov_ss=zeros(n,2);
            fit.draws.err_varcov=zeros(n,2,2);
            fit.draws.err_varcov(:,1,1)=exp(x(:));
            fit.draws.err_varcov(:,2,2)=exp(x(:));
            off=psyrat_diagnose_native_outputs(fit,'diff_sserr',localData(0),size(x,1),size(x,2));
            on=psyrat_diagnose_native_outputs(fit,'diff_sserr',localData(1),size(x,1),size(x,2));
            t.verifyEqual(localGate(off).out.conv.converged,1);
            t.verifyEqual(localGate(on).out.conv.converged,0);
            t.verifyFalse(any(contains(off.meta.diagnostics.permitted_constants,'[1,1]')));
        end
        function structuralMatricesAndInvalidControls(t)
            [fit,x]=localFit(t,'iid'); n=numel(x);
            fit.draws.chol_corrmat=zeros(n,2,2);
            fit.draws.chol_corrmat(:,1,1)=1;
            fit.draws.chol_corrmat(:,2,1)=tanh(x(:));
            fit.draws.chol_corrmat(:,2,2)=sqrt(1-tanh(x(:)).^2);
            good=psyrat_diagnose_native_outputs(fit,'sserr',localData(),size(x,1),size(x,2));
            t.verifyEqual(localGate(good).out.conv.converged,1);
            fit.draws.chol_corrmat(:,2,2)=0;
            bad=psyrat_diagnose_native_outputs(fit,'sserr',localData(),size(x,1),size(x,2));
            t.verifyEqual(localGate(bad).out.conv.converged,0);
            fit.draws=struct('cor_ID_Cell',ones(n,4,4));
            for j=1:4; fit.draws.cor_ID_Cell(:,j,j)=1+eps*sin(x(:)); end
            out=psyrat_diagnose_native_outputs(fit,'dod',localData(),size(x,1),size(x,2));
            t.verifyEqual(out.draws,fit.draws);
            t.verifyNumElements(out.meta.diagnostics.permitted_constants,4);
            fit.draws.cor_ID_Cell(1,1,1)=NaN;
            out=psyrat_diagnose_native_outputs(fit,'dod',localData(),size(x,1),size(x,2));
            t.verifyFalse(any(strcmp(out.meta.diagnostics.permitted_constants,'returned.cor_ID_Cell[1,1]')));
            fit.draws=struct('bad_shape',ones(n-1,2));
            out=psyrat_diagnose_native_outputs(fit,'sserr',localData(),size(x,1),size(x,2));
            t.verifyEqual(localGate(out).out.conv.converged,0);
        end
        function requiredFieldsShapesAndBound(t)
            [fit,x]=localFit(t,'iid'); data=localData();
            fit.draws=rmfield(fit.draws,'ind_sd');
            out=psyrat_diagnose_native_outputs(fit,'sserr',data,size(x,1),size(x,2));
            t.verifyTrue(any(strcmp(out.conv(:,1),'required_shape.ind_sd')));
            [fit,x]=localFit(t,'iid'); fit.draws.ind_sd=zeros(numel(x),0);
            out=psyrat_diagnose_native_outputs(fit,'sserr',data,size(x,1),size(x,2));
            t.verifyTrue(any(strcmp(out.conv(:,1),'required_shape.ind_sd')));
            t.verifyFalse(any(contains(out.meta.diagnostics.absent_fields,'ind_sd')));
            [fit,x]=localFit(t,'iid','dod'); fit.draws.cor_ID_Cell(:,1,1)=2;
            out=psyrat_diagnose_native_outputs(fit,'dod',data,size(x,1),size(x,2));
            t.verifyFalse(any(strcmp(out.meta.diagnostics.permitted_constants,'returned.cor_ID_Cell[1,1]')));
            t.verifyEqual(out.draws,fit.draws);
        end
        function dodViewsRetainSelectedDiagnosticEvidence(t)
            pd=PsyRATTestDataFactory.makeDoDRelData(); pd.rel.nchains=4;
            m=psyrat_diagnostic_metadata({}); m.format='modern_csv'; m.degraded=false;
            c={'name','n_eff','r_hat';'mu',400,1}; c{1,4}=3; c{1,5}=m;
            pd.rel.out.conv.data={c,c};
            vd=psyrat_dod_build_virtual_rel('psyrat_data',pd,'mode','single','groupidx',2,'eventidx',3);
            t.verifyEqual(vd.rel.out.conv.ndivergent,3);
            t.verifyEqual(vd.rel.out.conv.diagnostics{1}.format,'modern_csv');
            t.verifyEqual(vd.rel.out.conv.final_inference,'unassessed');
            t.verifySubstring(psyrat_convergence_notice(vd.rel),'3 divergent');
        end
        function chainInventoryAndSavedDisclosure(t)
            path=[tempname '.csv']; fid=fopen(path,'w');
            fprintf(fid,'# fixture chain\nlp__,mu,cor_id.1.1,cor_id.2.1\n'); fclose(fid);
            t.addTeardown(@()delete(path));
            [names,status]=psyrat_stan_expected_names({path},@(n)~strcmp(n,'lp__'));
            t.verifyEqual(names,{'mu','cor_id[1,1]','cor_id[2,1]'});
            t.verifyEqual(status,'chain_csv_headers');
            rows={'mu',400,1,400,450,.01}; info=psyrat_diagnostic_metadata({});
            info.format='modern_text'; info.precision_significant_digits=2;
            c=psyrat_diagnostic_table(rows,info,'parameters { real mu; }',struct(),names);
            fit=PsyRATNativeFit(struct(),c,'hmc',struct()); rel=localGate(fit);
            t.verifyEqual(rel.out.conv.converged,0);
            t.verifyEqual(rel.out.conv.invalid,{'cor_id[1,1]';'cor_id[2,1]'});
            msg=psyrat_convergence_notice(rel);
            t.verifySubstring(msg,'Early-workflow'); t.verifySubstring(msg,'Final-inference readiness: unassessed');
            t.verifySubstring(msg,'degraded');
            p=psyrat_provenance_lines(struct('rel',rel));
            t.verifyTrue(any(contains(p,'modern_text')));
            file=[tempname '.mat']; save(file,'rel'); t.addTeardown(@()delete(file));
            loaded=load(file); t.verifyEqual(loaded.rel,rel);
            pd=struct('rel',rel,'proc',struct('savepath',tempdir,'savename',[char(java.util.UUID.randomUUID) '.psyrat']));
            cfgfile=psyrat_write_runconfig(pd,psyrat_defaults()); t.addTeardown(@()delete(cfgfile));
            cfg=jsondecode(fileread(cfgfile));
            t.verifyEqual(cfg.diagnostics.final_inference,'unassessed');
            t.verifyEqual(cfg.diagnostics.diagnostics.format,'modern_text');
            t.verifyTrue(cfg.diagnostics.diagnostics.degraded);
            t.verifyEqual(cfg.diagnostics.diagnostics.precision_significant_digits,2);

        end
    end
end
function [fit,x]=localFit(t,name,key)
if nargin<3; key='sserr'; end
path=fullfile(t.projectRoot(),'tests','fixtures','rank_diagnostics',[name '.bin']);
fid=fopen(path,'rb','ieee-le'); x=fread(fid,[2000,4],'double'); fclose(fid);
c={'name','n_eff','r_hat';'moving',400,1};
m=psyrat_diagnostic_metadata({}); m.format='native_rank_fold'; m.degraded=false;
m.expected_names={'moving'}; c{1,7}=m;
%The independent source inventory, not the production schema, builds controls.
inv=jsondecode(fileread(fullfile(t.projectRoot(),'tests','fixtures','native-returned-field-inventory.json')));
spec=inv.keys.(key).fields; data=localData(); fields=fieldnames(spec); draws=struct();
for k=1:numel(fields)
    dims=spec.(fields{k});
    if iscell(dims)
        for j=1:numel(dims)
            if ischar(dims{j})
                if strcmp(dims{j},'K'); dims{j}=data.KDIM; else; dims{j}=data.(dims{j}); end
            end
        end
        dims=cell2mat(dims);
    end
    dims=reshape(dims,1,[]); draws.(fields{k})=repmat(reshape(exp(.1*x(:)),[],1),[1 dims]);
end
fit=PsyRATNativeFit(draws,c,'hmc',struct());
end
function data=localData(flag)
if nargin<1; flag=0; end
data=struct('KDIM',0,'NSUB',2,'C',4,'estimate_rescor',flag);
end
function rel=localGate(fit)
rel=struct('nchains',4,'engine','hmc','out',struct('conv',struct('data',{fit.conv})));
rel=psyrat_checkconv(rel);
end
