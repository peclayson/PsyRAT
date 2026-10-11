classdef TestBootstrapValidity < PsyRATTestBase
    % PBA-05: real-fit failure/control comparisons; no fitted-point padding.
    methods (Test)
        function warningSymmetryAndReplicateRecord(t)
            d = localData();
            c = localShadow(t,0,0); %#ok<NASGU>
            clean = psyrat_native_fitlme(struct('key','sing'),d,1979,4,struct('bootstrap_reps',10));
            localConfigure(2,0);
            warned = psyrat_native_fitlme(struct('key','sing'),d,1979,4,struct('bootstrap_reps',10));
            t.verifyEqual(warned.meta.bootstrap_ok,9);
            t.verifyEqual(warned.draws.mu,clean.draws.mu(2:end));
            t.verifyEqual(warned.draws.sig_u,clean.draws.sig_u(2:end));
            t.verifyEqual(warned.meta.bootstrap_success_mask,[false;true(9,1)]);
            t.verifyEqual(warned.meta.bootstrap_records(1).warning_id,'PsyRATTest:OptimizerNonConvergence');
            t.verifyFalse(warned.meta.bootstrap_records(1).valid);
            localConfigure(1,0);
            original = psyrat_native_fitlme(struct('key','sing'),d,1979,4,struct('bootstrap_reps',10));
            t.verifyFalse(original.meta.original.valid);
            t.verifyEqual(original.conv{2,3},Inf);
        end

        function independentStrataKeepEverySuccess(t)
            d=localData();
            data=struct('id_G1',d.id,'trl_G1',d.trl,'meas_G1',d.meas, ...
                'id_G2',d.id,'trl_G2',d.trl,'meas_G2',d.meas+.1);
            c=localShadow(t,0,0); %#ok<NASGU>
            spec=struct('key','sing_facet','suffix','G','nstrata',2);
            clean=psyrat_native_fitlme(spec,data,1979,4,struct('bootstrap_reps',10));
            % Failing refit does not consume fit RNG; simulation stream is identical.
            localConfigure(2,0);
            bad=psyrat_native_fitlme(spec,data,1979,4,struct('bootstrap_reps',10));
            t.verifyEqual(bad.draws.mu_G1,clean.draws.mu_G1(2:end));
            t.verifyEqual(bad.draws.mu_G2,clean.draws.mu_G2);
            t.verifyEqual([bad.meta.fits.bootstrap_ok],[9 10]);
        end

        function finiteRealAndZeroBoundaryContract(t)
            z=struct('mu',-3,'sig_u',0,'sig_trl',1,'sig_e',0);
            t.verifyTrue(psyrat_fitlme_validity(z,'','').valid);
            for v={NaN,Inf,1i}
                b=z; b.sig_u=v{1};
                t.verifyFalse(psyrat_fitlme_validity(b,'','').valid);
            end
            t.verifyFalse(psyrat_fitlme_validity(z,'stats:OptimizerNonConvergence','limit').valid);
            t.verifyTrue(psyrat_fitlme_validity(z,'stats:Informational','info').valid);
        end

        function raggedLayoutPreservesValuesAndMeans(t)
            x=psyrat_append_draws([], (1:9)');
            t.verifyTrue(isnumeric(x));
            x=psyrat_append_draws(x,(11:20)');
            t.verifyTrue(iscell(x));
            t.verifyEqual(psyrat_draw_column(x,1),(1:9)');
            t.verifyEqual(psyrat_draw_column(x,2),(11:20)');
            t.verifyEqual(mean(psyrat_draw_column(x,1)),5);
            t.verifyEqual(mean(psyrat_draw_column(x,2)),15.5);
            y=psyrat_append_draws((1:10)',(11:20)');
            t.verifySize(y,[10 2]);
        end

        function extractionRejectsNonfiniteAndKeepsZeroBoundary(t)
            d=localData();
            for v={NaN,Inf,-1,0}
                c=psyrat_test_fitlme_injection(0,0,2,v{1});
                f=psyrat_native_fitlme(struct('key','sing'),d,1979,4,struct('bootstrap_reps',10));
                if isequal(v{1},0)
                    t.verifyEqual(f.meta.bootstrap_ok,10);
                    t.verifyEqual(f.draws.sig_u(1),0);
                else
                    t.verifyEqual(f.meta.bootstrap_ok,9);
                    t.verifyFalse(f.meta.bootstrap_success_mask(1));
                    t.verifyEqual(f.meta.bootstrap_records(1).error_id,'PsyRAT:fitlmeInvalidEstimate');
                end
                clear c
            end
        end

        function caseFourSaveSummaryAndDisclosure(t)
            d=localData();
            id=arrayfun(@(x)sprintf('S%02d',x),[d.id;d.id;d.id;d.id],'UniformOutput',false);
            meas=[d.meas;d.meas+.1;d.meas+.2;d.meas+.3];
            group=[repmat({'G1'},80,1);repmat({'G2'},80,1)]; group=[group;group];
            event=[repmat({'A'},160,1);repmat({'B'},160,1)];
            c=localShadow(t,0,2); %#ok<NASGU>
            rel=psyrat_computevarcomp('data',table(id,meas,group,event),'engine',2, ...
                'bootstrap_reps',10,'chains',4,'sampling',10,'warmup',10, ...
                'seed',1979,'verbose',1,'showgui',0);
            t.verifyEqual(cellfun(@numel,rel.out.mu),[9 10 10 10]);
            t.verifyEqual([rel.engine_meta.fits.bootstrap_ok],[9 10 10 10]);
            t.verifyEqual(rel.engine_meta.bootstrap_ok,39);
            t.verifyEqual(rel.engine_meta.bootstrap_reps,40);
            t.verifyEqual(rel.engine_meta.fits(1).bootstrap_records(1).error_id,'PsyRATTest:InjectedRefitFailure');
            pd=struct('rel',rel,'proc',struct('interactive',false));
            [pd,~]=psyrat_relsummary('psyrat_data',pd,'analysis','sing','depcutoff',.1, ...
                'meascutoff',2,'depcentmeas',2,'interactive',false);
            t.verifyFalse(isempty(pd.relsummary));
            lines=psyrat_provenance_lines(pd);
            t.verifyTrue(any(contains(lines,'9 of 10 succeeded')));
            t.verifyTrue(contains(psyrat_engine_notice(rel),'A_;_G1'));
            dirpath=tempname; mkdir(dirpath); t.addTeardown(@()rmdir(dirpath,'s'));
            pd.proc.savepath=dirpath; pd.proc.savename='bootstrap.psyrat';
            save(fullfile(dirpath,'bootstrap.psyrat'),'pd','-mat');
            loaded=load(fullfile(dirpath,'bootstrap.psyrat'),'-mat');
            t.verifyEqual(loaded.pd.rel.out.mu,rel.out.mu);
            prefs=psyrat_defaults();
            prefs.proc.bootstrap_reps=10;
            config=psyrat_write_runconfig(pd,prefs);
            cfg=jsondecode(fileread(config));
            t.verifyEqual(cfg.estimation.bootstrap_reps,10);
            t.verifyEqual(cfg.estimation.bootstrap_requested_total,40);
            t.verifyEqual(cfg.estimation.bootstrap_counts_scope,'all independent fits');
            t.verifyEqual([cfg.estimation.bootstrap_fits.bootstrap_ok],[9 10 10 10]);
            t.verifyEqual(cfg.estimation.bootstrap_fits(1).bootstrap_success_mask,rel.engine_meta.fits(1).bootstrap_success_mask);
            t.verifyWarning(@()psyrat_checktraceplots(rel),'PsyRAT:bootstrapTraceUnavailable');
        end

        function crossedAndWeightedStrataReachSummaries(t)
            %Full case5 and case24 routing with two independent strata.
            %The warning changes validity only; all clean fitted values remain.
            rng(137); [id,time,trl]=ndgrid(1:12,1:3,1:4);
            id=id(:); time=time(:); trl=trl(:); n=numel(id);
            u=2*randn(12,1); o=randn(3,1); v=randn(4,1);
            meas=1+u(id)+o(time)+v(trl)+randn(n,1);
            ids=arrayfun(@(x)sprintf('S%02d',x),[id;id],'UniformOutput',false);
            group=[repmat({'G1'},n,1);repmat({'G2'},n,1)];
            T=table(ids,[meas;meas+.1],[time;time],group, ...
                'VariableNames',{'id','meas','time','group'});
            for weighted=[false true]
                input=T; args={};
                if weighted
                    input.weight=repmat(5+5*trl,2,1); args={'splits',3};
                end
                c=localShadow(t,0,0);
                clean=psyrat_computevarcomp('data',input,'engine',2,'bootstrap_reps',10, ...
                    'chains',4,'sampling',10,'warmup',10,'seed',1979,'verbose',1,'showgui',0,args{:});
                localConfigure(2,0);
                rel=psyrat_computevarcomp('data',input,'engine',2,'bootstrap_reps',10, ...
                    'chains',4,'sampling',10,'warmup',10,'seed',1979,'verbose',1,'showgui',0,args{:});
                t.verifyEqual(cellfun(@numel,rel.out.mu),[9 10]);
                t.verifyEqual(rel.out.mu{1},clean.out.mu(2:end,1));
                t.verifyEqual(rel.out.mu{2},clean.out.mu(:,2));
                t.verifyEqual([rel.engine_meta.fits.bootstrap_ok],[9 10]);
                if weighted
                    summary=psyrat_splits_summary(rel);
                    t.verifyEqual(numel(summary.strata),2);
                else
                    pd=struct('rel',rel,'proc',struct('interactive',false));
                    [pd,~]=psyrat_relsummary('psyrat_data',pd,'analysis','trt','relcutoff',.1, ...
                        'meascutoff',2,'relcentmeas',2,'gcoeff',1,'reltype',1,'interactive',false);
                    tbl=psyrat_trt_reloverallt('psyrat_data',pd,'gui',0);
                    t.verifyFalse(isempty(tbl));
                end
                clear c
            end
        end
    end
end

function d=localData()
rng(1979); id=repelem((1:10)',8); trl=repmat((1:8)',10,1);
u=2*randn(10,1); v=randn(8,1);
d=struct('id',id,'trl',trl,'meas',1+u(id)+v(trl)+randn(80,1));
end
function c=localShadow(~,w,f)
c=psyrat_test_fitlme_injection(w,f);
end
function localConfigure(w,f)
s=getappdata(0,'psyrat_test_fitlme_injection');
s.calls=0; s.warnAt=w; s.failAt=f;
setappdata(0,'psyrat_test_fitlme_injection',s);
end
