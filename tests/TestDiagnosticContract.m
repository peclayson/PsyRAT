classdef TestDiagnosticContract < PsyRATTestBase
    methods (Test)
        function independentSourceProofControls(t)
            root=fullfile(t.projectRoot(),'tests','fixtures','constant_source');
            for group={'expression_cases','zero_scope_cases','alias_scope_cases','parent_order_cases','string_legal_cases'}
                folder=fullfile(root,group{1});
                manifest=jsondecode(fileread(fullfile(folder,'manifest.json')));
                for k=1:numel(manifest.fixtures)
                    f=manifest.fixtures(k); if iscell(f); f=f{1}; end
                    src=fileread(fullfile(folder,f.file));
                    expected=cellstr(string(f.expected_bounded_rule_names(:)));
                    names=cellstr(string(f.probed_names(:)));
                    c=psyrat_constant_contract(src,f.data,names);
                    t.verifyEqual(c.names(:),expected(:),f.file);
                    rows=[names num2cell(nan(numel(names),5)); {'healthy',400,1,400,400,.01}];
                    tbl=psyrat_diagnostic_table(rows,psyrat_diagnostic_metadata({}),src,f.data);
                    rel=struct('nchains',4,'out',struct('conv',struct('data',{tbl})));
                    checked=psyrat_checkconv(rel);
                    t.verifyEqual(checked.out.conv.converged,double(numel(names)==numel(expected)),f.file);
                    % A healthy stochastic diagnostic remains an accepted control.
                    for j=1:numel(names)
                        if ~ismember(names{j},expected); rows(j,2:6)={400,1,400,400,.01}; end
                    end
                    rel.out.conv.data=psyrat_diagnostic_table(rows,psyrat_diagnostic_metadata({}),src,f.data);
                    checked=psyrat_checkconv(rel);
                    t.verifyEqual(checked.out.conv.converged,1,f.file);
                end
            end
        end
        function sourceProofVerdictsSurviveSaveAndSidecar(t)
            folder=fullfile(t.projectRoot(),'tests','fixtures','constant_source','parent_order_cases');
            manifest=jsondecode(fileread(fullfile(folder,'manifest.json')));
            output=tempname; mkdir(output); t.addTeardown(@()rmdir(output,'s'));
            % Actual .psyrat and JSON files for the temporary-parent failure
            % and the ready-parent valid control, using frozen source algebra.
            for k=1:2
                f=manifest.fixtures(k); if iscell(f); f=f{1}; end
                names=cellstr(string(f.probed_names(:)));
                expected=cellstr(string(f.expected_bounded_rule_names(:)));
                src=fileread(fullfile(folder,f.file));
                rows=[names num2cell(nan(numel(names),5)); {'healthy',400,1,400,400,.01}];
                tbl=psyrat_diagnostic_table(rows,psyrat_diagnostic_metadata({}),src,f.data);
                rel=struct('engine','cmdstan','nchains',4,'out',struct('conv',struct('data',{tbl})));
                rel=psyrat_checkconv(rel);
                psyrat_data=struct('rel',rel,'proc',struct('savepath',output,'savename',sprintf('proof%d.psyrat',k)));
                saved=fullfile(output,psyrat_data.proc.savename); save(saved,'psyrat_data','-mat');
                loaded=load(saved,'-mat'); checked=psyrat_checkconv(loaded.psyrat_data.rel);
                expectedGate=double(numel(names)==numel(expected));
                t.verifyEqual(checked.out.conv.converged,expectedGate,f.file);
                t.verifyEqual(checked.out.conv.data,rel.out.conv.data);
                cfg=jsondecode(fileread(psyrat_write_runconfig(loaded.psyrat_data,psyrat_defaults())));
                t.verifyEqual(cfg.diagnostics.converged,expectedGate,f.file);
                t.verifyEqual(cfg.diagnostics.final_inference,'unassessed');
                t.verifyEqual(sort(string(cfg.diagnostics.excluded(:))),sort(string(expected(:))));
                savedInfo=cfg.diagnostics.diagnostics;
                if iscell(savedInfo); savedInfo=savedInfo{1}; end
                t.assertNumElements(savedInfo,1);
                t.verifyEqual(sort(string(savedInfo.permitted_constants(:))),sort(string(expected(:))));
                t.verifyEqual(savedInfo.model_source_hash,tbl{1,8}.model_source_hash);
                t.verifyEqual(savedInfo.model_data_hash,tbl{1,8}.model_data_hash);
                t.verifyEqual(savedInfo.contract_status,tbl{1,8}.contract_status);
                if ~expectedGate; t.verifyTrue(any(strcmp(cfg.diagnostics.invalid,'alias[1,1]'))); end
            end
        end
        function complexDiagnosticTokensAreMalformed(t)
            hdr='Mean MCSE StdDev 5% 50% 95% ESS_bulk ESS_tail ESS_bulk/s R_hat';
            file=[tempname '.csv']; t.addTeardown(@()delete(file));
            for token={'.01','1i'}
                row=['mu 0 ' token{1} ' 1 -1 0 1 400 400 10 1'];
                [rows,info]=psyrat_parse_stansummary_text({hdr,row},@(n)strcmp(n,'mu'));
                tbl=psyrat_diagnostic_table(rows,info,'parameters { real mu; }',struct());
                rel=struct('nchains',4,'out',struct('conv',struct('data',{tbl})));
                out=psyrat_checkconv(rel);
                t.verifyEqual(out.out.conv.converged,double(strcmp(token{1},'.01')));
                fid=fopen(file,'w');
                fprintf(fid,'name,MCSE,ESS_bulk,ESS_tail,R_hat\n"mu",%s,400,400,1\n',token{1});
                fclose(fid);
                if strcmp(token{1},'1i')
                    t.verifyError(@()psyrat_parse_stansummary_csv(file,@(~)true),'PsyRAT:stansummaryCsvRow');
                else
                    rows=psyrat_parse_stansummary_csv(file,@(~)true);
                    t.verifyEqual(rows{1,6},.01);
                end
            end
        end
        function nestedAndChainedIndexWritesFailClosed(t)
            base=['data { array[1] int<lower=1,upper=2> index; } ' ...
                'parameters { corr_matrix[2] R; real z; } ' ...
                'model { R ~ lkj_corr(2); z ~ normal(0,1); } ' ...
                'generated quantities { matrix[2,2] alias=block(R,1,1,2,2); '];
            control=psyrat_constant_contract([base '}'],struct('index',1),{'alias[1,1]'});
            t.verifyEqual(control.names,{'alias[1,1]'});
            % With index=1 these assignments make the formerly unit entry z.
            for lhs={'alias[index[1],1]','alias[index[index[1]],1]','alias[index[1]][1]'}
                src=[base lhs{1} '=z; }'];
                c=psyrat_constant_contract(src,struct('index',1),{'alias[1,1]'});
                t.verifyEmpty(c.names,lhs{1});
                rows={'alias[1,1]',NaN,NaN,NaN,NaN,NaN;'z',400,1,400,400,.01};
                tbl=psyrat_diagnostic_table(rows,psyrat_diagnostic_metadata({}),src,struct('index',1));
                rel=struct('nchains',4,'out',struct('conv',struct('data',{tbl})));
                out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0,lhs{1});
            end
        end
        function malformedTextCannotBecomeAnExemptConstant(t)
            hdr='Mean MCSE StdDev 5% 50% 95% ESS_bulk ESS_tail ESS_bulk/s R_hat';
            good='mu 0 .01 1 -1 0 1 400 400 10 1';
            for token={'nan','-nan','+nan','broken'}
                row=['R[1,1] 1 0 0 1 1 1 ' token{1} ' ' token{1} ' nan ' token{1}];
                [rows,info]=psyrat_parse_stansummary_text({hdr,good,row},@(n)any(strcmp(n,{'mu','R[1,1]'})));
                tbl=psyrat_diagnostic_table(rows,info,'generated quantities { corr_matrix[2] R; }',struct());
                rel=struct('nchains',4,'out',struct('conv',struct('data',{tbl})));
                out=psyrat_checkconv(rel);
                t.verifyEqual(out.out.conv.converged,double(~strcmp(token{1},'broken')),token{1});
            end
            [rows,info]=psyrat_parse_stansummary_text({hdr,good,good},@(n)strcmp(n,'mu'));
            tbl=psyrat_diagnostic_table(rows,info,'parameters { real mu; }',struct());
            rel=struct('nchains',4,'out',struct('conv',struct('data',{tbl})));
            out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0,'Duplicate monitored rows are malformed.');
        end
        function typedMatricesAndBounds(t)
            src=['parameters { cholesky_factor_corr[2] L; cholesky_factor_cov[2] V; } ' ...
                'generated quantities { corr_matrix[2] R; matrix[2,2] cor_fake; }'];
            names={'R[1,1]','R[2,2]','R[1,2]','R[3,3]','L[1,1]','L[1,2]', ...
                'L[2,2]','L[2,1]','V[1,1]','cor_fake[1,1]','R[0,0]'};
            c=psyrat_constant_contract(src,struct(),names);
            t.verifyEqual(c.names,{'R[1,1]','R[2,2]','L[1,1]','L[1,2]'});
            t.verifyEqual(c.values,[1 1 1 0]);
        end
        function actualDodAndGammaAliases(t)
            dirpath=fullfile(t.projectRoot(),'tests','golden','stan');
            c=psyrat_constant_contract(fileread(fullfile(dirpath,'dod_event.stan')), ...
                struct('C',4),{'L_id[1,1]','L_id[1,8]','L_id[2,2]', ...
                'cor_ID_Cell[4,4]','cor_ID_Cell[5,5]','id_varcov[1,1]'});
            t.verifyEqual(c.names,{'L_id[1,1]','L_id[1,8]','cor_ID_Cell[4,4]'});
            c=psyrat_constant_contract(fileread(fullfile(dirpath,'ic_gamma_dodiff.stan')), ...
                struct(),{'cor_id[4,4]','cor_id[5,5]','cor_id[1,2]'});
            t.verifyEqual(c.names,{'cor_id[4,4]'});
        end
        function actualDisabledResidualGraph(t)
            src=fileread(fullfile(t.projectRoot(),'tests','golden','stan','ic_diff_group_event_sserr.stan'));
            names={'rho[2]','wp_cov_ss[2]','rescor[1]','err_varcov[1,2]', ...
                'err_varcov[2,1]','err_varcov[1,1]','rescor_mu','sd_rescor','z_rescor[1]','rho[3]'};
            off=psyrat_constant_contract(src,struct('estimate_rescor',0,'NSUB',2),names);
            on=psyrat_constant_contract(src,struct('estimate_rescor',1,'NSUB',2),names);
            missing=psyrat_constant_contract(src,struct(),names);
            t.verifyEqual(off.names,names(1:5)); t.verifyEqual(off.values,zeros(1,5));
            t.verifyEmpty(on.names); t.verifyEmpty(missing.names);
            % The stricter source/data proof retains the valid one-subject
            % declaration boundary and refuses an impossible zero dimension.
            oneNames={'rho[1]','wp_cov_ss[1]','rescor[1]'};
            one=psyrat_constant_contract(src,struct('estimate_rescor',0,'NSUB',1),oneNames);
            zero=psyrat_constant_contract(src,struct('estimate_rescor',0,'NSUB',0),oneNames);
            t.verifyEqual(one.names,oneNames); t.verifyEqual(one.values,[0 0 0]);
            t.verifyEmpty(zero.names);
        end
        function reassignedAliasesAndZeroGraphsFailClosed(t)
            src=['generated quantities { corr_matrix[2] parent; matrix[2,2] alias; ' ...
                'alias=block(parent,1,1,2,2); alias[1,1]=2; }'];
            c=psyrat_constant_contract(src,struct(),{'alias[1,1]'}); t.verifyEmpty(c.names);
            src=fileread(fullfile(t.projectRoot(),'tests','golden','stan','ic_diff_group_event_sserr.stan'));
            changed=strrep(src,'err_varcov = rep_matrix(0, 2, 2);','err_varcov = rep_matrix(rescor_mu, 2, 2);');
            c=psyrat_constant_contract(changed,struct('estimate_rescor',0,'NSUB',2),{'err_varcov[1,2]'});
            t.verifyEmpty(c.names);
            src=[src newline 'generated quantities { rescor[1]=normal_rng(0,1); }'];
            c=psyrat_constant_contract(src,struct('estimate_rescor',0,'NSUB',2),{'rescor[1]','rho[1]'});
            t.verifyEmpty(c.names);
        end
        function independentAliasFixtures(t)
            folder=fullfile(t.projectRoot(),'tests','fixtures','constant_aliases');
            manifest=jsondecode(fileread(fullfile(folder,'manifest.json')));
            for k=1:numel(manifest.fixtures)
                f=manifest.fixtures(k); if iscell(f); f=f{1}; end; data=struct();
                if isfield(f,'data')&&isstruct(f.data); data=f.data; end
                c=psyrat_constant_contract(fileread(fullfile(folder,f.file)),data,f.probed_names);
                expected=f.expected_bounded_rule_names;
                if isempty(expected); expected=cell(0,1); else; expected=cellstr(expected(:)); end
                t.verifyEqual(c.names(:),expected,f.file);
            end
        end
        function validConstantAndInvalidRows(t)
            c={'name','n_eff','r_hat';'R[1,1]',NaN,NaN;'mu',400,1};
            m=psyrat_diagnostic_metadata({});
            contract=psyrat_constant_contract('generated quantities { corr_matrix[2] R; }',struct(),c(2:end,1));
            m.permitted_constants=contract.names; m.contract_status=contract.status;
            c{1,5}=m;
            rel=struct('nchains',4,'out',struct('conv',struct('data',{c})));
            out=psyrat_checkconv(rel);
            t.verifyEqual(out.out.conv.converged,1);
            t.verifyEqual(out.out.conv.excluded,{'R[1,1]'});
            t.verifyEqual(out.out.conv.final_inference,'unassessed');
            for pair={[Inf,1],[400,0],[400,-1],[400,NaN],[400,1i]}
                rel.out.conv.data(3,2:3)=num2cell(pair{1});
                out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0);
            end
            m.expected_names={'mu','missing_stochastic'}; c{1,5}=m;
            rel.out.conv.data=c; out=psyrat_checkconv(rel);
            t.verifyEqual(out.out.conv.converged,0);
        end
        function malformedNamesAndNestedMissingFailClosed(t)
            a={'name','n_eff','r_hat';42,400,1;'mu',400,1};
            m=psyrat_diagnostic_metadata({}); m.expected_names={'mu','x','y'}; a{1,5}=m;
            b={'name','n_eff','r_hat';'mu',400,1}; b{1,5}=m;
            rel=struct('nchains',4,'out',struct('conv',struct('data',{{a,b}})));
            out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0);
            t.verifyTrue(all(ismember({'x','y','malformed_diagnostic_row'},out.out.conv.invalid)));
            rel.out.conv.data=42; out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0);
            rel.out.conv.data={b;b}; out=psyrat_checkconv(rel); t.verifyEqual(out.out.conv.converged,0);

        end
        function parserQualityAndMalformedRows(t)
            modern={'Mean MCSE StdDev 5% 50% 95% ESS_bulk ESS_tail ESS_bulk/s R_hat', ...
                'mu 0 .01 1 -1 0 1 140 276 10 1.0'};
            [r,m]=psyrat_parse_stansummary_text(modern,@(n)strcmp(n,'mu'));
            t.verifyEqual(r{1,2},140); t.verifyEqual(r{1,4},140); t.verifyEqual(r{1,5},276);
            t.verifyEqual(m.format,'modern_text'); t.verifyTrue(m.degraded);
            t.verifyEqual(m.precision_significant_digits,2);
            legacy={'Mean MCSE StdDev 5% 50% 95% N_Eff N_Eff/s R_hat', ...
                'mu 0 .01 1 -1 0 1 140 10 1.0'};
            [r,m]=psyrat_parse_stansummary_text(legacy,@(n)strcmp(n,'mu'));
            t.verifyEqual(r{1,2},140); t.verifyTrue(isnan(r{1,4})&&isnan(r{1,5}));
            t.verifyEqual(m.format,'legacy_text');
            f=[tempname '.csv']; t.addTeardown(@()delete(f));
            fid=fopen(f,'w'); fprintf(fid,'name,ESS_bulk,ESS_tail,R_hat\n"mu",400,400,1\n"sig_id",400\n'); fclose(fid);
            t.verifyError(@()psyrat_parse_stansummary_csv(f,@(~)true),'PsyRAT:stansummaryCsvRow');
        end
    end
end
