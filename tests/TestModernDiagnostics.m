classdef TestModernDiagnostics < PsyRATTestBase
    %PBA-03: exact arrays and primary expectations were generated independently
    %by R posterior 1.6.0 BEFORE production comparison. Never fit the oracle.
    methods (Test)
        function trustedReferenceAgreement(t)
            dirpath=fullfile(t.projectRoot(),'tests','fixtures','rank_diagnostics');
            expected=readtable(fullfile(dirpath,'posterior_binary_reference.csv'),'TextType','string');
            guards={'same_constant','different_constants','mixed_one_constant','nonfinite'};
            for k=1:height(expected)
                name=char(expected.name(k));
                fid=fopen(fullfile(dirpath,[name '.bin']),'rb','ieee-le');
                x=fread(fid,[expected.iterations(k),expected.chains(k)],'double'); fclose(fid);
                [~,c]=psyrat_hmc_sample('diagnostics',x,{'quantity'},.7*ones(1,size(x,2)));
                if ismember(name,guards)
                    t.verifyEqual(c{2,3},Inf,name);
                    t.verifyEqual(c{2,2},0,name);
                    continue;
                end
                t.verifyEqual(c{2,3},expected.rhat(k),'AbsTol',1e-10,name);
                t.verifyEqual(c{2,4},expected.ess_bulk(k),'RelTol',1e-10,name);
                t.verifyEqual(c{2,5},expected.ess_tail(k),'RelTol',1e-10,name);
                t.verifyEqual(c{2,6},expected.mcse_mean(k),'RelTol',1e-10,name);
            end
        end

        function invalidAndTooShortDiagnosticsFail(t)
            for n=1:5
                [~,c]=psyrat_hmc_sample('diagnostics',randn(n,4),{'x'},.7*ones(1,4));
                t.verifyEqual(c{2,3},Inf); t.verifyEqual(c{2,2},0);
            end
            for a={NaN,Inf,-.1,1.1,.7+.1i}
                [~,c]=psyrat_hmc_sample('diagnostics',randn(100,4),{'x'},[.7 .7 .7 a{1}]);
                t.verifyEqual(c{2,3},Inf);
            end
        end

        function stochasticConstantFailsBesideHealthyRow(t)
            c={'name','n_eff','r_hat';'sig_id',NaN,NaN;'mu',400,1};
            rel=struct('nchains',4,'out',struct('conv',struct('data',{c})));
            out=psyrat_checkconv(rel);
            t.verifyEqual(out.out.conv.converged,0);
            t.verifyEmpty(out.out.conv.excluded);
        end
    end
end
