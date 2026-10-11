classdef TestDodPluginCell < PsyRATTestBase
    % Bounded follow-up to independent 2B review: freeze the retained floating
    % operation order and exercise nonfinite component reasons directly.
    methods (Test)
        function testLegacyBoundaryIsNotRoundedAway(testCase)
            bp=diag([1 0 0 0]); bt=zeros(4); bt(1:2,1:2)=[1 2;2 4];
            % Independent scalar evaluation of the legacy diagonal-then-pair
            % order. The exact algebraic value is zero; no tolerance policy
            % is inferred from that identity. This frozen IEEE double result
            % pins the actual value, rather than reading it back from output.
            n1=29; n2=58;
            legacy=1/n1+4/n2-4/(2/(1/n1+1/n2));
            testCase.assertEqual(legacy,-1.3877787807814457e-17);
            s=psyrat_dod_plugin_cell(bp,bt,ones(1,4),[n1 n2 1 1],1);
            testCase.verifyEqual(s.bt_scaled,legacy);
            testCase.verifyFalse(s.admissible);
            testCase.verifyTrue(contains(s.reason,'between-trial contrast variance negative'));
            testCase.verifyEqual(s.value,n2/(3*n2+3),'AbsTol',1e-15);
            % An exactly representable zero remains a valid boundary.
            control=psyrat_dod_plugin_cell(bp,bt,ones(1,4),[1 2 1 1],1);
            testCase.verifyEqual(control.bt_scaled,0);
            testCase.verifyTrue(control.admissible);
        end

        function testNonfiniteComponentsAreRetainedAndNamed(testCase)
            bp=diag([1 0 0 0]); bt=eye(4); er=ones(1,4);
            control=psyrat_dod_plugin_cell(bp,bt,er,ones(1,4),1);
            testCase.verifyTrue(control.admissible);
            bpn=bp; bpn(1,1)=NaN;
            s=psyrat_dod_plugin_cell(bpn,bt,er,ones(1,4),1);
            testCase.verifyTrue(isnan(s.u)); testCase.verifyTrue(isnan(s.value));
            testCase.verifyFalse(s.admissible);
            testCase.verifyTrue(contains(s.reason,'universe-score variance not finite'));
            btn=bt; btn(1,1)=Inf;
            s=psyrat_dod_plugin_cell(bp,btn,er,ones(1,4),2);
            % The G ratio is still finite because trial variance is outside
            % its denominator; the established evaluation record still flags it.
            testCase.verifyEqual(s.value,1/5);
            testCase.verifyEqual(s.bt_scaled,Inf);
            testCase.verifyFalse(s.admissible);
            testCase.verifyTrue(contains(s.reason,'between-trial contrast variance not finite'));
            testCase.verifyTrue(contains(s.reason,'not in the generalizability denominator'));
            ern=er; ern(1)=NaN;
            s=psyrat_dod_plugin_cell(bp,bt,ern,ones(1,4),1);
            testCase.verifyTrue(isnan(s.er_scaled));
            testCase.verifyTrue(isnan(s.denominator));
            testCase.verifyFalse(s.admissible);
            testCase.verifyTrue(contains(s.reason,'residual contrast variance not finite'));
        end
    end
end
