classdef PsyRATTestLmeProxy
    %Test-only proxy: replace one extracted covariance, keep the real fit.
    properties
        Fit
        Value
    end
    methods
        function o=PsyRATTestLmeProxy(f,v); o.Fit=f; o.Value=v; end
        function [psi,mse]=covarianceParameters(o)
            [psi,mse]=covarianceParameters(o.Fit);
            psi{1}(1,1)=o.Value;
        end
        function fe=fixedEffects(o); fe=fixedEffects(o.Fit); end
    end
end
