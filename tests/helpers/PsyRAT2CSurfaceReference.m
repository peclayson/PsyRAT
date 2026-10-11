function t=PsyRAT2CSurfaceReference(alpha,sd,scale,scale0,scaleSD,covariance,facets,ni,no,kind,target)
% Independent planted surface reference. No production converters or ratios.
% sd order P,O,I,PO,PI,OI. scale=1 lognu; scale=2 LS.
% Frozen algebra in the 2C SCIENTIFIC_CONTRACT and independent Decimal oracle.
assert(numel(sd)==6 && ismember(scale,[1 2]) && ismember(facets,[1 2]));
if nargin<11; target='population_marginal'; end
k=PsyRAT2CMeanReference(alpha,sd(1),sd(2),sd(3),sd(4),sd(5),sd(6));
V=sum(sd.^2); W=V-sd(1)^2;
if scale==1
    if strcmp(target,'population_marginal')
        C=2*exp(2*alpha+2*V-scale0+scaleSD^2/2-2*covariance);
    else
        assert(strcmp(target,'typical_person')); C=2*exp(2*alpha+2*W-scale0);
    end
else
    if strcmp(target,'population_marginal')
        C=exp(2*scale0+2*scaleSD^2);
    else
        assert(strcmp(target,'typical_person')); C=exp(2*scale0);
    end
end
P=k(1); O=k(2); I=k(3); PO=k(4); PI=k(5); OI=k(6);
if facets==1
    E=PI+C; signal=P; relative=E/ni; absolute=relative+I/ni;
else
    E=k(7)+C;
    switch kind
        case 1; signal=P+PO/no; relative=PI/ni+E/(ni*no); extra=I/ni+OI/(ni*no);
        case 2; signal=P+PI/ni; relative=PO/no+E/(ni*no); extra=O/no+OI/(ni*no);
        case 3; signal=P; relative=PI/ni+PO/no+E/(ni*no); extra=I/ni+O/no+OI/(ni*no);
        otherwise; error('Invalid allocation');
    end
    absolute=relative+extra;
end
t=struct('G',signal/(signal+relative),'D',signal/(signal+absolute),...
    'residual',E,'observation',C,'components',k,'sem_rel',sqrt(relative),...
    'sem_abs',sqrt(absolute),'surface_target',target);
end
