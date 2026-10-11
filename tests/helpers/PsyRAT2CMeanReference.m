function k=PsyRAT2CMeanReference(alpha,sp,so,si,spo,spi,soi)
% Independent support-union reference frozen before 2C production comparison.
% Column order is P,O,I,PO,PI,OI,PIO. Test-only: no production helper calls.
% Derivation: the 2C science record (local development evidence, not shipped).
v=[sp(:).^2 so(:).^2 si(:).^2 spo(:).^2 spi(:).^2 soi(:).^2];
B=exp(2*alpha(:)+sum(v,2)); q=expm1(v);
masks=[1 2 4 3 5 6]; bySupport=zeros(size(v,1),7);
for subset=1:63
    term=B; support=0;
    for factor=1:6
        if bitget(subset,factor)
            support=bitor(support,masks(factor));
            term=term.*q(:,factor);
        end
    end
    bySupport(:,support)=bySupport(:,support)+term;
end
k=bySupport(:,[1 2 4 3 5 6 7]);
end
