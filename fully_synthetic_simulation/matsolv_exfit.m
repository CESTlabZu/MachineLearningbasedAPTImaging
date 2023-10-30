function vals = matsolv_exfit(pars, x, Delta)


vals=zeros(length(x),1);

fs=pars(1);
ksw=pars(2);
R2b=pars(3);
res_off=Delta(4);

for ii=1:1:length(x)
    w1=Delta(1);
    offset=x(ii);
p=fs.*ksw.*w1.^2./(w1.^2+ksw.*(R2b+ksw)+(offset-res_off).^2.*ksw./(R2b+ksw));
 vals(ii)=p;
end


nnnn=length(x);

 
 
for ii=1:1:length(x)
    w1=Delta(2);
    offset=x(ii);
p=fs.*ksw.*w1.^2./(w1.^2+ksw.*(R2b+ksw)+(offset-res_off).^2.*ksw./(R2b+ksw));
vals(nnnn+ii)=p;
end


 
for ii=1:1:length(x)
    w1=Delta(3);
    offset=x(ii);
p=fs.*ksw.*w1.^2./(w1.^2+ksw.*(R2b+ksw)+(offset-res_off).^2.*ksw./(R2b+ksw));
vals(2*nnnn+ii)=p;
end




vals=vals';
    end
