function [cc1]=conc(cci,cco)
cc1=cco;
c=0;
if (numel(cco)==1)
y=1;
z=numel(cci);
else
    y=numel(cco)+1;
    z=y-1+numel(cci);
end
for i=y:z
    c=c+1;
    cc1(i)=cci(c);
end
end