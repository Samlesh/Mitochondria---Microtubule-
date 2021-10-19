function [tt1]=time_stor(tti,tto)
tt1=tto;
c=0;
if (numel(tto)==1)
y=1;
z=numel(tti);
else
    y=numel(tto)+1;
    z=y-1+numel(tti);
end
for i=y:z
    c=c+1;
    tt1(i)=tti(c);
end
end