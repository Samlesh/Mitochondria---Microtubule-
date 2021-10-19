function [lm1]=len_stor(lmi1,lo1)
lm1=lo1;
c=0;
if (numel(lo1)==1)
y=1;
z=numel(lmi1);
else
    y=numel(lo1)+1;
    z=y-1+numel(lmi1);
end
for i=y:z
    c=c+1;
    lm1(i)=lmi1(c);
end
end