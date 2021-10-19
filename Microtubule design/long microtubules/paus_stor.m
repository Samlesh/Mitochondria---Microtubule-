function [pa1]=paus_stor(pai,pao)
pa1=pao;
c=0;
if (numel(pao)==1)
y=1;
z=numel(pai);
else
    y=numel(pao)+1;
    z=y-1+numel(pai);
end
for i=y:z
    c=c+1;
    pa1(i)=pai(c);
end
end