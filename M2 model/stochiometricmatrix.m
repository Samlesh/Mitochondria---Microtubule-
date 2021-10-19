function [stochiomatrixchk,random_x]=stochiometricmatrix()
stochiomatrixfission(100,20)=0;
for i=1:20
    stochiomatrixfission(1,i)=i;
end
k=numel(stochiomatrixfission(1,:))+1;%k=61
n=1;e=1;prev=2;
%for j=2:numel(stochiomatrix(:,1))
for i=1:(numel(stochiomatrixfission(1,:))-1)
    m=1;
d=stochiomatrixfission(1,k-i);%stores 20
tempmatrix(1,20)=0;
if (d>3)
if (mod(d,2)==0)
c=d/2;
else
c=d/2-1;
c=c+1;
end
for j=e:c
    n=n+1;
    if(j~=d/2)
    tempmatrix(j,d)=-1;
    tempmatrix(j,m)=1;
    tempmatrix(j,d-m)=1;   
    else
      tempmatrix(j,d)=-1;   
      tempmatrix(j,m)=2;
    end
    m=m+1;
end
stochiomatrixfission(prev:n,:)=tempmatrix;
prev=n+1;
else
stochiomatrixfission(prev,1:20)=0;
if(d==3)
stochiomatrixfission(prev,d)=-1;
stochiomatrixfission(prev,d-1)=1;
stochiomatrixfission(prev,d-2)=1;
else
stochiomatrixfission(prev,d)=-1;
stochiomatrixfission(prev,d-1)=2;    
end
prev=prev+1;
end
clear('tempmatrix');
end
stochiomatrixfusion(numel(stochiomatrixfission(:,1)),20)=0;c=0;
for j=1:numel(stochiomatrixfission(:,1))
stochiomatrixfusion(j,:)=-1.*stochiomatrixfission(end-c,:);
c=c+1;
end
% stochiomatrix(numel(stochiomatrixfission(:,1)),20)=0;
stochiomatrix=[stochiomatrixfission;stochiomatrixfusion];

stochiomatrixchk=stochiomatrix(2:end,:);
random_x = stochiomatrix(randperm(size(stochiomatrixchk, 1)), :);
end
%end