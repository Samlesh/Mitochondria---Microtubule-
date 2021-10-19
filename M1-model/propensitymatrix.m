function[propensity]=propensitymatrix(n0,parameter,random_x)
%[stochiomatrixchk,random_x]=stochiometricmatrix();
nzero(1)=0;
kb1=parameter(1);%dissimilar fission
kb2=parameter(2);%similar fission
ka1=parameter(3);%fusion of different lengths
ka2=parameter(4);%fusion of similar lengths
n=n0;
propensity(numel(random_x(:,1)))=0;
for i=1:numel(random_x(:,1))
temp(20)=0; 
temp(:)=random_x(i,:);
nzero=find(temp);
[~,lzero]=find(temp<0);
if (numel(lzero)==1 && temp(lzero)==-1) %fission
h=find(temp==2);
l=isempty(h);
if l==1
    propensity(i)=kb1*n(lzero)*((lzero-1));%2is not there, normal fission
else
    propensity(i)=kb2*n(lzero)*((lzero-1));%2 is there similar fisiion
end      
elseif(numel(lzero)==1 && temp(lzero)==-2)%fusion of similar lengyth
propensity(i)=ka2*n(lzero(1))*(n(lzero(1))-1);
elseif(numel(lzero)==2)%fusion of different length
propensity(i)=ka1*n(lzero(1))*(n(lzero(2)));   
end
clear('lzero');
end
end

