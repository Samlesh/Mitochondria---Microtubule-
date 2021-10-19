function [propensity]= propensitymatrix_mt(ct,con,nGTP,kgr,khy,kde)
tubestate=ct;
propensity(3)=0;
if tubestate==1  %%the tip is in GTP state 
propensity(1)=kgr*con;
else
    propensity(1)=0;
end
if nGTP~=0       %%GTP molecules are present 
propensity(2)=khy;
else
    propensity(2)=0;
end
if tubestate==2 
propensity(3)=kde;
else
    propensity(3)=0;
end
end