function[pt]=pt_cal(nGTP,khy,chkcat)
if(nGTP>1)
if (chkcat==1)
count=10*nGTP;
count=round(count);
elseif (chkcat==0) 
count=58;
end
pt=0;
for i=1:count
r1=rand;
tau=(1/khy)*log(1/r1);
pt=pt+tau;
end
else
    pt=0;
end
end
