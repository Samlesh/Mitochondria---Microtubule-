function vsto=dist_realisation(sl1,tts1,nume)
k1=1;k2=0;vsto(1)=0;lls(1)=0;
count=0;f=0;ll=1;
for u=1:numel(nume)
    k2=k2+nume(u);
    XX=tts1(k1:k2);%stores time for a particular realisation
    YY=sl1(k1:k2);%stores length for a particular realisation
    c=0;cz(1)=0;
    vsim(1)=0;
    lenstor(1)=0;
    f=sum(YY==0);
    if (f==1)
        vsim(1)=YY(1,end)/XX(1,end);
    end
    count=1;h=0;
    for j=2:numel(YY)
        if (round(YY(j),1)==0)
            h=h+1;
            l=YY(j-1)-YY(count);
            t=XX(j-1)-XX(count);
            vsim(h)=l/t;
            lenstor(h)=l;
            count=j;
        elseif (j==numel(YY)&&(round(YY(end),1)~=0))
            h=h+1;
            l=YY(j)-YY(count);
            t=XX(j)-XX(count);
            vsim(h)=l/t;
            lenstor(h)=l;
            count=j;
        end
    end
 vsto(u)=mean(vsim);
k1=k2+1;
end
