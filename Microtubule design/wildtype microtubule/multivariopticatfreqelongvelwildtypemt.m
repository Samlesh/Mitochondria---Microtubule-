function out=multivariopticatfreqelongvelwildtypemt(x)
a=x(1);b=x(2);k_depoly=x(3);k_depolyedge=x(4);
%a-growth rate constant %b-hydrolysis rate constant
%%k_depolyedge-depolymerisation rate constant
[sl1,sll1,sl2,sll2,sl3,sll3,sl4,sll4,tts1,tts2,tts3,tts4,ttts1,ttts2,ttts3,ttts4,cc1,nume1,nume2,nume3,nume4,cat_l,catfreq1,catfreq2,catfreq3,catfreq4,ttv1,ttv2,ttv3,ttv4,pa1,pa2,pa3,pa4,lengths_m1,lengths_m2,lengths_m3,lengths_m4,vavgs_m1,vavgs_m2,vavgs_m3,vavgs_m4]=pri_mic_code_wmt(a,b,k_depolyedge);
%take input function takes length vs time as input and generates elongation velocity , 
% catastrophe frequency and shortening velocity   
[Vavg_length_m1,~,mcat_m1,~,~,vavge_m1,lengthe_m1]=take_input(sl1,tts1,nume1,pa1,catfreq1);
[Vavg_length_m2,~,mcat_m2,~,~,vavge_m2,lengthe_m2]=take_input(sl2,tts2,nume2,pa2,catfreq2);
[Vavg_length_m3,~,mcat_m3,~,~,vavge_m3,lengthe_m3]=take_input(sl3,tts3,nume3,pa3,catfreq3);
[Vavg_length_m4,~,mcat_m4,~,~,vavge_m4,lengthe_m4]=take_input(sl4,tts4,nume4,pa4,catfreq4);
[vavgsm1]=average_velocity(lengths_m1,vavgs_m1);[vavgsm2]=average_velocity(lengths_m2,vavgs_m2);
[vavgsm3]=average_velocity(lengths_m3,vavgs_m3);[vavgsm4]=average_velocity(lengths_m4,vavgs_m4);
ve=[vavge_m1;vavge_m2;vavge_m3;vavge_m4];vef=mean(ve,1);
vs=[vavgsm1;vavgsm2;vavgsm3;vavgsm4];vsf=mean(vs,1);
refl=0.1:0.1:5.8;
[le1,~]=histc(lengthe_m1,refl);[ls1,~]=histc(lengths_m1,refl);le1=le1./sum(le1);ls1=ls1./sum(ls1);
[le2,~]=histc(lengthe_m2,refl);[ls2,~]=histc(lengths_m2,refl);le2=le2./sum(le2);ls2=ls2./sum(ls2);
[le3,~]=histc(lengthe_m3,refl);[ls3,~]=histc(lengths_m3,refl);le3=le3./sum(le3);ls3=ls3./sum(ls3);
[le4,~]=histc(lengthe_m4,refl);[ls4,~]=histc(lengths_m4,refl);le4=le4./sum(le4);ls4=ls4./sum(ls4);
le=[le1;le2;le3;le4];lef=mean(le,1);
ls=[ls1;ls2;ls3;ls4]; lsf=mean(ls,1);
m1ve=vavge_m1.*le1;m1vef=sum(m1ve);m1vs=vavgsm1.*ls1;m1vsf=sum(m1vs);
m2ve=vavge_m2.*le2;m2vef=sum(m2ve);m2vs=vavgsm2.*ls2;m2vsf=sum(m2vs);
m3ve=vavge_m3.*le3;m3vef=sum(m3ve);m3vs=vavgsm3.*ls3;m3vsf=sum(m3vs);
m4ve=vavge_m4.*le4;m4vef=sum(m4ve);m4vs=vavgsm4.*ls4;m4vsf=sum(m4vs);
venew=vef.*lef;vsnew=vsf.*lsf;
vechk=sum(venew);vschk=sum(vsnew);
%  [pavg1,finalp1,finalprob1]=pause_time(sll1,ttts1); [pavg2,finalp2,finalprob2]=pause_time(sll2,ttts2);
%  [pavg3,finalp3,finalprob3]=pause_time(sll3,ttts3); [pavg4,finalp4,finalprob4]=pause_time(sll4,ttts4);
Ve1=Vavg_length_m1;Ve2=Vavg_length_m2;Ve3=Vavg_length_m3;Ve4=Vavg_length_m4;Veavg=(Ve1+Ve2+Ve3+Ve4)/4;
%Vs1=Vavg_length1_m1;Vs2=Vavg_length1_m2;Vs3=Vavg_length1_m3;Vs4=Vavg_length1_m4;Vsavg=(Vs1+Vs2+Vs3+Vs4)/4;
%P1=pavg1;P2=pavg2;P3=pavg3;P4=pavg4;Pavg=(P1+P2+P3+P4)/4;
C1=mcat_m1;C2=mcat_m2;C3=mcat_m3;C4=mcat_m4;Cavg=(C1+C2+C3+C4)/4;
[length1,time1,index1]=nucleation(sll1,ttts1);[length2,time2,index2]=nucleation(sll2,ttts2);
[length3,time3,index3]=nucleation(sll3,ttts3);[length4,time4,index4]=nucleation(sll4,ttts4);
timestop1=find(~time1);timestop2=find(~time2);timestop3=find(~time3);timestop4=find(~time4);
avml1=averagelength(length1,time1,index1);avml2=averagelength(length2,time2,index2);avml3=averagelength(length3,time3,index3);
avml4=averagelength(length4,time4,index4);
avml=[avml1,avml2,avml3,avml4];
lchk=mean(avml);
m1=0.05; %m1 is the experimental elongation velocity
m2=2.62; %m2 is the experimental observed mean length
%m3=28;
m4=0.14;  %m4 is the experimental shrinkage velocity
%m3=8.5/60;
out=sqrt(((vechk-m1)/m1)^2+((lchk-m2)/m2)^2+((vschk-m4)/m4)^2);
%out represents the difference between experiments and simulated values and
%is the output of this function
end