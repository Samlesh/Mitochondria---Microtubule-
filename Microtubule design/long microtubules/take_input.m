function [Vavg_length,mpe,mcat,vsto,Vavg_sd,vavg,lls]=take_input(sl1,tts1,nume,pal,catfreq)
%elongation velocity
[vstor,lls]=velelong(sl1,tts1,nume);
[vavg]=average_velocity(lls,vstor);
[vsto]=dist_realisation(sl1,tts1,nume);
Vavg_length=sum(vavg*(0.1/5.8));Vavg_sd=std(vavg);
%shortening velocity
% [vshort,lls1]=velshort(sl1,tts1,ttv1,nume);
% [vavg1]=average_velocity(lls1,vshort);
% [vsto1]=dist_realisation1(sl1,tts1,ttv1,nume);
% Vavg_length1=sum(vavg1*(0.1/5.8));Vavg_sd1=std(vavg1);
mcat=mean(catfreq);
[pause]=p_chk(sl1,pal,nume);
mpe=pause;
end