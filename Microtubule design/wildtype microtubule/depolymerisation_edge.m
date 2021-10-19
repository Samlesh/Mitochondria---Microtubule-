function [t,t_keep,l_keep]=depolymerisation_edge(k_depolyedge,l,trec)
t=0;
t_keep(2)=0;l_keep(2)=0;
for i=1:58
r1=rand;
tau=(1/k_depolyedge)*log(1/r1);
t=t+tau;
l=l-0.1;
trec=trec+tau;
t_keep(i)=trec;
l=round(l,1);
l_keep(i)=l;
end
end