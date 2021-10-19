function [t,t_matrix,l_matrix]=depolymerisation(fr,k_depoly,trec,l)
t=0;
t_matrix(1)=0;l_matrix(1)=0;
for i=1:fr 
r1=rand;
tau=(1/k_depoly)*log(1/r1);
t=t+tau;
trec=trec+tau;
l=l-0.1;
t_matrix(i)=trec;
l=round(l,1);
l_matrix(i)=l;
end
end