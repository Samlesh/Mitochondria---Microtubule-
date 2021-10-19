function [t,t_matrix,l_matrix]=depolymerisation_GTPseed(k_depoly,trec,l)
t_matrix(1)=0;l_matrix(1)=0;
count=10*(l);
count=round(count);
for i=1:count
r1=rand;

tau=(1/k_depoly)*log(1/r1);

% t=t+tau;
trec=trec+tau;
l=l-0.1;
t_matrix(i)=trec;
l=round(l,1);
l_matrix(i)=l;
end
t=trec;
end