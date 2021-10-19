%optimising the output of multivariopticatfreqelongvelshortmt using genetic
%algorithm
options=optimoptions('ga','Display','iter','UseParallel','always');
[x,y]=ga(@multivarioptielongationvelocityshortmt,4,[],[],[],[],[0.01 0.01 1 5],[1 1 5 10],[],options);
save('resultshortmtelongationvel.mat','x','y');