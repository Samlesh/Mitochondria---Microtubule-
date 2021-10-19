%optimising the output of multivariopticatfreqelongvelwildtypemt using genetic
%algorithm
options=optimoptions('ga','PlotFcn','gaplotbestf','Display','iter','UseParallel','always','MaxGenerations',200);
[x,y]=ga(@multivariopticatfreqelongvelwildtypemt,4,[],[],[],[],[0.01 0.01 0.01 0.01],[10 10 10 10],[],options);
save('wildtypemtGDPincludedmodified.mat','x','y');