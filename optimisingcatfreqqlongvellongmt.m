%optimising the output of multivariopticatfreqelongvellongmt using genetic
%algorithm
options=optimoptions('ga','PlotFcn','gaplotbestf','Display','iter','UseParallel','always','MaxGenerations',200);
[x,y]=ga(@multivariopticatfreqelongvellongmt,4,[],[],[],[],[0.01 0.01 0.01 0.01],[10 10 10 10],[],options);
save('longmtGDPincluded1.mat','x','y');