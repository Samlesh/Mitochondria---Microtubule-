tic
Almt=load('lengthEvo_G3B.mat');%load the raw long microtubule data
%Almt stores the mitochondrial length values bserved in experiments for 
%long microtubule cells 
[lmin]=0.6;%findmin will find the minimum mitochondria in long 
           %microtubule case
lint=lmin;%quatising step
stormatrix(20,20,15)=0;%initial bining and generation of initial conditon
% bining is the function that takes input mitochondrial lengths and
% generates the number of mitochondria present in each of the 20 length
% fragments
for k=1:15
    %there are 15 cells in the experiments having long microtubules 
    temp=Almt.lengthEvo_G3B(k,:);
    [new2]=binning(temp,lmin,lint);
    %new2 stores the number of mitochondria present in each length fragment
    %for all time points
    stormatrix(:,:,k)=new2';   
    %stormatrix(20X20x15) is a 3D matrix where rows are are the time points
    %(20 time points separated by 12 seconds in experiments), the columns
    % are number of mitochondria present in each of the 20 length
    % fragments and the 3rd dimension k represent the cell whose length is
    % being stored
    Mlilsws(:,:,k)=lmt_without_sector(stormatrix(:,:,k));
    %lmt_without_sector returns the number of mitochondria at every time
    %obtained from simulation using the time=0 data from the
    %stormatrix(:,:,k) for the kth cell
toc
end
save('longmtwithoutsector.mat','Mlilsws');
%Mlilsws saves the resulatant mitochondrial number with time in the matrix
%longmtwithoutsector.mat