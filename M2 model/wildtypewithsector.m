tic
Awmt=load('/home/samlesh/Documents/matlab_R2020b_glnxa64/samlesh/data folder/lengthEvo.mat');
[lmin]=0.6;%findmin will find the minimum mitochondria in wild type
lint=lmin;%quatising step
stormatrix(20,20,21)=0;%initial bining and generation of initial conditon
for k=1:21   
    temp=Awmt.lengthEvo(k,:);
    [new2]=binning(temp,lmin,lint);
    stormatrix(:,:,k)=new2';   
    Mwilswiths(:,:,k)=wmt_with_sector(stormatrix(:,:,k));
toc
end
save('wildtypemtwithsector.mat','Mwilswiths');
