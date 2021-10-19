tic
Awmt=load('lengthEvo.mat');
%Awmt stores the mitochondrial length values observed in experiments for 
%short microtubule cells
%microtubule lengths for sector 1 and 2
load('wmtsector1.mat');
load('wmtsector2.mat');
[lmin]=0.6;%findmin will find the minimum mitochondria in short
           %microtubule case
lint=lmin;%quatising step
stormatrix(20,20,21)=0;%initial bining and generation of initial conditon
% bining is the function that takes input mitochondrial lengths and
% generates the number of mitochondria present in each of the 20 length
% fragments
for k=1:21
    %there are 21 cells in the experiments having wildtype microtubules 
    temp=Awmt.lengthEvo(k,:);
    [new2]=binning(temp,lmin,lint);
    %new2 stores the number of mitochondria present in each length fragment
    %for all time points
    stormatrix(:,:,k)=new2';
    %stormatrix(20X20x21) is a 3D matrix where rows are are the time points
    %(20 time points separated by 12 seconds in experiments), the columns
    % are number of mitochondria present in each of the 20 length
    % fragments and the 3rd dimension k represent the cell whose length is
    % being stored
end
%76 is the mean number of mitochondria acros all cells at time t=0
initialwmt=stormatrix(:,:,4);%stores the time evolution of cell 4
n0=initialwmt(:,1);%stores the initial condition
picks=find(n0);stores(1:numel(picks))=0;
%picks is the matrix to store non-zero mitochondria present corresponding
%to respective length fragments in the initial state
for i=1:numel(picks)
    r=randi(2,1);
    stores(i)=r;
end
%stores is the matrix which consists the sector in which the mitochondria
%are present 
C=unique(stores);
%C represents the number of mitochondria present in sector 1 and sector 2
initialc(1:numel(n0),1:numel(C))=0;
%%segregation of initial condition if the entire cell into intial condition
%%of sector 1 and sector 2
for j=1:numel(stores)
    f111=stores(j);
    f222=picks(j);
    for i=1:numel(C)
        if (f111==C(i))
            initialc(f222,i)=n0(f222);
        end
    end
end
%initialc is the matrix which consists of initial condition of both the
%sectors
factor1=76;factor2=76*76;
species(1:1000,1:20,1:10,1:numel(C))=0;%keep the mer data at each counter
time(1:1000,1:1000,1:numel(C))=0;fusion(2,2)=0;fission(2,2)=0;
fu(1,1)=0.017/(2*(76)^2);fu(1,2)=0.017/((76)^2);%stores the fusion 
                       %frequencies in the presence of microtubule(wild-type mt)
fu(2,1)=0.011/(2*(76)^2);fu(2,2)=0.011/((76)^2);%stores the fusion 
                         %frequencies in the absence of microtubule(no mt)
fi(1,1)=0.018/76;fi(1,2)=2*0.018/76;%stores the fission frequencies
                                   %in the presence of microtubles(wild-type mt)
fi(2,1)=0.027/76;fi(2,2)=2*0.027/76;%stores the fission frequencies
                                    %in the absence of microtubules(no mt)
fu=fu*factor2;fi=fi*factor1;
%factor 2 and factor 1 are multiplied to operate at experimentally observed
%rate constants of fission and fusion
[stochiomatrixchk,random_stochio]=stochiometricmatrix();
%stochiomatrichk is the stochiomatix matrix of all the 201 reactions
toc
for s=1:numel(initialc(1,:))
    n0=initialc(:,s);
     %n0 is the initial condition of a particular sector 
    sectornumber=C(s);
    %sectornumber is the sector in which the 1000 realisations below are
    %carried out
    if sectornumber==1
        mt_ldata=a1w;
        %a1w-length of microtubules in sector 1
    elseif sectornumber==2
        mt_ldata=a2w;
        %a2w-length of microtubules in sector 2
    end
   timemt=0:0.01:1000;
    for i=1:500
        lengthsim=mt_ldata(1,:);
        %lengthsim stores the length of microtubules
        timesim=timemt;
        %timesim stores the time steps of microtubules
        count=1;
        fusion =fu;
        fission=fi;
        cafi(1,1)=0;
        %cafi stores the mitochondria which takes part in fission-fusion
        %reactions in a particular simulation
        tempt(1:1000)=0;temps(1:400,1:20)=0;

        temps(1,:)=n0';
        test=1;
        %test is the counter of while loop
        countin=0;
        while (test>0)
            speciessend=temps(count,:);
            %speciessend is the state vector at every time step
            kb1=fission(2,2);
            kb2=fission(2,1);
            ka1=fusion(2,2);
            ka2=fusion(2,1);
            parameter=[kb1,kb2,ka1,ka2];
            %parameter matrix is the rate constant matrix
            [propensity]=propensitymatrix(speciessend,parameter,stochiomatrixchk);
            %propensitymatrix computes the propensities of all 201 reactions at
       %every time step
            p0=sum(propensity);
            ra=rand(1,2);
            tau=-log(ra(1))/p0;
             %generation of next time instant tau
            tempt(count+1)=tempt(count)+tau;
            %here the current state is updated by execution of the reaction 
        %corresponding to mu
            if(tempt(count+1)>500)
                break;
            end
            [~,chkbint]=histc(tempt(count+1),timesim);
            checklength=lengthsim(chkbint);
            %if the value inside the bracket is 0 then there would be no
            %microtubule and sum would be 1 and chkmt would be 1
            %checklength refers to length of microtubules
            if (checklength~=0)
                ml=checklength;
                countin=countin+1;
                kb1=fission(1,2);
                kb2=fission(1,1);
                ka1=fusion(1,2);
                ka2=fusion(1,1);
                parameter=[kb1,kb2,ka1,ka2];
                %parameter matrix is the rate constant matrix
                [propensity1,cafi]=propensitypresence(speciessend,parameter,stochiomatrixchk);
                %propensitypresence computes the propensities of all 201 reactions at
                %every time step
                if (numel(cafi(1,:))>1)
                    mtp=(lmin.*cafi(:,2));
                %mtp computes the length of the mitochondria participating
                %in fission-fusion reactions
                mtf=find(mtp<ml);
                %mtf find the mitochondria having lengths smaller than the
                %microtubule
                chkelement=isempty(mtf);
                %chkelement checks if mtf is empty or not--if empty then
                %checkelement will be 1
                if (chkelement~=1)
                    propensity1(cafi(mtf,1))=0;
                    %putting propensity of fission for mitochondria having smaller lengths to be zero 
                    p01=sum(propensity1);
                    %recomputing propensity sum in p0
                    if (p01==0)
                        %if no reaction is possible after fission
                        %propensity is zero
                        temps(count+1,:)=temps(count,:);
                        count=count+1;
                        clear('cafi','mtp','mtf');
                        continue;
                    end
                    [~,mu]=histc(ra(2)*p01, [0;cumsum(propensity1(:))]);
                    %mu represents the reaction carried out every time step
                    temps(count+1,:)=temps(count,:)+stochiomatrixchk(mu,:);
                    %here the current state is updated by execution of the reaction 
                    %corresponding to mu
                    count=count+1;
                    %store the updated species and time
                    continue;
                end
                clear('cafi','mtp','mtf');
                %if a shorter microtubule is present and is not in contact
                %with the mitochondria
                [~,mu]=histc(ra(2)*p0, [0;cumsum(propensity(:))]);
                temps(count+1,:)=temps(count,:)+stochiomatrixchk(mu,:);
                a=temps;
                count=count+1;
                continue;
                end
            end
            %if the microtubule length is zero
            [~,mu]=histc(ra(2)*p0, [0;cumsum(propensity(:))]);
            temps(count+1,:)=temps(count,:)+stochiomatrixchk(mu,:);
            a=temps;
            count=count+1;
        end
        if (mod(i,10)==0)
            toc
        end
        species(1:count,:,i,s)=temps(1:count,:);
     %   tempt=tempt';
        time(i,1:count,s)=tempt(1:count);
        clear('temps','tempt');
        clear('lengthsim','timesim');
    end
    clear('mt_ldata','mt_tdata','mt_tmarkdatastart','mt_tmarkdataend');
end
%
if (numel(C)~=1)
    %if only one sector is present
counter=0;
for i=1:numel(time(1,1,:))
    for j=1:numel(time(:,1,1))
        time(j,end+1,i)=0;
        catcht=find(~time(j,2:end,i));
        counter=counter+1;
        time(j,catcht(1),i)=1000000;
        storetime(counter)=time(j,catcht(1),i);
        clear('catcht');
    end
end

T0=0.01:0.01:1000;
%T0 is time matrix at which the above system is investigated
%%interpolation of the state vector at time steps corresponding to T0
for i=1:numel(species(1,1,1,:))
    ffef(1:numel(T0),1:20)=0;
    %In rows of ffef stores ,the number of mitochondria corresponding to each
    %time step in Ta are stored 
    f1s=species(:,:,:,i);
    f1t=time(:,:,i);
    c=0;
    for j=1:numel(f1s(1,1,:))
        f2s=f1s(:,:,j);
        f2t=[0 f1t(j,(find(f1t(j,:))))];
        ta=[f2t(1:end)];
        if (numel(ta)==1)
            continue;
        end
        for k=1:20
            sp=f2s(:,k);
            spi=interp1(ta,sp(1:numel(ta)),T0,'previous');
            spef(1:numel(spi),k)=spi;
        end
        ffef=ffef+spef;
        rr=find(spef(:,9)==1);
        rrnum=numel(rr);
        if mod(j,100)==0
            toc
        end
        c=c+1;
        clear('spef','spi','ta');
    end

    ffef=ffef./c;
    rowlimit=numel(ffef(:,1));
    ffef1(1,:)=f2s(1,:);
    ffef1((2:(rowlimit+1)),:)=ffef;
    ffes(1:numel(ffef1(:,1)),1:numel(ffef1(1,:)),i)=ffef1;
    clear('ffef1');
end
Mwtc4=sum(ffes,3);
%Mwtc4 stores the averaged number of mitochondria present in each fragment
%corresponding to vaious time steps and the averaging is done across the
%realizations
toc
save('cell4sectorwmt','Mwtc4');
else
  %if more than one sector is present     
    for i=1:numel(time(:,1))  
        time(i,end+1)=0;
      catcht=find(~time(i,2:end));
        %counter=counter+1;
        
      
        time(i,catcht(1))=1000000;
        %storetime(counter)=time(j,catcht(1),i);
        clear('catcht');
    end
T0=0.01:0.01:1000;
%for i=1:numel(species(1,1,1,:))
    ffef(1:numel(T0),1:20)=0;
    f1s=species;
    f1t=time;
    c=0;
    for j=1:numel(f1s(1,1,:))
        f2s=f1s(:,:,j);
        f2t=[0 f1t(j,(find(f1t(j,:))))];
        ta=[f2t(1:end)];
        if (numel(ta)==1)
            continue;
        end
        for k=1:20
            sp=f2s(:,k);
            spi=interp1(ta,sp(1:numel(ta)),T0,'previous');
            spef(1:numel(spi),k)=spi;
        end
        ffef=ffef+spef;
        rr=find(spef(:,9)==1);
        rrnum=numel(rr);
        if mod(j,100)==0
            toc
        end
        c=c+1;
        clear('spef','spi','ta');
    end
   
    ffef=ffef./c;
    rowlimit=numel(ffef(:,1));
    ffef1(1,:)=f2s(1,:);
    ffef1((2:(rowlimit+1)),:)=ffef;
    ffes(1:numel(ffef1(:,1)),1:numel(ffef1(1,:)))=ffef1;
    clear('ffef1');

Mwtc4=sum(ffes,3);
%Mwtc4 stores the averaged number of mitochondria present in each fragment
%corresponding to vaious time steps and the averaging is done across the
%realizations
toc
save('cell4sectorwmt','Mwtc4' );
end