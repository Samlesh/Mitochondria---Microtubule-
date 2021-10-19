function [mitochondria]=lmt_with_sector(principalmat)
%principalmat -- the matrix which consists of the number of mitochondria
            %present in each fragement at various time steps in experiments
initial=principalmat;%stores the time evolution of respective cell
n0=initial(:,1);%stores the initial condition
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
factor1=3;factor2=3*3;
species(1:1000,1:20,1:10,1:numel(C))=0;%keep the mer data at each counter
time(1:10,1:1000,1:numel(C))=0;fusion(2,2)=0;fission(2,2)=0;
fu(1,1)=0.011/(2*(3)^2);fu(1,2)=0.011/((3)^2);%stores the fusion 
                       %frequencies in the presence of microtubule(long mt)
fu(2,1)=0.011/(2*(3)^2);fu(2,2)=0.011/((3)^2);%stores the fusion 
                         %frequencies in the absence of microtubule(no mt)
fi(1,1)=0.006*0.01/3;fi(1,2)=2*0.006*0.01/3;%stores the fission frequencies
                                   %in the presence of microtubles(long mt)
fi(2,1)=0.027/3;fi(2,2)=2*0.027/3;%stores the fission frequencies
                                    %in the absence of microtubules(no mt)
fu=fu*factor2;fi=fi*factor1;
%factor 2 and factor 1 are multiplied to operate at experimentally observed
%rate constants of fission and fuion
[stochiomatrixchk,random_stochio]=stochiometricmatrix();
%stochiomatrichk is the stochiomatix matrix of all the 201 reactions
for s=1:numel(initialc(1,:))
    n0=initialc(:,s);
    %n0 is the initial condition of a particular sector 
    sectornumber=C(s);
    %sectornumber is the sector in which the 1000 realisations below are
    %carried out
for i=1:1000
    count=1;
    fusion =fu;
    fission=fi;
    tempt=zeros(2,1);temps=zeros(2,20);
    temps(1,:)=n0';
    test=1;stop=0;
    %test is the counter of while loop
    while (test>0)
        stop=stop+1;
        speciessend=temps(count,:);
        %speciessend is the state vector at every time step
        kb1=fission(1,2);
        kb2=fission(1,1);
        ka1=fusion(1,2);
        ka2=fusion(1,1);
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
        if(tempt(count+1)>900)
            break;
        end
        [~,mu]=histc(ra(2)*p0, [0;cumsum(propensity(:))]);
        %mu represents the reaction carried out every time step
        temps(count+1,:)=temps(count,:)+stochiomatrixchk(mu,:);
        %here the current state is updated by execution of the reaction 
        %corresponding to mu
        count=count+1;
        
    end
    
    if (mod(i,100)==0)
        toc
    end
      species(1:count,:,i,s)=temps(1:count,:);
      time(i,1:count,s)=tempt(1:count);
      clear('temps','tempt');  
end
end

for i=1:numel(time(1,1,:))
    counter=0;
    for j=1:numel(time(:,1,1))
        catcht=find(~time(j,2:end,i));
        counter=counter+1;
        time(j,catcht(1),i)=1000;
        storetime(counter)=time(j,catcht(1),i);
        clear('catcht');
    end
end
%storetime stores the counter at which last time step of aparticular
%simulation is stored

if (numel(C)~=1)
    %if only one sector is present
T0=0.1:0.1:600;
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
        f2t=[0 f1t(j,(find(f1t(j,:)))) ];
        ta=[f2t(1:end)];
       if (numel(ta)==0 || numel(ta)==1 ||  ta(2)==1000 )
            continue;
        end
        for k=1:20
            sp=f2s(:,k);
            spi=interp1(ta,sp(1:numel(ta)),T0,'previous');
            spef(1:numel(spi),k)=spi;
            clear('sp','spi')
        end
        ffef=ffef+spef;
        rr=find(spef(:,9)==1);
        rrnum=numel(rr);
        if mod(j,100)==0
            toc
        end
        c=c+1;
        clear('spef','ta','f2t');
    end

    ffef=ffef./c;
    rowlimit=numel(ffef(:,1));
    ffef1(1,:)=f2s(1,:);
    ffef1((2:(rowlimit+1)),:)=ffef;
    ffes(1:numel(ffef1(:,1)),1:numel(ffef1(1,:)),i)=ffef1;
    clear('ffef1','ffef');
end
Mltc=sum(ffes,3);
mitochondria=Mltc;
%mitochondria stores the averaged number of mitochondria present in each fragment
%corresponding to vaious time steps and the averaging is done across the
%realizations
else
%if more thanone sector is present    
T0=0.1:0.1:600;
    ffef(1:numel(T0),1:20)=0;
    f1s=species;
    f1t=time;
    c=0;
    for j=1:numel(f1s(1,1,:))
        f2s=f1s(:,:,j);
        f2t=[0 f1t(j,(find(f1t(j,:))))];
        
        ta=[f2t(1:end)];
        
       if (numel(ta)==0 || numel(ta)==1 || ta(2)==1000)
            continue;
        end
        for k=1:20
            sp=f2s(:,k);
            spi=interp1(ta,sp(1:numel(ta)),T0,'previous');
            spef(1:numel(spi),k)=spi;
            clear('sp','spi')
        end
        ffef=ffef+spef;
        rr=find(spef(:,9)==1);
        rrnum=numel(rr);
        if mod(j,100)==0
            toc
        end
        c=c+1;
          clear('spef','ta','f2t');
    end

    ffef=ffef./c;
    rowlimit=numel(ffef(:,1));
    ffef1(1,:)=f2s(1,:);
    ffef1((2:(rowlimit+1)),:)=ffef;
    ffes(1:numel(ffef1(:,1)),1:numel(ffef1(1,:)))=ffef1;
    clear('ffef1');
    mitochondria=ffes;
%mitochondria stores the averaged number of mitochondria present in each fragment
%corresponding to vaious time steps and the averaging is done across the
%realizations

end
end