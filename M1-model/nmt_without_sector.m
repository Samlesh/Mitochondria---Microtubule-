function [ffes]=nmt_without_sector(principalmat)
%principalmat -- the matrix which consists of the number of mitochondria
            %present in each fragement at various time steps in experiments

initial=principalmat;%stores the time evolution of cell 1
n0=initial(:,1);%stores the initial condition
factor1=5;factor2=5*5;
species(1:1000,1:20,1:10)=0;%keep the mer data at each counter
time(1:1000,1:1000)=0;fusion(2,2)=0;fission(2,2)=0;
fu(1,1)=0.017/(2*(5)^2);fu(1,2)=0.017/((5)^2);%stores the fusion frequencies in the presence of microtubule(wild mt)
fu(2,1)=0.011/(2*(5)^2);fu(2,2)=0.011/((5)^2);%stores the fusion 
                         %frequencies in the absence of microtubule(no mt)
fi(1,1)=0.018/5;fi(1,2)=2*0.018/5;%stores the fission frquencies in the presence of microtubles(wmt)
fi(2,1)=0.027/5;fi(2,2)=2*0.027/5;%stores the fission frequencies
                                    %in the absence of microtubules(no mt)
fu=fu*factor2;fi=fi*factor1;
%factor 2 and factor 1 are multiplied to operate at experimentally observed
%rate constants of fission and fuion
[stochiomatrixchk,random_stochio]=stochiometricmatrix();
%stochiomatrichk is the stochiomatix matrix of all the 201 reactions
for i=1:1000
    count=1;
    fusion =fu;
    fission=fi;
    tempt=zeros(2,1);temps=zeros(2,20);
    temps(1,:)=n0';
   %n0' is the initial condition in a cell at t=0 amd is stored in the 
   %first row of species matrix temps    
    test=1;stop=0;
    %test is the counter of while loop
    while (test>0)
        stop=stop+1;
        speciessend=temps(count,:);
        %speciessend is the state vector at every time step
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
        if(tempt(count+1)>1000)
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
    species(1:numel(temps(:,1)),:,i)=temps;
    %species stores the number of mitochondria in all the length fragments
    %at respective timestep for ith realisation
    time(i,1:numel(tempt'))=tempt';
    %time stores the time steps for all the realisations 
    clear('temps','tempt');
    
end

Ta=0.1:0.1:1000;
%Ta is time matrix at which the above system is investigated
ffef(1:numel(Ta),1:20)=0;
%In rows of ffef stores ,the number of mitochondria corresponding to each
%time step in Ta are stored 
c=0;
%%interpolation of the state vector at time steps corresponding to Ta
for j=1:numel(species(1,1,:))
    f2s=species(:,:,j);
    f2t=[0 time(j,(find(time(j,:))))];
    ta=[f2t(1:end-1)];
    if (numel(ta)==0)
        continue;
    end
    for k=1:20
        sp=f2s(:,k);
        spi=interp1(ta,sp(1:numel(ta)),Ta,'previous');
        spef(1:numel(spi),k)=spi;
    end
    ffef=ffef+spef;
    %
    if mod(j,100)==0
        toc
    end
    c=c+1;
    clear('spef','spi','ta','f2s','f2t','sp');
end
ffef=ffef./c;
ffes(1:numel(ffef(:,1)),1:numel(ffef(1,:)))=ffef;
%ffes stores the averaged number of mitochondria present in each fragment
%corresponding to vaious time steps and the averaging is done across the
%realizations
end