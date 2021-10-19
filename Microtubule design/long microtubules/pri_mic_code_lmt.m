function [sl1,sll1,sl2,sll2,sl3,sll3,sl4,sll4,tts1,tts2,tts3,tts4,ttts1,ttts2,ttts3,ttts4,cc1,nume1,nume2,nume3,nume4,cat_l,catfreq1,catfreq2,catfreq3,catfreq4,ttv1,ttv2,ttv3,ttv4,pa1,pa2,pa3,pa4,ls1,ls2,ls3,ls4,simv1,simv2,simv3,simv4]=pri_mic_code_lmt(a,b,k_depolyedge)
%sl1 and sll1 are used to store the lengths of microtubules 
% with sl1 does not include the pasue lengths like
%similarly sl2,sll2 , sl3,sll3 and sl4,sll4 for 2,3 and 4 for four microtubule bundles
%four microtubule bundles are simulated and each bundle has its own time stamps
%tts stores the time stamps of growth for microtubule bundles 1 2 3 4 without the pausetime 
%ttts saves the time matrix with pausetime
%time arrays to store the hydrolysis events of microtubule bundles 1,2,3,4
%cc1 stores the concentration of tubulin in every time step
%pa matrices store the pausetimes of repective microtubules overall 
% combining all realisation pause times
%ttv1,ttv2,ttv3,ttv4 saves time to compute velocities or growth rate in micrometer/s
%catfreq1,catfreq2,catefreq3,catfreq4  are the catastrophe frequencies on every realisation
%simv1,simv2,simv3,simv4 store the shrinkage velocities of microtubule 1,2,3,4
%ls1,ls2,ls3,ls4 stores the lengths of microtubule 1,2,3,4
%a--k_growth , b--k_hydrolysis , K_depolyedge--depolymerisation constant
%nume1,nume2,nume3 and nume4 matrices keep the number of counter experienced at each realisation
tthy1(1:2)=0;
tthy2(1:2)=0;
tthy3(1:2)=0;
tthy4(1:2)=0;

sl1(1)=0;sl2(1)=0;sl3(1)=0;sl4(1)=0;
sll1(1)=0; sll2(1)=0;sll3(1)=0;sll4(1)=0;


tts1(1)=0;tts2(1)=0;tts3(1)=0;tts4(1)=0;
ttts1(1)=0;ttts2(1)=0;ttts3(1)=0;ttts4(1)=0;
ttv1(1)=0;ttv2(1)=0;ttv3(1)=0;ttv4(1)=0;

cat_l(1)=0;

cc1(1)=0;

nume1(1)=0;nume2(1)=0;nume3(1)=0;nume4(1)=0;

pa1(1)=0;pa2(1)=0;pa3(1)=0;pa4(1)=0;
%outer loop works for realisations and inner loop is for simulation
load('longmtld.mat');
lengthcount=numel(initiallongmt);
initiallongmt=round(initiallongmt,1);
tic
for i=1:1000
    cc(1)=3;
    stl=randi([1 lengthcount],4,1);
    stl1=initiallongmt(stl(1));stl2=initiallongmt(stl(2));
    stl3=initiallongmt(stl(3));stl4=initiallongmt(stl(4));
    %cc stores the concentration at every counter
    counter=2;
    %these variables are used to move to the next counter;however these
    %these values are needed to be updated every realisation for storing
    %the repective lengths and time in the correct format
    a11=2;a12=2;a13=2;a14=2;
    af1=1;af2=1;af3=1;af4=1;
    al1=2;al2=2;al3=2;al4=2;
    afl1=1;afl2=1;afl3=1;afl4=1;
    %%
    %concentration of tubuliin 3 micromolar
    con=3;
    %kgrowth and khy are passed as parameters
    
    len1=stl1;len2=stl2;len3=stl3;len4=stl4;
    %len 1,len 2 , len 3 stores the microtubule lengths at every time
    %instant
    t=0;
    %t is a dummy variable
    tc1=1;tc2=1;tc3=1;tc4=1;
    %tc gave the counter to store the lengths for every microtubule
    p1(1)=0;p2(1)=0;p3(1)=0;p4(1)=0;
    %these p's store the pausetime individual realisation
    w1=0;w2=0;w3=0;w4=0; %number of catastrophes
    kt1=0;kt2=0;kt3=0;kt4=0;%time of catastrophes
    %w stores the total number of catastrophes in a realistion
    catfreq1(1)=0;catfreq2(1)=0;catfreq3(1)=0;catfreq4(1)=0;
    %catfreq--catastrophe frequency calcuated by dividing the total number w by the
    %time stored last in tc
    ttgr1(1)=0;t1=0;ttgr2(1)=0;t2=0;ttgr3(1)=0;t3=0;ttgr4(1)=0;t4=0;
    %tgr stores the instantaneous growth time  and is cleaered while ttgr
    %keeps track of time of an entire realisation and t1 stores the growth
    %at every counter
    tv1(1)=0;tvel1=0;tv2(1)=0;tvel2=0; tv3(1)=0;tvel3=0;tv4(1)=0;tvel4=0;
    %tvel is the value of time for velocity and tv is the matrix of those times
    tgr1(1)=0;tgr2(1)=0;tgr3(1)=0;tgr4(1)=0;
    %tgr stores the growth time at every realisation before getting cleared
    %and handling the matrix over to tts
    thy1(1)=0;thy2(1)=0;thy3(1)=0;thy4(1)=0;
    %thy store the hydrolysis time at eery realisation
    l1_de(1)=0;l2_de(1)=0;l3_de(1)=0;l4_de(1)=0;
    %stores the depolymerisation lengths l_des
    t1_de(1)=0;t2_de(1)=0;t3_de(1)=0;t4_de(1)=0;
    %stores the depolymrisation times t_des
    pos1_de(1)=0;pos2_de(1)=0;pos3_de(1)=0;pos4_de(1)=0;
    %stores the position of depolymerisation--meaning it is a continuous
    %length array so the position where to add the depolymerised lengths is
    %given by this pos_de
    ll1(1,:)=stl1;ll2(1,:)=stl2;ll3(1,:)=stl3;ll4(1,:)=stl4;
    %ll gives the lengths without depolymeried lengths series like 1 2 3
    %then 0 not 3 2 1 0 ;
    lll1(1,:)=stl1;lll2(1,:)=stl2;lll3(1,:)=stl3;lll4(1,:)=stl4;
    chkcatfreq1=0; chkcatfreq2=0; chkcatfreq3=0; chkcatfreq4=0;
    %however %lll gives 0 1 2  3 3 2 1 0
    %1- GTP %2 -GDP
    currenttubstate1=1;currenttubstate2=1;currenttubstate3=1;currenttubstate4=1;
    %tipstates are initiated with 1
      %number of GTP molecules
    nGTPmt1=len1*10;nGTPmt2=len2*10;nGTPmt3=len3*10;nGTPmt4=len4*10;
    countshort1=0;countshort2=0;countshort3=0;countshort4=0;
    simv1(1)=0;simv2(1)=0;simv3(1)=0;simv4(1)=0;
    ls1(1)=0;ls2(1)=0;ls3(1)=0;ls4(1)=0;
    if i==1
        sl1(1)=stl1;sl2(1)=stl2;sl3(1)=stl3;sl4(1)=stl4;
        sll1(1)=stl1; sll2(1)=stl2;sll3(1)=stl3;sll4(1)=stl4;
    end
    while (counter>0)
        %microtubules are chosen to grow randomly from the nucleus
        %overall 4 microtubule bundles are present hence 1 4 randi is used
        p0=0;propensity(3)=0;
        d=randi([1 4]);
        if(d==1)%microtubule 1
            kgr=a;khy=b;kde=k_depolyedge;
            if (round(len1,1)<5.8)
                [propensity]= propensitymatrix_mt(currenttubstate1,con,nGTPmt1,kgr,khy,kde);
                p0=sum(propensity);
                ra=rand(1,2);
                tau1=-log(ra(1))/p0;
                tc1=tc1+1;
                ttgr1(tc1)=ttgr1(tc1-1)+tau1;
                [~,mu]=histc(ra(2)*p0,[0;cumsum(propensity(:))]);
                if mu ==1  %%growth event
                    len1=len1+0.1;
                    con=con-0.00614;
                    nGTPmt1=nGTPmt1+1;
                    pausetime1=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==2 %%hydrolysis event
                    nGTPmt1=nGTPmt1-1;
                    pausetime1=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==3 %% depolymerisation
                    pausetime1=ttgr1(tc1-1)-ttgr1(tc1-2);
                    p1(tc1)=pausetime1;
                    w1=w1+1;
                    kt1=kt1+ttgr1(tc1-1);
                    con=con+len1*10*0.00614;                   
                   % l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                    tsend=ttgr1(tc1-1);
                    [~,tmat,lmat]=depolymerisation_GTPseed(kde,tsend,len1);
                    countshort1=countshort1+1;
                    ls1(countshort1)=ll1(tc1-1);
                    simv1(countshort1)=(ll1(tc1-1)-0)/(tmat(end)-ttgr1(tc1-1));
                    ll1=[ll1(1:(tc1-1)) lmat];
                    ttgr1=[ttgr1(1:(tc1-1)) tmat];
                    tc1=numel(ttgr1);
                    currenttubstate1=1;
                    chkcatfreq1=1;
                    len1=0;
                    cc(counter)=con;

                    o=max(ttgr1);p=max(ttgr2);q=max(ttgr3);efour=max(ttgr4);
                    clear('l_mat','t_mat','pos_m','tde','l1');
                    if(o >1000 && p>1000 && q>1000 && efour>1000)
                        break;
                    end
                    counter=counter+1;
                    continue
                end
                p1(tc1)=pausetime1;
                [l1_de]=[l1_de l_mat];
                [t1_de]=[t1_de t_mat];
                [pos1_de]=[pos1_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll1(tc1)=len1;
                %check tipstate
                if nGTPmt1>0
                    currenttubstate1=1;
                elseif nGTPmt1==0
                    currenttubstate1=2;
                end
                %depolymerisation at edge of the cell
            elseif (round(len1,1)==5.8)
                pausetime1=pt_cal(nGTPmt1,khy,chkcatfreq1);
                l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                tsend=ttgr1(tc1)+pausetime1;
                [tadd,tde,l1]=depolymerisation_GTPseed(kde,tsend,len1);
                l_mat=[l_mat  len1 len1 l1];
                t_mat=[t_mat ttgr1(tc1) tsend tde];
                pos_m=[pos_m tc1];
                con=con+len1*10*0.00614;
                countshort1=countshort1+1;
                ls1(countshort1)=5.8;
                simv1(countshort1)=5.8/(tadd-tsend);
                len1=0;
                tc1=tc1+1;
                ttgr1(tc1)=tadd;
                w1=w1+1;
                kt1=kt1+ttgr1(tc1-1)+pausetime1;
                currenttubstate1=1;
                p1(tc1)=pausetime1;
                [l1_de]=[l1_de l_mat];
                [t1_de]=[t1_de t_mat];
                [pos1_de]=[pos1_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll1(tc1)=len1;
                chkcatfreq1=1;
                nGTPmt1=0;
            end
        elseif (d==2)%microtubule 2           
            kgr=a;khy=b;kde=k_depolyedge;
            if(round(len2,1)<5.8)
                [propensity]= propensitymatrix_mt(currenttubstate2,con,nGTPmt2,kgr,khy,kde);
                p0=sum(propensity);
                ra=rand(1,2);
                tau2=-log(ra(1))/p0;
                tc2=tc2+1;
                ttgr2(tc2)=ttgr2(tc2-1)+tau2;
                [~,mu]=histc(ra(2)*p0,[0;cumsum(propensity(:))]);
                if mu ==1  %%growth event
                    len2=len2+0.1;
                    con=con-0.00669;
                    nGTPmt2=nGTPmt2+1;
                    pausetime2=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==2
                    nGTPmt2=nGTPmt2-1;
                    pausetime2=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==3
                    pausetime2=ttgr2(tc2-1)-ttgr2(tc2-2);
                    p2(tc2)=pausetime2;
                    w2=w2+1;
                    kt2=kt2+ttgr2(tc2-1);
                    con=con+len2*10*0.00614;                   
                   % l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                    tsend=ttgr2(tc2-1);
                    [~,tmat,lmat]=depolymerisation_GTPseed(kde,tsend,len2);
                    countshort2=countshort2+1;
                    ls2(countshort2)=ll2(tc2-1);
                    simv2(countshort2)=(ll2(tc2-1)-0)/(tmat(end)-ttgr2(tc2-1));
                    ll2=[ll2(1:(tc2-1)) lmat];
                    ttgr2=[ttgr2(1:(tc2-1)) tmat];
                    tc2=numel(ttgr2);
                    currenttubstate2=1;
                    len2=0;
                    chkcatfreq2=1;
                    cc(counter)=con;
                    o=max(ttgr1);p=max(ttgr2);q=max(ttgr3);efour=max(ttgr4);
                    clear('l_mat','t_mat','pos_m','tde','l1');
                    if(o >1000 && p>1000 && q>1000 && efour>1000)
                        break;
                    end
                    counter=counter+1;
                    continue
                end
                p2(tc2)=pausetime2;
                [l2_de]=[l2_de l_mat];
                [t2_de]=[t2_de t_mat];
                [pos2_de]=[pos2_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll2(tc2)=len2;
                if nGTPmt2>0
                    currenttubstate2=1;
                elseif nGTPmt2==0
                    currenttubstate2=2;
                end
            elseif (round(len2,1)==5.8)
                pausetime2=pt_cal(nGTPmt2,khy,chkcatfreq2);
                l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                tsend=ttgr2(tc2)+pausetime2;
                [tadd,tde,l1]=depolymerisation_GTPseed(kde,tsend,len2);
                l_mat=[l_mat  len2 len2  l1];
                t_mat=[t_mat ttgr2(tc2) tsend tde];
                pos_m=[pos_m tc2];
                con=con+len2*10*0.00614;
                countshort2=countshort2+1;
                ls2(countshort2)=5.8;
                simv2(countshort2)=5.8/(tadd-tsend);
                len2=0;
                tc2=tc2+1;
                ttgr2(tc2)=tadd;
                w2=w2+1;
                kt2=kt2+ttgr2(tc2-1)+pausetime2;
                currenttubstate2=1;
                p2(tc2)=pausetime2;
                [l2_de]=[l2_de l_mat];
                [t2_de]=[t2_de t_mat];
                [pos2_de]=[pos2_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll2(tc2)=len2;
                chkcatfreq2=1;
                nGTPmt2=0;
            end
        elseif (d==3) %microtubule 3
            kgr=a;khy=b;kde=k_depolyedge;
            if(round(len3,1)<5.8)
                [propensity]= propensitymatrix_mt(currenttubstate3,con,nGTPmt3,kgr,khy,kde);
                p0=sum(propensity);
                ra=rand(1,2);
                tau3=-log(ra(1))/p0;
                tc3=tc3+1;
                ttgr3(tc3)=ttgr3(tc3-1)+tau3;
                [~,mu]=histc(ra(2)*p0,[0;cumsum(propensity(:))]);
                if mu ==1  %%growth event
                    len3=len3+0.1;
                    con=con-0.00614;
                    nGTPmt3=nGTPmt3+1;
                    pausetime3=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==2
                    nGTPmt3=nGTPmt3-1;
                    pausetime3=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==3
                    pausetime3=ttgr3(tc3-1)-ttgr3(tc3-2);
                    p3(tc3)=pausetime3;
                    w3=w3+1;
                    kt3=kt3+ttgr3(tc3-1);
                    con=con+len3*10*0.00614;
                    tsend=ttgr3(tc3-1);
                    [~,tmat,lmat]=depolymerisation_GTPseed(kde,tsend,len3);
                    countshort3=countshort3+1;
                    ls3(countshort3)=ll3(tc3-1);
                    simv3(countshort3)=(ll3(tc3-1)-0)/(tmat(end)-ttgr3(tc3-1));
                    ll3=[ll3(1:(tc3-1)) lmat];
                    ttgr3=[ttgr3(1:(tc3-1)) tmat];
                    tc3=numel(ttgr3);
                    currenttubstate3=1;
                    chkcatfreq3=1;
                    len3=0;
                    cc(counter)=con;
                    o=max(ttgr1);p=max(ttgr2);q=max(ttgr3);efour=max(ttgr4);
                    clear('l_mat','t_mat','pos_m','tde','l1');
                    if(o >1000 && p>1000 && q>1000 && efour>1000)
                        break;
                    end
                    counter=counter+1;
                    continue
                end
                p3(tc3)=pausetime3;
                [l3_de]=[l3_de l_mat];
                [t3_de]=[t3_de t_mat];
                [pos3_de]=[pos3_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll3(tc3)=len3;
                if nGTPmt3>0
                    currenttubstate3=1;
                elseif nGTPmt3==0
                    currenttubstate3=2;
                end
            elseif(round(len3,1)==5.8)
                pausetime3=pt_cal(nGTPmt3,khy,chkcatfreq3);
                l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                tsend=ttgr3(tc3)+pausetime3;
                [tadd,tde,l1]=depolymerisation_GTPseed(kde,tsend,len3);
                l_mat=[l_mat  len3 len3 l1];
                t_mat=[t_mat ttgr3(tc3) tsend tde];
                pos_m=[pos_m tc3];
                con=con+len3*10*0.00614;
                countshort3=countshort3+1;
                ls3(countshort3)=5.8;
                simv3(countshort3)=5.8/(tadd-tsend);
                len3=0;
                tc3=tc3+1;
                ttgr3(tc3)=tadd;
                w3=w3+1;
                kt3=kt3+ttgr3(tc3-1)+pausetime3;
                currenttubstate3=1;
                p3(tc3)=pausetime3;
                [l3_de]=[l3_de l_mat];
                [t3_de]=[t3_de t_mat];
                [pos3_de]=[pos3_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll3(tc3)=len3;
                chkcatfreq3=1;
                nGTPmt3=0;
            end
        elseif (d==4)%microtubule 4
            kgr=a;khy=b;kde=k_depolyedge;
            if(round(len4,1)<5.8)
                [propensity]= propensitymatrix_mt(currenttubstate4,con,nGTPmt4,kgr,khy,kde);
                p0=sum(propensity);
                ra=rand(1,2);
                tau4=-log(ra(1))/p0;
                tc4=tc4+1;
                ttgr4(tc4)=ttgr4(tc4-1)+tau4;
                [~,mu]=histc(ra(2)*p0,[0;cumsum(propensity(:))]);
                if mu ==1  %%growth event
                    len4=len4+0.1;
                    con=con-0.00614;
                    nGTPmt4=nGTPmt4+1;
                    pausetime4=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==2
                    nGTPmt4=nGTPmt4-1;
                    pausetime4=0;
                    l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                elseif mu ==3
                    pausetime4=ttgr4(tc4-1)-ttgr4(tc4-2);
                    p4(tc4)=pausetime4;
                    w4=w4+1;
                    kt4=kt4+ttgr4(tc4-1);
                    con=con+len4*10*0.00614;                   
                    tsend=ttgr4(tc4-1);
                    [~,tmat,lmat]=depolymerisation_GTPseed(kde,tsend,len4);
                    countshort4=countshort4+1;
                    ls4(countshort4)=ll4(tc4-1);
                    simv4(countshort4)=(ll4(tc4-1)-0)/(tmat(end)-ttgr4(tc4-1));
                    ll4=[ll4(1:(tc4-1)) lmat];
                    ttgr4=[ttgr4(1:(tc4-1)) tmat];
                    tc4=numel(ttgr4);
                    currenttubstate4=1;
                    chkcatfreq4=1;
                    len4=0;
                    cc(counter)=con;
                    o=max(ttgr1);p=max(ttgr2);q=max(ttgr3);efour=max(ttgr4);
                    clear('l_mat','t_mat','pos_m','tde','l1');
                    if(o >1000 && p>1000 && q>1000 && efour>1000)
                        break;
                    end
                    counter=counter+1;
                    continue
                end
                p4(tc4)=pausetime4;
                [l4_de]=[l4_de l_mat];
                [t4_de]=[t4_de t_mat];
                [pos4_de]=[pos4_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll4(tc4)=len4;
                if nGTPmt4>0
                    currenttubstate4=1;
                elseif nGTPmt4==0
                    currenttubstate4=2;
                end
            elseif(round(len4,1)==5.8)
                pausetime4=pt_cal(nGTPmt4,khy,chkcatfreq4);
                l_mat(1)=0;t_mat(1)=0;pos_m(1)=0;
                tsend=ttgr4(tc4)+pausetime4;
                [tadd,tde,l1]=depolymerisation_GTPseed(kde,tsend,len4);
                l_mat=[l_mat  len4  len4  l1];
                t_mat=[t_mat ttgr4(tc4) tsend tde];
                pos_m=[pos_m tc4];
                con=con+len4*10*0.00614;
                countshort4=countshort4+1;
                ls4(countshort4)=5.8;
                simv4(countshort4)=5.8/(tadd-tsend);
                len4=0;
                tc4=tc4+1;
                ttgr4(tc4)=tadd;
                w4=w4+1;
                kt4=kt4+ttgr4(tc4-1)+pausetime4;
                currenttubstate4=1;
                p4(tc4)=pausetime4;
                [l4_de]=[l4_de l_mat];
                [t4_de]=[t4_de t_mat];
                [pos4_de]=[pos4_de pos_m];
                clear('l_mat','t_mat','pos_m','tde','l1');
                ll4(tc4)=len4;
                chkcatfreq4=1;
                nGTPmt4=0;
            end
        end
        cc(counter)=con;
        o=max(ttgr1);p=max(ttgr2);q=max(ttgr3);efour=max(ttgr4);
        if(o >1000 && p>1000 && q>1000 && efour>1000)
            break;
        end
        %o, p and q stores the maximumtine in a counter ideally 228 s is seen
        %in experiments we are going for 4000/2000 seconds while linking
        %mitochondria to microtubule for being on the safe limit we keep it to 7000
        counter=counter+1;
    end
    lll1=ll1;lll2=ll2;lll3=ll3;lll4=ll4;
    tgr1=ttgr1;tgr2=ttgr2;tgr3=ttgr3;tgr4=ttgr4;
    %complete gives polymerised +depolymerised lengths with the pause times
    %included in growth time arrays
    [lll1,tgr1]=complete(lll1,tgr1,l1_de,t1_de,pos1_de);
    [lll2,tgr2]=complete(lll2,tgr2,l2_de,t2_de,pos2_de);
    [lll3,tgr3]=complete(lll3,tgr3,l3_de,t3_de,pos3_de);
    [lll4,tgr4]=complete(lll4,tgr4,l4_de,t4_de,pos4_de);
    %numel will give you the limits one realisation has gone
    nume1(i)=tc1;nume2(i)=tc2;nume3(i)=tc3;nume4(i)=tc4;
    catfreq1(i)=w1/kt1;catfreq2(i)=w2/kt2;catfreq3(i)=w3/kt3;catfreq4(i)=w4/kt4;
    %sl stores the length without depolymerisation lengths
    [sl1]=len_stor(ll1,sl1);[sl2]=len_stor(ll2,sl2);[sl3]=len_stor(ll3,sl3);[sl4]=len_stor(ll4,sl4);
    %sll stores the length with depolymerisation lengths
    [sll1]=len_stor(lll1,sll1);[sll2]=len_stor(lll2,sll2);[sll3]=len_stor(lll3,sll3);[sll4]=len_stor(lll4,sll4);
    %ttgr stores the growth time without depolymerisation times
    [tts1]=time_stor(ttgr1,tts1);[tts2]=time_stor(ttgr2,tts2);[tts3]=time_stor(ttgr3,tts3);[tts4]=time_stor(ttgr4,tts4);
    %tts stores the growth time without depolymeriation times
    [ttts1]=time_stor(tgr1,ttts1);[ttts2]=time_stor(tgr2,ttts2);[ttts3]=time_stor(tgr3,ttts3);[ttts4]=time_stor(tgr4,ttts4);
    %ttts stores the growth time with depolymerisation time
    %[ttv1]=time_stor(tv1,ttv1);[ttv2]=time_stor(tv2,ttv2);[ttv3]=time_stor(tv3,ttv3);[ttv4]=time_stor(tv4,ttv4);
    %ttv stores the velocity time at respective counters
    %[tthy1]=time_stor(thy1,tthy1);[tthy2]=time_stor(thy2,tthy2);[tthy3]=time_stor(thy3,tthy3);[tthy4]=time_stor(thy4,tthy4);
    %tthy stores the hydrolysis time
    [cc1]=conc(cc,cc1);
    %cc1 stores the concentration at respective realisations
    [pa1]=paus_stor(p1,pa1);[pa2]=paus_stor(p2,pa2);[pa3]=paus_stor(p3,pa3);[pa4]=paus_stor(p4,pa4);
    %pa stores the pause time matrix at the end of each realisation
    clear('ll1','lll1','lll2','lll3','lll4');
    clear('ll2');
    clear('ll3');
    clear('ll4');
    clear('ttgr1','tgr1','tgr2','tgr3','tgr4','ttgr4');clear('ttgr2');clear('ttgr3');
    clear('cc');
    clear('catas');
    %  clear('tv1');clear('tv2');clear('tv3','tv4');
    % clear('thy1','thy2','thy3','thy4');
    clear('l1_de','l2_de','l3_de','l4_de');clear('t1_de','t2_de','t3_de','t4_de');
    clear('pos1_de','pos2_de','pos3_de','pos4_de');
end
end

