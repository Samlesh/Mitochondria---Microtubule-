function[new2]=binning(g,lmin,lint)
new2(20,20)=0;bincounts(20)=0;
for i=1:20
    bincounts(i)=lmin;
    lmin=lmin+lint;
end
for i=1:20
    mermatrix(1:20)=0;
    e=g{1,i}(:,1);
    f=numel(e);
    for j=1:f
        if (e(j)>0.6 && e(j)<11.6)
        [~,ind]=histc(e(j),bincounts(:));
        mermatrix(ind)=mermatrix(ind)+1;
        end
    end
    new2(i,:)=mermatrix(:);
end
end


% a=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31 a32 a33 a34 a35 a36 a37 a38 a39 a40 a41 a42];
% for count=1:numel(a)
%     new2(i,count)=a(count);
% end
%summatrix(i)=sum(new(i,:));
%         if(e(j)>=lmin && e(j)<lmin+t)
%             a1=a1+1;
%         elseif(e(j)>=lmin+t &The printing of smart ID cards for other members whose name is not listed in attached sheet is in progress. Once the printed cards will be received by us, it will be informed and delivered to respective departments only.
%& e(j)<lmin+2*t)
%             a2=a2+1;
%         elseif(e(j)>=lmin+2*t && e(j)<lmin+3*t)
%             a3=a3+1;
%         elseif(e(j)>=lmin+3*t && e(j)<lmin+4*t)
%             a4=a4+1;
%         elseif(e(j)>=lmin+4*t && e(j)<lmin+5*t)
%             a5=a5+1;
%         elseif(e(j)>=lmin+5*t && e(j)<lmin+6*t)
%             a6=a6+1;
%         elseif(e(j)>=lmin+6*t && e(j)<lmin+7*t)
%             a7=a7+1;
%         elseif(e(j)>=lmin+7*t && e(j)<lmin+8*t)
%             a8=a8+1;
%         elseif(e(j)>=lmin+8*t && e(j)<lmin+9*t)
%             a9=a9+1;
%         elseif(e(j)>=lmin+9*t && e(j)<lmin+10*t)
%             a10=a10+1;
%         elseif(e(j)>=lmin+10*t && e(j)<lmin+11*t)
%             a11=a11+1;
%         elseif(e(j)>=lmin+11*t && e(j)<lmin+12*t)
%             a12=a12+1;
%         elseif(e(j)>=lmin+12*t && e(j)<lmin13*t)
%             a13=a13+1;
%         elseif(e(j)>=lmin+13*t && e(j)<lmin+14*t)
%             a14=a14+1;
%         elseif(e(j)>=lmin+14*t && e(j)<lmin+15*t)
%             a15=a15+1;
%         elseif(e(j)>=lmin+15*t && e(j)<lmin+16*t)
%             a16=a16+1;
%         elseif(e(j)>=lmin+16*t && e(j)<lmin+17*t)
%             a17=a17+1;
%         elseif(e(j)>=lmin+17*t && e(j)<lmin+18*t)
%             a18=a18+1;
%         elseif(e(j)>=lmin+18*t && e(j)<lmin+19*t)
%             a19=a19+1;
%      a=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31 a32 a33 a34 a35 a36 a37 a38 a39 a40 a41 a42];
% for count=1:numel(a)
%     new2(i,count)=a(count);
% end
%summatrix(i)=sum(new(i,:));    elseif(e(j)>=lmin+19*t && e(j)<lmin+20*t)
%             a20=a20+1;
%         elseif(e(j)>=lmin+20*t && e(j)<lmin+21*t)
%             a21=a21+1;
%         elseif(e(j)>=lmin+21*t && e(j)<lmin+22*t)
%             a22=a22+1;  
%         elseif(e(j)>=lmin+22*t && e(j)<lmin+23*t)
%             a23=a23+1;
%         elseif(e(j)>=lmin+23*t && e(j)<lmin+24*t)
%             a24=a24+1;
%         elseif(e(j)>=lmin+24*t && e(j)<lmin+25*t)
%             a25=a25+1;
%         elseif(e(j)>=lmin+25*t && e(j)<lmin+26*t)
%             a26=a26+1;
%         elseif(e(j)>=lmin+26*t && e(j)<lmin+27*t)
%             a27=a27+1;
%         elseif(e(j)>=lmin+27*t && e(j)<lmin+28*t)
%             a28=a28+1;
%         elseif(e(j)>=lmin+28*t && e(j)<lmin+29*t)
%             a29=a29+1;
%         elseif(e(j)>=lmin+29*t && e(j)<lmin+30*t)
%             a30=a30+1;
%         elseif(e(j)>=lmin+30*t && e(j)<lmin+31*t)
%             a31=a31+1;
%         elseif(e(j)>=lmin+31*t && e(j)<lmin+32*t)
%             a32=a32+1;
%         elseif(e(j)>=lmin+32*t && e(j)<lmin33*t)
%             a33=a33+1;
%         elseif(e(j)>=lmin+33*t && e(j)<lmin+34*t)
%             a34=a34+1;
%         elseif(e(j)>=lmin+34*t && e(j)<lmin+35*t)
%             a35=a35+1;
%         elseif(e(j)>=lmin+35*t && e(j)<lmin+36*t)
%             a36=a36+1;
%         elseif(e(j)>=lmin+36*t && e(j)<lmin+37*t)
%             a37=a37+1;
%         elseif(e(j)>=lmin+37*t && e(j)<lmin+38*t)
%             a38=a38+1;
%         elseif(e(j)>=lmin+38*t && e(j)<lmin+39*t)
%             a39=a39+1;
%         elseif(e(j)>=lmin+39*t && e(j)<lmin+40*t)
%             a40=a40+1;
%         elseif(e(j)>=lmin+40*t && e(j)<lmin+41*t)
%             a41=a41+1;
%         elseif(e(j)>=lmin+41*t && e(j)<lmin+42*t)
%             a42=a42+1;        
%         elseif(e(j)>=lmin+42*t && e(j)<lmin+43*t)
%             a43=a43+1;
%         elseif(e(j)>=lmin+43*t && e(j)<lmin+44*t)
%             a44=a44+1;
%         elseif(e(j)>=lmin+44*t && e(j)<lmin+31*t)
%             a31=a31+1;
%         elseif(e(j)>=lmin+31*t && e(j)<lmin+32*t)
%             a32=a32+1;
%         elseif(e(j)>=lmin+32*t && e(j)<lmin33*t)
%             a33=a33+1;
%         elseif(e(j)>=lmin+33*t && e(j)<lmin+34*t)
%             a34=a34+1;
%         elseif(e(j)>=lmin+34*t && e(j)<lmin+35*t)
%             a35=a35+1;
%         elseif(e(j)>=lmin+35*t && e(j)<lmin+36*t)
%             a36=a36+1;
%         elseif(e(j)>=lmin+36*t && e(j)<lmin+37*t)
%             a37=a37+1;
%         elseif(e(j)>=lmin+37*t && e(j)<lmin+38*t)
%             a38=a38+1;
%         elseif(e(j)>=lmin+38*t && e(j)<lmin+39*t)
%             a39=a39+1;
%         elseif(e(j)>=lmin+39*t && e(j)<lmin+40*t)
%             a40=a40+1;
%         elseif(e(j)>=lmin+40*t && e(j)<lmin+41*t)
%             a41=a41+1;
%         elseif(e(j)>=lmin+41*t && e(j)<lmin+42*t)
%             a42=a42+1; 