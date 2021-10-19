function [set]=p_chk(sl1,pal,nume)
cach_z(2)=0;
cal(2)=0;
cap(2)=0;
e=0;
templ(1)=0;
c=0;d=0;f=0;
d(1)=0;
for i=1:numel(nume)
    % c=c+nume(i);
    if(i~=1)
        f=f+1;
        d(f)=d(f-1)+nume(i-1);
    else
        f=f+1;
        d(f)=1;
    end
end
for i=1:numel(sl1)
    m=0;
    m=sum(d==i);
    if (round(sl1(i),1)==0 && m==0)
        e=e+1;
        cal(e)=sl1(i-1);
        cap(e)=pal(i);%pal matrix stores the pause time of individual microtubules
    end
end

X=0.1:0.1:5.7;h=0;s=0;
set(1:2,1:2)=0;t=0;
for l=1:numel(X)
    w=0;
    temp=round(X(l),1);
    w=sum(round(cal,1)==temp);
    t=0;
    for k=1:numel(cap)
        if (round(cal(k),1)==round(temp,1))
            t=t+cap(k);
        end
    end
    set(1,l)=temp;
    if(t==0 && w==0)
        set(2,l)=0;
    else
        set(2,l)=t/w;
    end
end
end