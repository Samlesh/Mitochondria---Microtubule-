function [ll,tt]=complete(ll,tt,l_de,t_de,pos_de)
o=find(pos_de);%o represents how many depolymerisation events hae occured
for i=1:numel(o)
   c=0;
    for j=1:numel(tt)
        c=0;
        if j==pos_de(o(i))
            for k=1:numel(t_de)
                if k~=numel(t_de)
                    d=t_de(k+1);
                else
                    d=(t_de(1,end));
                end
                if t_de(k)~=0
                    c=c+1;
                    tempt(c)=t_de(k);
                    templ(c)=l_de(k);
                    t_de(k)=0;
                    l_de(k)=0;
                    if d==0
                        break;
                    end
                end
            end
            if  j~=numel(tt)
                w=j+1;
            else
                w=numel(tt);
            end
            z=j-1;
            if pos_de(o(i))~=numel(tt)
                tt=[tt(1:z) tempt(1:end-1) tt(w:end)];
                ll=[ll(1:z) templ(1:end-1) ll(w:end)];
                pos_de(o)=pos_de(o)+numel(tempt)-2;
            else
                tt=[tt(1:z) tempt];
                ll=[ll(1:z) templ];
            end
            clear ('tempt');clear('templ');
            break;
        end
    end
end
