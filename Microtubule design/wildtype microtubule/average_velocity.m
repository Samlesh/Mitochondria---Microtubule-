function[vavg]=average_velocity(lls,vstor)
lenref=0.1:0.1:5.8;
vavg(58)=0;
for k=1:numel(lenref)
    temp=round(lenref(k),1);
    sum_vel=0;
    h=0;
    for i=1:numel(lls)
        if (round(lls(i),1)==temp)
            h=h+1;
            sum_vel=sum_vel+vstor(i);
        end        
    end
if sum_vel==0 && h==0
    vavg(k)=0;
else
vavg(k)=sum_vel/h;    
end
end