      %*****建立速度模型********
function [Vp,Vs,Rho]=MakeModel(x,z,flag)
% xn=length(x);_50_180
% zn=length(z);
xn=x;
zn=z;
Vp=zeros(zn,xn);
Vs=zeros(zn,xn);
Rho=zeros(zn,xn);
if flag==1
%********************带小方块高速体介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=150
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000; 
        elseif i<=450+round(1/12*j)
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2000;
        else
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end


for j=1:xn
    for i=1:zn
        if 350<i && i<400 && 550-round(2/3*i)<j&&j<600-round(2/3*i)
            Vp(i,j)=500;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end



for j=1:xn
    for i=1:zn
        if 500<i && i<550 && 10+round(1/3*i)<j&&j<60+round(1/3*i)
            Vp(i,j)=1000;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end

end

if flag==2

% ********************倾斜层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000; 
        elseif i<=150+round(1/12*j)
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2000;
        else
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end

end

if flag==3
 %********************复杂层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000; 
        elseif i<=70+round(1/12*j)
             Vp(i,j)=1700;Vs(i,j)=0;Rho(i,j)=1100;
        elseif i<=180+round(-1/6*j)
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2500;
        else
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end

end
if flag==4
%********************小偏移距复杂层状介质模型**********************
 for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000; 
        elseif i<=70+round(1/12*j)
             Vp(i,j)=1700;Vs(i,j)=0;Rho(i,j)=1100;
        elseif i<=150+round(1/10*j)
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2500;
        else
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
 end
end

if flag==5

% *****************不含反射层*****************************
for j=1:xn
    for i=1:zn
        Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000;
    end
end

end

if flag==6
%********************3层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000;
        elseif i>50&&i<=70
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=1500;
        elseif i>70&&i<=120
             Vp(i,j)=2200;Vs(i,j)=0;Rho(i,j)=1700;
        elseif i>120&&i<=140
             Vp(i,j)=3000;Vs(i,j)=0;Rho(i,j)=2000;
        elseif i>140&&i<=200
             Vp(i,j)=4500;Vs(i,j)=0;Rho(i,j)=3000;
        elseif i>200&&i<=250
            Vp(i,j)=2800;Vs(i,j)=0;Rho(i,j)=2100;
        else 
            Vp(i,j)=4000;Vs(i,j)=0;Rho(i,j)=2800;
        end
    end
end

end

if flag==7
%********************1层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1900;Vs(i,j)=0;Rho(i,j)=1000;
        else
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2500;
        end
    end
end
end
if flag==8
%********************水平层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1900;Vs(i,j)=0;Rho(i,j)=1000; 
        elseif i<=150
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2500;  
        else
             Vp(i,j)=2200;Vs(i,j)=0;Rho(i,j)=4000;
        end
    end
end

end

if flag==9
% ********************5层状介质模型**********************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000;
        elseif i>50&&i<=100
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=2300;
        elseif i>100&&i<=180
             Vp(i,j)=3000;Vs(i,j)=0;Rho(i,j)=4000;
        else
             Vp(i,j)=2000;Vs(i,j)=0;Rho(i,j)=2500;

        end
    end
end
end

if flag==10
% ********************5层中的1层状介质模型*****************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000;
        else
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=2300;

        end
    end
end
end

if flag==11


% ********************5层中的1,2层状介质模型*****************
for j=1:xn
    for i=1:zn
        if i<=50
             Vp(i,j)=1500;Vs(i,j)=0;Rho(i,j)=1000;
        elseif i>50&&i<=100
             Vp(i,j)=2500;Vs(i,j)=0;Rho(i,j)=2300;
        else 
             Vp(i,j)=3000;Vs(i,j)=0;Rho(i,j)=4000;  
        end
    end
end
end