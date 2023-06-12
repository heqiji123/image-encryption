%% 求解四翼混沌
% x' = ax+byz
% y' = cy+dxz
% z' = exy+kz+mxw
% w' = ny
function fy=Chaos_output(x0,y0,z0,w0,iterations,KK)
% 微分方程求解
% opt = odeset('Mass',@mass);
Y=[x0;y0;z0;w0];
[~,y] = ode45(@ode,-50:100/(iterations):50,Y);
fy=y;
    function f = ode(~,y)
        y1=y(1);
        y2=y(2);
        y3=y(3);
        y4=y(4);
        f = [8*y1-y2*y3;-40*y2+y1*y3;2*y1*y2+KK*y3+y1*y4;-2*y2];
    end
end