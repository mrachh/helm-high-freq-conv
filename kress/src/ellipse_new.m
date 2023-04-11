function [srcinfo,h] = ellipse_new(a,b,n)
    srcinfo = zeros(6,n);
    ts = (0:n-1)/n*2*pi;
    srcinfo(1,:) = a*cos(ts);
    srcinfo(2,:) = b*sin(ts);
    srcinfo(5,:) = sqrt((b*cos(ts)).^2 + (a*sin(ts)).^2);
    srcinfo(3,:) = b*cos(ts)./srcinfo(5,:);
    srcinfo(4,:) = a*sin(ts)./srcinfo(5,:);
    d2xdt2 = -a*cos(ts);
    d2ydt2 = -b*sin(ts);
    dxdt = -a*sin(ts);
    dydt = b*cos(ts);
    srcinfo(6,:) = (dxdt.*d2ydt2 - dydt.*d2xdt2)./srcinfo(5,:).^3; 
    h = 2*pi/n;
end
