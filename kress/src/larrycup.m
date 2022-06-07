function [src_info] = larrycup(a,b,n,m)
    src_info = zeros(8,n);
    nhalf = ceil(m/2);
    s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
    r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
    c = a; %*(1-b/pi);  % is theta rounding scale
    sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
    th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
    rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
    z = z*1.2;  % vert stretch! makes ellipse cavity
    Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve
    zhat = fft(Z(:))/m;
    t1 = (0:(n-1))/n;
    h = 1.0/n;
    xy = fourierZ(zhat,t1);
    dxydt = fourierZp(zhat,t1);
    d2xydt = fourierZpp(zhat,t1);
    d3xydt = fourierZppp(zhat,t1);
    src_info(1,:) = real(xy);
    src_info(2,:) = imag(xy);
    src_info(3,:) = real(dxydt)/2/pi;
    src_info(4,:) = imag(dxydt)/2/pi;
    src_info(5,:) = real(d2xydt)/4/pi/pi;
    src_info(6,:) = imag(d2xydt)/4/pi/pi;
    src_info(7,:) = real(d3xydt)/8/pi/pi/pi;
    src_info(8,:) = imag(d3xydt)/8/pi/pi/pi;

end
