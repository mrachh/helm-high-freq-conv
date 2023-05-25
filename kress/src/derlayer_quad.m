function quad = derlayer_quad(n)
% return quadrature in 2n points

j = 0 : 2*n - 1 ;
t = j(2:2:end)*pi/n;

quad          = zeros(2*n,1);
quad(1:2:end) = 0;
quad(2:2:end) = 1./(2 * n * sin(t/2).^2);
quad(1)       = -n/2;
