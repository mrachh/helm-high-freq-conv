function emat = get_emat(n,m)
    tm = 0:pi/m:2*pi-pi/m;
    fs = -(n-1):1:n;
    emat = exp(1j*tm.'*fs);
end