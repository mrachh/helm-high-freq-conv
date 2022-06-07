function qmat = get_qmat(n)
    fs = -(n-1):1:n;
    ts = 0:pi/n:2*pi-pi/n;
    qmat = 1/(2*n)*exp(-1j*fs.'*ts);
end