function [fl, lownoise, fh, highnoise] = NLNM(unit)

%returns amplitude spectral densities

PL = [0.1 0.17 0.4 0.8 1.24 2.4 4.3 5 6 10 12 15.6 21.9 31.6 45 70 ...
    101 154 328 600 10000];
AL = [-162.36 -166.7 -170 -166.4 -168.6 -159.98 -141.1 -71.36 -97.26 ...
    -132.18 -205.27 -37.65 -114.37 -160.58 -187.5 -216.47 -185 ...
    -168.34 -217.43 -258.28 -346.88];
BL = [5.64 0 -8.3 28.9 52.48 29.81 0 -99.77 -66.49 -31.57 36.16 ...
    -104.33 -47.1 -16.28 0 15.7 0 -7.61 11.9 26.6 48.75];

PH = [0.1 0.22 0.32 0.8 3.8 4.6 6.3 7.9 15.4 20 354.8 10000];
AH = [-108.73 -150.34 -122.31 -116.85 -108.48 -74.66 0.66 -93.37 73.54 ...
    -151.52 -206.66 -206.66];
BH = [-17.23 -80.5 -23.87 32.51 18.08 -32.95 -127.18 -22.42 -162.98 ...
    10.01 31.63 31.63];

fl = 1./PL;
fh = 1./PH;
lownoise = 10.^((AL + BL.*log10(PL))/20);
highnoise = 10.^((AH + BH.*log10(PH))/20);


if unit==1 %displacement
    lownoise = lownoise.*(PL/(2*pi)).^2;
    highnoise = highnoise.*(PH/(2*pi)).^2;
end

if unit==2 %speed
    lownoise = lownoise.*(PL/(2*pi));
    highnoise = highnoise.*(PH/(2*pi));
end

if unit==3 %acceleration
    lownoise = lownoise;
    highnoise = highnoise;
end