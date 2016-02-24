function [ideg,irad] = interp_phi(sig)

%% normalize
sig = normalize_data(sig);
%% find extremas
[~,imax,~,~] = extrema(sig);

imax(sig(imax) < .5) = [];
%% get local maxima
lm1 = sort(imax);
%% create vector with phase values corresponding to local extremas and zero crossings
n = length(lm1)*2;%
phi = zeros(1,n);
phi(1) = 0;
phi(2) = pi;

c = 2*pi;
k = 1;
for it  =1:length(lm1)-1
    k = k+2;
    phi(k) = phi(k-2)+c;
end;
k = 2;
for it  =1:length(lm1)-1
    k = k+2;
    phi(k) = phi(k-2)+c;
end;
%% interpolate indexes for local minima
v = sort(lm1);
x = phi(1:2:end);

xq = phi(2:2:end);

vq = interp1(x,v,xq,'linear');
vq = fix(vq);

v = phi;
x = sort([lm1 vq]);

v(end) = [];
x(end) = [];
%% do linear interpolation of phase
%x = sort([lm1 zcn' lm2 zcp']);

xq = 1:length(sig);
vq = interp1(x,v,xq,'linear');
%% convert back to deg and rad
deg = vq.*180/pi;
[ideg] = mod(deg,360);

rad = ideg.*pi/180;
[irad] = rad-pi;
%% check consistency

if round(max(ideg)-min(ideg)*1)/1 == 360 ~=1 ||round(2*pi) == round(max(irad)-min(irad)) ~=1
    error('phase values are out of range');
end;
% code by F. Roux, Sept 2015