%%
f1 = 11;
f2 = 11;

Fs = 1000;
t = 0:1/Fs:1;
sig = sin(2*pi*f1.*t);
sig = sig+0.15*max(sig)*randn(1,length(t));
%%
sig2 = [zeros(1,length(sig)*2) sig zeros(1,length(sig)*2)];
sig2 = conv(sig2,gausswin(floor(1/f2*250)),'same');
%%
sig = sig2(length(sig)*2+1:end-length(sig)*2);

sig = (sig-min(sig))./(max(sig)-min(sig));
%% get local maxima
[~,imax,~,~] = extrema(sig);

imax(sig(imax) < .5) = [];

x1 = find(sig > max(sig)*0.8);
x1 = intersect(x1,imax);
lm1 = x1;
%% get local minima
[~,~,~,imin] = extrema(sig);

imin(sig(imin) > .5) = [];

x2 = find(sig < max(sig)*0.2);
x2 = intersect(x2,imin);
lm2 = x2;
%% get zero crossings
sig3 = sig-mean(sig);
t1 = sig3(1:length(sig3)-1);
t2 = sig3(2:length(sig3));
tt = t1.*t2;
x3 = find(sign(tt)==-1);

k1 = 0;
k2 = 0;
zcn = zeros(length(x3),1);
zcp = zeros(length(x3),1);
for it = 1:length(x3)
    
    d = diff(sig([x3(it)-1 x3(it)]));
    
    if sign(d)==1
        k1 = k1+1;
        zcp(k1) = x3(it);
    else
        k2 = k2+1;
        zcn(k2) = x3(it);
    end;
    
end;
zcp(k1+1:end) = [];
zcn(k2+1:end) = [];

zcp(diff(zcp)==1) = [];
zcn(diff(zcn)==1) = [];
%%
figure;
subplot(311);
hold on;
plot(t,sig);
plot(t(lm1),sig(lm1),'ro');
plot(t(lm2),sig(lm2),'bo');
plot(t(zcn),sig(zcn),'mo');
plot(t(zcp),sig(zcp),'go');
ylabel('Amplitude [a.u.]');
xlabel('Time [s]');
axis tight;
xlim([0 .5]);

%%
phi = [];
phi(1) = 0;
phi(2) = pi/2;
phi(3) = pi;
phi(4) = 3/2*pi;

k = 1;
for it  =2:length(lm1)
    k = k+4;
    phi(k) = phi(k-4)+2*pi;
end;
k = 2;
for it  =2:length(zcn)
    k = k+4;
    phi(k) = phi(k-4)+2*pi;
end;
k = 3;
for it  =2:length(lm2)
    k = k+4;
    phi(k) = phi(k-4)+2*pi;
end;
k = 4;
for it  =2:length(zcp)
    k = k+4;
    phi(k) = phi(k-4)+2*pi;
end;
%%
x = sort([lm1 zcn' lm2 zcp']);
v = phi;

xq = 1:length(t);
vq = interp1(x,v,xq,'linear');
%%
deg = vq.*180/pi;
deg = mod(deg,360);
rad = deg.*pi/180;
rad = rad-pi;
%%
subplot(412);
hold on;
plot(t,vq,'mo');
plot(t(x),phi,'b.-');
axis tight;
ylabel('Phase [deg]');
xlabel('Time [s]');

subplot(413);
plot(t,deg,'b.');
axis tight;

subplot(414);
plot(t,rad,'b.');
axis tight;