%%
restoredefaultpath;
addpath(genpath('~/froux/chronux_2_10/'));
addpath('~/froux/mFiles/shadedErrorBar');
addpath(genpath('~/froux/Rodent_data/mcode/'));
%% load the data
p2f = '/bcbl/home/home_a-f/froux/Rodent_data/Contest_Simon/';
eeg_dat = load([p2f,'rawlfp_spk.mat']);
pos_dat = load([p2f,'positiondata.mat']);
%% get position and time stamps
pos_timeStamps = pos_dat.position.time;
pos_x = pos_dat.position.x;
pos_y = pos_dat.position.y;

sp_dat.ts = eeg_dat.data.time(eeg_dat.data.spk)';
%% convert postition time to lfp time - returns one index of LFP time per spatial bin
[tidx1] = post2lfpt(pos_timeStamps,eeg_dat.data.time);% get position times corresponding to lfp samples
%% convert lfp time to position time - returns indexes of LFP for each position bin
[tidx2] = lfpt2post(eeg_dat.data.time,pos_timeStamps);% get lfp samples corresponding to position times
%% compute bearing
[B] = compute_bearing(pos_x,pos_y);
%% separate positions in quandrants
[st_idx,spc,Fs] = spike_time_index(pos_timeStamps,sp_dat.ts);
[nposx] = normalize_data(pos_x);
[nposy] = normalize_data(pos_y);

nposx = nposx(st_idx);
nposy = nposy(st_idx);

idx = {};
idx{1} = intersect(find(nposy>.5),find(nposx>.5));% upper right quadrant
idx{2} = intersect(find(nposy>.5),find(nposx<.5));% upper left quadrant
idx{3} = intersect(find(nposy<.5),find(nposx>.5));% lower right quadrant
idx{4} = intersect(find(nposy<.5),find(nposx<.5));% lower left quadrant
%% compute the directional FR in each quadrant
dfr = struct;
for it = 1:length(idx)
    [dfr(it).T,dfr(it).R] = compute_directional_FR(B(st_idx(idx{it})),0:pi/12:2*pi);
end;
%% compute the 2D firing map
[FRM,x_bin,y_bin] = compute_rateMap(pos_x,pos_y,pos_timeStamps,sp_dat.ts,3);
%% compute the gridness scores
m = ones(size(FRM)).*mean(mean(FRM));
sd = ones(size(FRM)).*mean(mean(FRM));

Z = (FRM-m)./sd;
Z = Z.*(Z>2);

th = [0:6:180];
[spAC] = rotate_autocorr(FRM,th);

%attempt to calculate gridness but probably needs correction for
%ellipticity
a = [min(min(spAC(find(th==60),:,:))) min(min(spAC(find(th==120),:,:)))];
b = [max(max(spAC(find(th==30),:,:))) max(max(spAC(find(th==90),:,:))) max(max(spAC(find(th==150),:,:)))];
G = min(a)-max(b); 
%% band-pass filter LFP
Fs = eeg_dat.data.fs; 
sig = eeg_dat.data.lfp';
n  = length(sig);
pad = [zeros(1,n) sig zeros(1,n)];

[filt_data,Sf,ff] = band_pass_filter_LFP(pad,Fs,[5 11]);
filt_data = filt_data(n+1:end-n);

[phi] = angle(hilbert(filt_data));

phi_est = sin(phi)+cos(phi);

[phi] = angle(hilbert(phi_est));

[ideg,irad] = interp_phi(phi);
%% phase spike histo
[phspH,pbins] = phase_spike_histogram(phi,sp_dat.ts,eeg_dat.data.time);
%% spike train autocorrelation
Fs = eeg_dat.data.fs;
N = eeg_dat.data.fs*.3;

sp_tr = zeros(1,length(eeg_dat.data.time));
sp_tr(eeg_dat.data.spk) = 1;
sp_tr = conv(sp_tr,gausswin(floor(0.05*Fs)),'same');

[c,lags] = xcorr(sp_tr,N,'coeff');

% do a shifting of spike train to destroy temporal order
rdc = zeros(1000,size(c,2));
for it = 1:1000;    
    
    [shufx]= shift_signal(sp_tr,[3:.25:100],Fs,1);
    [rdc(it,:)] = xcorr(sp_tr,shufx,N,'coeff');
end;

c(c > mean(c)+0.25*std(c)) = 0;

lags = lags./Fs;
%% spike triggered potential
sel_idx = eeg_dat.data.spk;
[~,s_idx] = sort(pos_x(tidx2(sel_idx)'));
sel_idx = sel_idx(s_idx);

[stP1,stP_idx1] = spike_triggered_P(sel_idx,eeg_dat.data.lfp,eeg_dat.data.fs,.3);
parfor it = 1:size(stP1,1)
    stP1(it,:) = normalize_data(stP1(it,:));
end;

deg = rad2deg(phi+pi);
[stP2,stP_idx2] = spike_triggered_P(sel_idx,deg,eeg_dat.data.fs,.3);
%% spatio temporal phase-firing map
sel_idx = st_idx;
[phiXY] = phi2xy(pos_x(sel_idx),pos_y(sel_idx),pos_timeStamps(sel_idx),phi,eeg_dat.data.time);
[stFRM,b1,b2] = spatio_temporal_fRM(phiXY.pos_x,phiXY.phi);
%% laps through field in lower quadrant
sel_idxSP = [idx{4}];% indexes for spikes that occured at coords in ll quad
sel_idxSP = sort(sel_idxSP);% sort according to time

sel_idxP= st_idx(sel_idxSP);%indexes of coordinates in llquad for spikes

if max(nposx(idx{4})) >.5
    error('coordinate out of range');
end;

if length(sel_idxSP) ~= length(sel_idxP)
    error('wrong index assignment');
end;

[fm1,fm2] = find(FRM == max(max(FRM)));

fCx = x_bin(fm2);% x coordinate of field with max firing
fCy = y_bin(fm1);

[fC_spix,fC_lap_idx] = compute_spikes_per_lap(pos_x,fCx,sel_idxP,pos_timeStamps,st_idx);

sel_idx = eeg_dat.data.spk(sel_idxSP(fC_spix));
[stP3,stP_idx3] = spike_triggered_P(sel_idx,eeg_dat.data.lfp,eeg_dat.data.fs,.3);
parfor it = 1:size(stP3,1)
    stP3(it,:) = normalize_data(stP3(it,:));
end;

[phiXY2] = phi2xy(pos_x(sel_idxP(fC_spix)),pos_y(sel_idxP(fC_spix)),pos_timeStamps(sel_idxP(fC_spix)),phi,eeg_dat.data.time);
[stFRM2,b3,b4] = spatio_temporal_fRM(phiXY2.pos_x,phiXY2.phi);
%% measure theta phase precession
M = [phiXY.phi phiXY.pos_x];% spike phase x spike position

[Rsq,Rx] = phase_precession(M);%compute the R²
%% compute spike rate
[mfit] = locfit(sp_dat.ts,'family','rate','h',.05,'nn',0.3,'kern','gauss');
%% spectra and spike field coherence
[S_lfp,f_lfp,Serr_lfp] = compute_power_spectrum_LFP(eeg_dat.data.lfp,eeg_dat.data.fs);%spectrum LFP
[S_sp,f_sp,Serr_sp] = compute_power_spectrum_spikes(sp_dat.ts,eeg_dat.data.fs);%spectrum spikes
[C,~,Sxy,Sx,Sy,t,f,zerosp] = compute_spike_field_locking(eeg_dat.data.lfp,sp_dat.ts,[1 .5],eeg_dat.data.fs);%spike field coherence
%% Cross frequency
[TFR,t2,f2] = compute_TFR(eeg_dat.data.lfp,eeg_dat.data.fs);
[CFC] = compute_CFC(TFR,'pow-pow');
%%
indat.sig = eeg_dat.data.lfp;
indat.Fs = eeg_dat.data.fs;
[PAC] = compute_CFC(indat,'PAC');
clear indat;
%%
pad = [zeros(1,n) eeg_dat.data.lfp' zeros(1,n)];
[amp_data,~,~] = band_pass_filter_LFP(pad,eeg_dat.data.fs,[100 110]);
amp_data = amp_data(n+1:end-n);

[~,imax,~,~] = extrema(amp_data);
imax = sort(imax);
[stP3,stP_idx3] = spike_triggered_P(imax,eeg_dat.data.lfp,eeg_dat.data.fs,.15);
%% visuzlize spiking field
figure;
subplot(411);
hold on;
plot(pos_x,pos_y,'Color',[.25 .25 .25]);
plot_firing_field(pos_x(st_idx),pos_y(st_idx));
axis tight;
subplot(412);
plot(x_bin,mean(FRM,1),'k','LineWidth',3);
axis tight;
xlabel('Position [a.u.]');
ylabel('Firing rate [Hz]');
subplot(413);
plot_firing_map(x_bin,y_bin,FRM);
axis tight;
subplot(414);
plot_spatial_autocorr(Z);
axis off;
%%
figure;
subplot(4,4,1:3);
hold on;
plot(nposx,nposy,'Color',[.25 .25 .25]);
plot(nposx(idx{1}),nposy(idx{1}),'bo');

subplot(4,4,4);
plot_dfr_hist(dfr(1),90,'b');

subplot(4,4,5:7);
hold on;
plot(nposx,nposy,'Color',[.25 .25 .25]);
plot(nposx(idx{2}),nposy(idx{2}),'ro');

subplot(4,4,8);
plot_dfr_hist(dfr(2),90,'r');

subplot(4,4,9:11);
hold on;
plot(nposx,nposy,'Color',[.25 .25 .25]);
plot(nposx(idx{3}),nposy(idx{3}),'mo');

subplot(4,4,12);
plot_dfr_hist(dfr(3),90,'m');

subplot(4,4,13:15);
hold on;
plot(nposx,nposy,'Color',[.25 .25 .25]);
plot(nposx(idx{4}),nposy(idx{4}),'go');

subplot(4,4,16);
plot_dfr_hist(dfr(4),90,'g');
%% visualize spectral data
figure;
subplot(331);
a = gca;
hold on;
shadedErrorBar(f_lfp(f_lfp>=1),20*log10(S_lfp(f_lfp>=1)+1),squeeze(log10([Serr_lfp(2,(f_lfp>=1),:);Serr_lfp(1,(f_lfp>=1),:)]+1)),'b',1);
axis tight;xlim([0 50]);
title('LFP spectrum');
subplot(332);
a  = [a gca];
hold on;
shadedErrorBar(f_sp(f_sp>=1),20*log10(S_sp(f_sp>=1)+1),squeeze(log10([Serr_sp(2,(f_sp>=1),:);Serr_sp(1,(f_sp>=1),:)]+1)),'r',1);
axis tight;xlim([0 50]);
title('Spike spectrum');
subplot(333);
a  = [a gca];
plot(f(f>=1),mean(C(zerosp==0,(f>=1)),1));
axis tight;xlim([0 50]);
title('Spike-field coherence');

subplot(3,3,4:9);
a  = [a gca];
imagesc(t,f(f>=1),20*log10(Sx(:,f>=1)'+eps));
axis xy;

for it = 1:length(a)-2
    xlabel(a(it),'Frequency [Hz]');
    ylabel(a(it),'Power [dB]');
end;

for it = 3
    xlabel(a(it),'Frequency [Hz]');
    ylabel(a(it),'Coherence [a.u.]');
end;

for it = 4
    xlabel(a(it),'Time [s]');
    ylabel(a(it),'Frequency [Hz]');
end;
%% visualize CFC
figure;
subplot(3,2,1);
plot_CFC(CFC',f2,f2);
xlabel('Frequency for power [Hz]');
ylabel('Frequency for power [Hz]');

subplot(3,2,3);
hold on;
idx1 = find(f2 >= 5 & f2 < 11);
pow1 = mean(TFR(:,idx1),2);

idx2 = find(f2 >= 14 & f2 < 20);
pow2 = mean(TFR(:,idx2),2);

idx3 = find(f2 >= 80 & f2 < 120);
pow3 = mean(TFR(:,idx3),2);

h(1) = plot(pow1,pow2,'ko');
h(2) = plot(pow1,pow3,'o','Color',[.75 .75 .75]);
legend(h,'Beta','Gamma');

axis tight;
xlabel('Theta power');
ylabel('Beta power');

subplot(3,2,2);
pf = 4:2:20;
af = 30:4:120;
pcolor(PAC.pf,PAC.af,PAC.mi);
shading interp;lighting phong;
xlabel('Frequency for phase [Hz]');
ylabel('Frequency for amplitude [Hz]');

subplot(3,2,4),
idx1 = find(PAC.pf >= 6 & PAC.pf < 10);
idx2 = find(PAC.af >= 80 & PAC.af < 120);

Y = squeeze(mean(mean(PAC.pac(idx2,idx1,:),2),1));
Y = normalize_data(Y);
X = 0:pi/9:2*pi;
X = X.*180/pi;
bar(X,Y,1,'k','EdgeColor',[.75 .75 .75]);
axis tight;
xlabel('Theta phase [deg]');
ylabel('Gamma power [a.u.]');

subplot(326);
plot_stP(stP3,.3,[]);
%%
figure;
subplot(321);
hold on;
h(1) = plot(eeg_dat.data.time,eeg_dat.data.lfp,'k');
h(2) = plot(eeg_dat.data.time,filt_data,'Color',[.75 .75 .75]);
xlim([54 56]);
xlabel('Time [s]');
ylabel('Amplitude');
legend(h,'LFP','5-11Hz');

subplot(322);
plot_phase_spike_histogram(pbins,phspH)

subplot(323);
plot_stP(stP1,.3,[]);

subplot(324);
plot_stP(stP2,.3,[]);

subplot(325);
hold on;
bar(lags,c,'k');
plot(lags,ones(1,length(lags))*mean(mean(rdc,1))+mean(2*std(rdc,0,1)),'r--','LineWidth',3);
axis tight;
xlabel('Time lag [s]');
ylabel('Correlation');
%%
figure;
subplot(311);
plot(x_bin,mean(FRM,1),'k','LineWidth',3);
axis tight;
xlabel('Horizontal Position [a.u.]');
ylabel('Average firing rate [Hz]');
subplot(312);
plot_phi2xy(phiXY.pos_x,phiXY.phi);
subplot(313);
plot_spatio_temporal_FF(stFRM,b1,b2);
%%
figure;
subplot(421);
hold on;
plot(pos_x,pos_y,'k-');
plot_firing_field_box(fCx,fCy);
plot(pos_x(sel_idxP),pos_y(sel_idxP),'ro');
xlabel('Horiz. pos. [a.u.]');
ylabel('Vert. pos. [a.u.]');
axis(gca,'tight');

subplot(422);
hold on;
[dx] = abs(pos_x(sel_idxP)-fCx);

plot(sort(dx));
plot([1 length(dx)],[trsh trsh],'r--');
axis tight;
ylabel('Distance from field [a.u.]');
xlabel('Sorted coordinate #');
axis(gca,'tight');

subplot(423);
hold on;
plot(pos_x,pos_y,'k-');

C = sel_idxP(fC_spix);
pix = pos_x(C(1)-200:C(1)+200);
piy = pos_y(C(1)-200:C(1)+200);
plot(pix,piy,'g-','LineWidth',3);
plot(pix(1),piy(1),'wo','MarkerSize',9,'MarkerFaceColor','g');

C = sel_idxP(fC_spix);
pix = pos_x(C(end)-200:C(end)+200);
piy = pos_y(C(end)-200:C(end)+200);
plot(pix,piy,'y-','LineWidth',3);
plot(pix(1),piy(1),'ko','MarkerSize',9,'MarkerFaceColor','y');

plot_firing_field_box(fCx,fCy);
plot(pos_x(C),pos_y(C),'ro');

xlabel('Horiz. pos. [a.u.]');
ylabel('Vert. pos. [a.u.]');
axis(gca,'tight');

subplot(424);
for it = 1:length(fC_lap_idx);
    k1 = it;    
    k2 = it+1;
    
    for jt = 1:length(fC_lap_idx{it});
        x = pos_x(sel_idxP(fC_spix(lap_idx{it}(jt))));
        hold on;plot([x x],[k1 k2],'k-');
    end;
end;
axis(gca,'tight');
xlim(gca,[min(pos_x) max(pos_x)]);
xlabel(gca,'Horizontal position [a.u.]');
ylabel(gca,'Lap number');

subplot(425);
a = gca;
x = deg(eeg_dat.data.spk(sel_idxSP(fC_spix)));
[~,s_idx] = sort(x);
plot_stP(stP3(s_idx,:),.3,[]);
ylabel(a,'Spike sorted by phase');

subplot(426);
a = gca;
plot_phi2xy(phiXY2.pos_x,phiXY2.phi);
set(a,'Xlim',[min(pos_x) max(pos_x)]);

subplot(427);
plot(0:360,Rsq);
axis tight;
xlabel('Phase rotation [deg]');
ylabel('R^2');

subplot(428);
[ixM] = find(Rsq ==max(Rsq));
TPPc = squeeze(Rx(ixM(2),:,:));
X = TPPc(:,1)+360;
Y = TPPc(:,2);

[b,bint,r,rint] = regress(Y,[ones(size(X)) X]);
yp = b(1)+b(2)*X;

hold on;
plot(X,yp,'r');
plot(X,Y,'ko');
axis tight;
xlabel('Rotated theta phase [deg]');
ylabel('Rotated position [a.u.]');
%%
for it= 1:gcf
    set(it,'Color','w');
    set(it,'PaperPositionMode','auto');
end;
%%
savepath = '/bcbl/home/home_a-f/froux/Rodent_data/figures/';
fig_name = 'Birmingham_Contest_Sept2015_fig';
for it = 1:gcf
    print(it,'-r300','-dtiff',[savepath,fig_name,num2str(it),'.tif']);
end;
