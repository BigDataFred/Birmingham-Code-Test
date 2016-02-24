function [filt_data,S,f] = band_pass_filter_LFP(data,Fs,Bpf)
%% uses the ft_preproc toolxbox
addpath('~/froux/fieldtrip-20150801/preproc/');
%%
N = 6*fix(Fs/Bpf(1));
[filt_data] = ft_preproc_bandpassfilter(data,Fs,Bpf,N,'fir','twopass','no',[],[],[],[],'yes');

y = fft(data,512)./Fs;
y = fftshift(y);
y = y(round(length(y)/2)+1:end);
S.y = real(y).^2+imag(y).^2;

y2 = fft(filt_data,512)./Fs;
y2 = fftshift(y2);
y2 = y(round(length(y2)/2)+1:end);
S.y2 = real(y2).^2+imag(y2).^2;

f = Fs/2*linspace(0,1,256);

return;
%code by F.Roux, Sept 2015
