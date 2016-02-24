function [S,f,Serr] = compute_power_spectrum_LFP(data,Fs)
%%
params = [];
params.Fs = Fs;
params.fpass = [0 params.Fs/2];
params.pad = 0;
params.tapers = [25 49];
params.err = [1 0.01];

[S,f,Serr] = mtspectrumc(data,params);