function [S,t,f] = compute_TFR(data,Fs)
%%
params = [];
params.Fs = Fs;
params.fpass = [0 params.Fs/2];
params.pad = 2;
params.tapers = [3 5];
params.err = [0 0.01];

movingwin = [1.5 .5];

[S,t,f] = mtspecgramc(data,movingwin,params);

%code by F.Roux, Sept 2015