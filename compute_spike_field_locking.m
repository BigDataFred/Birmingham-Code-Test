function [C,phi,Sxy,Sx,Sy,t,f,zerosp] = compute_spike_field_locking(data1,data2,movingwin,Fs)
%%
if isempty(movingwin)
    movingwin = [.5 .05];
end;

params = [];
params.Fs = Fs;
params.fpass = [0 params.Fs/2];
params.pad = 6;
params.tapers = [3 5];
params.trialave = 1;

[C,phi,Sxy,Sx,Sy,t,f,zerosp]=cohgramcpt(data1,data2,movingwin,params,[]);