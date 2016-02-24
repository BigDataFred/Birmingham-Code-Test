function [stP,idx] = spike_triggered_P(sp_idx,lfp,Fs,tw)

N = fix(Fs*tw/2);
stP = zeros(length(sp_idx),2*(N)+1);
idx = zeros(length(sp_idx),1);

k = 0;
for it = 1:length(sp_idx)
    
    if sign(sp_idx(it)-N) ==1
        k = k+1;
        stP(k,:) = lfp(sp_idx(it)-N:sp_idx(it)+N);
        idx(k) = it;
    end;
    
end;
stP(k+1:end,:) = [];
idx(k+1:end) = [];

return;

%code by F. Roux, Sept 2015