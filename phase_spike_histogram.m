function [phspH,pbin] = phase_spike_histogram(phi,sp_times,t)
%%
k = 0;
idx = zeros(1,length(t)-1);
for it = 1:length(t)-1
   
    ix = find(sign(sp_times-t(it))~=-1 & sign(sp_times-t(it+1))==-1);
    if ~isempty(ix)
        k = k+1;
        idx(k) = it;
    end;
    
end;
idx(k+1:end) = [];
phi = phi(idx);
%%
phi = phi+pi;
[deg] = rad2deg(phi);
deg = round(deg*1)/1;
%%
pbin = 0:24:360;
phspH = zeros(1,length(pbin));
for kt = 1:length(pbin)-1
    
    [idx] = find(sign(deg-pbin(kt))~=-1 & sign(deg-pbin(kt+1))==-1);
    if isempty(idx)
        error(['empty bin detected @ ', num2str(pbin(kt)),' rad']);
    else
        phspH(kt) = length(idx);
    end;
end;
phspH(end) = phspH(1);


return
% code by F.Roux, Sept 2015