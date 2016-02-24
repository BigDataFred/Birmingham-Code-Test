function [st_idx,spc,Fs] = spike_time_index(pos_time,spike_times)
%%
Fs = 1/mean(diff(pos_time));

sel_idx = zeros(length(pos_time),1);
for it = 1:length(pos_time)-1
    
    idx = find(spike_times>= pos_time(it) & spike_times < pos_time(it+1));%+1/Fs
    if ~isempty(idx)
       sel_idx(it) = length(idx);
    end;
end;

x = zeros(length(sel_idx),1);for it = 1:length(sel_idx);x(it) = sel_idx(it);end;
if sum(x) ~= length(spike_times); 
    error('vectors must be of same length');
end;

st_idx = find(x>0);
spc = x(st_idx);