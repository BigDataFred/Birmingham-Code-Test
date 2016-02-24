function [fC_spix,fC_lap_idx] = compute_spikes_per_lap(pos_x,fCx,sel_idxP,pos_t,st_idx)

%%
Fs  = 1/pos_t(2)-pos_t(1);% sampling rate of LED cam
trsh = 30;% arbitrary threshold 
[dx] = abs(pos_x(sel_idxP)-fCx);%measure distance of each spike position from field center

[fC_spix] = find(dx <trsh);%fix is an index vector that indexes the positions associated with a spike within
                       %the field with center coordinate fCx

[dt] = diff(pos_t(st_idx(fC_spix)));%get the difference in the time steps associated with each position
dt = (dt>1/Fs);%threshold for refresh rate

[idx] = find(dt==1);% each cut in the time stamp corresponds to a new lap?

fC_lap_idx = cell(1,length(idx));% index marks the number of putative laps
for it = 1:length(fC_lap_idx)
    
    if it == 1% on the first lap

        a = [1 idx(it)];%get all time stamps before cut
        if sign(diff(a))==-1 %check 
            a(2) = [];
        end;
        fC_lap_idx{it} = 1:max(a);% save number of spikes associated with lap

    % all other laps
    elseif it < length(fC_lap_idx)
        a = [idx(it-1)+1 idx(it)-1];
        if sign(diff(a))==-1
            a(2) = [];
        end;
        fC_lap_idx{it} = min(a):max(a);

    else% make exception for last lap

        a = [idx(it) length(fC_spix)];
        if sign(diff(a))==-1
            a(2) = [];
        end;
        fC_lap_idx{it} = min(a):max(a);

    end;
end;

return;

%code by F.Roux, Sept 2015