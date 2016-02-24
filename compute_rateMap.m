function [FRM,x_bin,y_bin] = compute_rateMap(pos_x,pos_y,pos_time,spike_times,bin_size)
%%
[st_idx,spc,Fs] = spike_time_index(pos_time,spike_times);
%% create bins for spatial rate map
x_bin = floor(min(pos_x)):bin_size:ceil(max(pos_x));
y_bin = floor(min(pos_y)):bin_size:ceil(max(pos_y));

xc = double(pos_x);
yc = double(pos_y);

xcs = xc(st_idx);
ycs = yc(st_idx);
%%
if length(ycs) ~= length(spc)
    error('vector length must be equal');
end;
%%
SC = zeros(length(y_bin),length(x_bin));
for it = 1:length(x_bin)
    
    if it < length(x_bin)
        sx = find(xcs >= x_bin(it) & xcs < x_bin(it+1));
    else
        sx = find(xcs >= x_bin(it));
    end;
    
    for jt = 1:length(y_bin)
        if jt < length(y_bin)
            idx = intersect(sx,find(ycs >= y_bin(jt) & ycs < y_bin(jt+1)));                
        else
            idx = intersect(sx,find(ycs >= y_bin(jt)));     
        end;
        SC(jt,it) = length(idx);
    end;
end;
%% calculate time spent at each location
tt = zeros(length(y_bin),length(x_bin));
for it = 1:length(x_bin)
    
    if it < length(x_bin)
        sx = find(xc >= x_bin(it) & xc < x_bin(it+1));
    else
        sx = find(xc >= x_bin(it));
    end;
    
    for jt = 1:length(y_bin)
        if jt < length(y_bin)
            idx = intersect(sx,find(yc >= y_bin(jt) & yc < y_bin(jt+1)));                
        else
            idx = intersect(sx,find(yc >= y_bin(jt)));     
        end;
        tt(jt,it) = length(idx)*1/Fs;
    end;
end;
%% normalize spike count by time spent at each location
SC = SC./tt;
SC(isnan(SC)) = 0;
SC(isinf(SC))= 0;
%% smooth map using 2-D gaussian with width = 2SDs
gkw = gaussian2d(bin_size,2);
FRM = conv2(SC,gkw);
d1 = size(FRM,1)-size(SC,1);
d2 = size(FRM,2)-size(SC,2);

FRM = FRM(round(d1/2)+1:end-floor(d1/2),:);
FRM = FRM(:,round(d2/2)+1:end-floor(d2/2));

return;
%code by F.Roux, Sept 2015