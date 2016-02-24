function [B] = compute_bearing(pos_x,pos_y)
%% normalize spatial coordinates
[x] = normalize_data(pos_x);
[y] = normalize_data(pos_y);
%%
B = zeros(1,length(x-1));
for it = 2:length(x)

    t= atan2(y(it)-y(it-1),x(it)-x(it-1));
    if t ~= 0 || ~isnan(t)
        B(it) = t;
    else
        B(it-1) = B(it-e);
    end;
end;
B = B+pi;

return;
%code by F.Roux, Sept 2015
