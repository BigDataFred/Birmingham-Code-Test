function [phiXY] = phi2xy(pos_x,pos_y,pos_t,phi,lfp_t)
%%
idx = cell(length(pos_t),1);
parfor it = 1:length(pos_t)
    %if it <length(pos_t)
    %    idx{it} = find(lfp_t >= pos_t(it) & lfp_t < pos_t(it+1));
    %else
        %idx{it} = find(lfp_t >= pos_t(it));
    %end;

    idx{it} = find(abs(lfp_t-pos_t(it)) == min(abs(lfp_t - pos_t(it))));
end;
if length(idx) ~= length(pos_t)
    error('number of elements does not match');
end;
%%
phi = phi+pi;
deg = rad2deg(phi);
deg2 = zeros(length(idx),1);
for it= 1:length(idx)
    deg2(it) = mean(deg(idx{it}));
end;

phiXY.phi = deg2;
phiXY.pos_t = pos_t;
phiXY.pos_x = pos_x;
phiXY.pos_y = pos_y;

return;
%code by F.Roux, Sept 2015