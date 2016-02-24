function [stFRM,pos_bin,p_bin] = spatio_temporal_fRM(pos,phi)

pos_bin = linspace(floor(min(pos)),round(max(pos)),100);
p_bin = 0:8:360;

stFRM = zeros(length(p_bin),length(pos_bin));
k=0;
for it = 1:length(pos_bin)-1
    
    if it < length(pos_bin)
        idx1 = find(pos >= pos_bin(it) & pos < pos_bin(it+1));
    else
        idx1 = find(pos >= pos_bin(it));
    end;
    
    for jt = 1:length(p_bin)
        
        if jt < length(p_bin)
            idx2 = find(phi >= p_bin(jt) & phi < p_bin(jt+1));        
        else
            idx2 = find(phi >= p_bin(jt));
        end;
        
        idx3 = intersect(idx1,idx2); 
        if ~isempty(idx3)
            k = k+1;
            n(k) = length(idx3);
            stFRM(jt,it) = n(k);%/(n*0.02);
        else
            stFRM(jt,it) = 0;

        end;
        
    end;
end;

[stFRM] = smooth2D(stFRM,4);
return;
%code by F.Roux, Sept 2015