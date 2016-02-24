function [T,R] = compute_directional_FR(B,pbins)
    
if isempty(pbins)
    pbins = 0:pi/6:2*pi;
end;

[R,T] = hist(B,pbins);
R([1 end]) = sum(R([1 end]));


return;
%code by F.Roux, Sept 2015
