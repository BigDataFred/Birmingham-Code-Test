function [spAC] = rotate_autocorr(FRM,theta)

dim = [size(FRM,1),size(FRM,2)];
spAC = zeros(length(theta),max(dim),max(dim));
for jt = 1:length(theta)

    [R] = imrotate(FRM,theta(jt),'nearest','crop');%probably this can be optimized
    
   
    spAC(jt,:,:) = corr(FRM,R);

end;
spAC(isnan(spAC)) = 0;

return

%code by F.Roux, Sept. 2015