function [Rsq,Rx] = phase_precession(M)

X(:,1) = M(:,1);%spike phase
X(:,2) = M(:,2);%spike position

theta = 0:360;
Rsq = zeros(length(theta),1);
Rx = zeros(length(theta),size(X,1),2);
for jt = 1:length(theta)    
    
    R=R2d(theta(jt));
    x = X*R;
    dim(jt,:) = size(x);
    [~,~,~,~,stats] = regress(x(:,1),[ones(size(x(:,2))) x(:,2)]);%compute R²
    Rsq(jt) = stats(1);
    Rx(jt,:,:) = x;
end;


return;

%code by F.Roux, Sept. 2015