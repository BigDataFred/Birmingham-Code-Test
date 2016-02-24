function [Y] = smooth2D(X,bin_size)
gkw = gaussian2d(bin_size,2);
Y = conv2(X,gkw);
d1 = size(Y,1)-size(X,1);
d2 = size(Y,2)-size(X,2);

Y = Y(round(d1/2)+1:end-floor(d1/2),:);
Y = Y(:,round(d2/2)+1:end-floor(d2/2));

%code by F.Roux, Sept 2015