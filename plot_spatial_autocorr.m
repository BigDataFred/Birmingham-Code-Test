function plot_spatial_autocorr(M,spAC)

C = xcorr2(M);

Z = normalize_data(C);
m = ones(size(Z,1),size(Z,2))*mean(Z(:));
sd = ones(size(Z,1),size(Z,2))*std(Z(:));
Z = (Z-m)./sd;

Z = Z.*(Z>1.25);

[ex,ey] = fit_ellipse_2D(Z);

hold on;
imagesc(Z);
plot(ex,ey,'w');
axis xy;
axis off;
return

%code by F. Roux, Sept. 2015