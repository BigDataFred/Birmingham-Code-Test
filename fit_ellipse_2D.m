function [ex,ey] = fit_ellipse_2D(Z)


[i1,i2] = find(Z ~=0);

[oy,ox] = find(Z == max(max(Z)));

[y1] = min(i1);
[y2] = max(i1);

[x1] = min(i2);
[x2] = max(i2);


b = y2-y1;%vertical radius
a = x2-x1;%horizontal radius

t = -pi:0.01:pi;
ex = ox+a/2*cos(t);
ey = oy+b/2*sin(t);

return;
% code by F.Roux, Sept. 2015