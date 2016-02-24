function plot_firing_map(x_bin,y_bin,FRM)
%% spatial coordinates
if length(y_bin) < length(x_bin)
    x_u = [min(x_bin) max(x_bin);
        min(x_bin) max(x_bin);
        min(x_bin) min(x_bin)
        max(x_bin) max(x_bin)];
        
     y_u = [min(y_bin) min(y_bin);
        max(y_bin) max(y_bin);
        min(y_bin) max(y_bin);
        min(y_bin) max(y_bin);
        ];
    
        r(1) = max(x_bin);
        r(2) = max(y_bin);
else
    d = max(x_bin)*2;
    r = [d/2 d/2];
    th = 0:pi/360:2*pi;
    x_u = r*cos(th)+0;
    y_u = r*sin(th)+0;
end;
%% calculate time spent at each location
imagesc(x_bin,y_bin,FRM);
axis(gca,'tight');
axis(gca,'xy');

hold on;
plot(x_u,y_u,'k','LineWidth',3);


xlabel(gca,'Horizontal position [a.u.]');
ylabel(gca,'Vertical position [a.u.]');
%caxis(gca,[0 200]);
% set(gca,'XLim',[-r(1)-5 r(1)+5]);
% set(gca,'YLim',[-r(2)-5 r(2)+5]);