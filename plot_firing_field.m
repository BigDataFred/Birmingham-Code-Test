function plot_firing_field(pos_x,pos_y)
%% spatial coordinates

y_bin = unique(sort(round(pos_y*1)/1));
x_bin = unique(sort(round(pos_x*1)/1));

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
    
        r(1) = max(abs(x_bin));
        r(2) = max(abs(y_bin));
else
    d = max(x_bin)*2;
    r = [d/2 d/2];
    th = 0:pi/360:2*pi;
    x_u = r*cos(th)+0;
    y_u = r*sin(th)+0;
end;

%% calculate time spent at each location
hold on;
plot(x_u,y_u,'k','LineWidth',3);
plot(pos_x,pos_y,'r.');

xlabel(gca,'Horizontal position [a.u.]');
ylabel(gca,'Vertical position [a.u.]');

% set(gca,'XLim',[-r(1)-5 r(1)+5]);
% set(gca,'YLim',[-r(2)-5 r(2)+5]);