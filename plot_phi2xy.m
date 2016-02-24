function plot_phi2xy(pos_c,phi)

hold on;
plot([pos_c;pos_c],[phi;phi+361],'k.');
plot([min(pos_c) max(pos_c)],[mean(phi) mean(phi)],'r--');
plot([min(pos_c) max(pos_c)],[mean(phi)+360 mean(phi)+360],'r--');
axis tight;
xlabel('Position [cm]');
ylabel('Phase [deg]');
set(gca,'YTick',[0:180:360*2]);

%code by F.Roux, Sept 2015