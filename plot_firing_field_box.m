function plot_firing_field_box(fCx,fCy)

plot([fCx-20 fCx+30],[fCy+20 fCy+20],'b-','LineWidth',3);
plot([fCx-20 fCx+30],[fCy-20 fCy-20],'b-','LineWidth',3);
plot([fCx-20 fCx-20],[fCy-20 fCy+20],'b-','LineWidth',3);
plot([fCx+30 fCx+30],[fCy-20 fCy+20],'b-','LineWidth',3);

% code by F. Roux,Sept. 2015