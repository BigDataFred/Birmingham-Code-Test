function plot_phase_spike_histogram(pbins,phspH)
%%
bar([pbins pbins+361],[phspH phspH],26,'k');axis tight;box off;
xlabel('Phase [deg]');ylabel('Spike count');
set(gca,'XTick',0:180:360*2);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);