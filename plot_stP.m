function plot_stP(stP,tw,flag)
t = linspace(-tw/2,tw/2,size(stP,2));

hold on;
imagesc(t,1:size(stP,1),stP);
plot([0 0],[1 size(stP,1)],'Color',[.75 .75 .75],'LineWidth',3);

if flag ==2
    for it = 1:size(stP,1)
        idx = find(stP(it,:) ==max(stP(it,:)));
        d= abs(t(idx)-0);
        idx = idx(find(d==min(d))); 
        plot(t(idx),it,'k.');
    end;
elseif flag ==3
    for it = 1:size(stP,1)
        idx = find(stP(it,:) ==min(stP(it,:)));
        d= abs(t(idx)-0);
        idx = idx(find(d==min(d)));
        plot(t(idx),it,'k.');
    end;
end;
axis xy;axis tight;
ylabel('Spike # sorted by position');
xlabel('Time [s]');

p = get(gca,'position');
axes('position',[p(1) p(2)+.025 p(3) p(4)/2])

plot(t,mean(stP,1),'k','LineWidth',3);
axis tight;axis off;


return;
%code by F.Roux, Sept 2015
