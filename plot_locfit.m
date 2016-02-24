function plot_locfit(dsp,fit)
% Take in spike times (in dsp) and the fit (in fit) and generate a plot
% with error bars
% 
N=length(dsp);
xfit = lfmarg(fit);
% placing 'band','g' before varargin{:} ensures that
% user-provided 'band' has precedence.
ypp = predict(fit,xfit,'band','g');
yfit = ypp{1};
se = ypp{2};
bands = ypp{3};

data = fit.data;
xdata = data.x;
p = size(xdata,2);
cv = 1.96;
fali = fit.fit_points.family_link;
cl = invlink(bands(:,1),fali);
cu = invlink(bands(:,2),fali);
yfit = invlink(yfit,fali);

plot(xfit{1},yfit/N,'Linewidth',1.5);
hold on;
plot(xfit{1},cu/N,':','Linewidth',1.5);
plot(xfit{1},cl/N,':','Linewidth',1.5);
yy=get(gca,'ylim');
yy=linspace(yy(1),yy(2),N);

for n=1:N;
    plot(dsp(n).times,ones(length(dsp(n).times),1),'rd', 'Markersize',1.5);
end;