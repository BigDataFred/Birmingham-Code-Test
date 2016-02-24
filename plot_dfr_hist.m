function plot_dfr_hist(dfr,binsize,C)

if isempty(binsize)
    binsize = 90;
end;

r = 0:.1:1;g = 0:.1:1;b=0:.1:1;

T = cell(length(dfr),1);
R = cell(length(dfr),1);
h = zeros(length(dfr),1);
hold on;
for it = 1:length(dfr)
    
    T{it} = dfr(it).T;
    R{it} = dfr(it).R;   
    
end;
    
yL = [max(unique(sort([R{:}]))) max(unique(sort([R{:}])))];
yL = [-max(yL) max(yL)];
plot(yL,[0 0],'k');
plot([0 0],yL,'k');

for it = 1:length(dfr)  
    
    h(it) = polar(T{it},R{it});
    
    if isempty(C)
        ix =randperm(length(g));
        set(h(it),'Color',[r(ix(1)) g(ix(6)) b(ix(11))]);
    else
        set(h(it),'Color',C);
    end;
    
end;

axis tight;axis off;


return;
%code by F.Roux, Sept 2015