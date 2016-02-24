function plot_spatio_temporal_FF(stFRM,b1,b2)

imagesc(b1,[b2 b2+361],repmat([stFRM./0.02],[2 1]));axis xy;axis tight;
xlabel('Position [cm]');ylabel('Phase [deg]');

%code by F.Roux,Sept 2015