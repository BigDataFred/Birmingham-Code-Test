function [ndata] = normalize_data(data)
%%
ndata = (data-min(min(data)))./(max(max(data))-min(min(data)));
