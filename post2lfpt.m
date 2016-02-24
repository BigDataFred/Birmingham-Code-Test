function [idx] = post2lfpt(post,lfpt)
%returns indexes of position time bins corresponding to lfp time bins
idx = cell(1,length(lfpt));
k=0;
for it = 1:length(lfpt)
    
    ix = [];
    if it < length(lfpt)
        ix = find(post >= lfpt(it) & post  < lfpt(it+1));
    else
        ix = find(post >= lfpt(it));
    end;
    
    if ~isempty(ix)
        k=k+1;
        idx{k} = it;
    end;
    
end;
idx(k+1:end) = [];

%check consistency
if length(idx) ~= length(post) || corr(lfpt([idx{:}])',post) < .99
    error('sorted values mismatch');
end;
return
%code by F.Roux, Sept 2015