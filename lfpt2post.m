function [idx] = lfpt2post(lfpt,post)

Fs1  = 1/(lfpt(2)-lfpt(1));
Fs2 = 1/(post(2)-post(1));

if sign(Fs2-Fs1)==1
    error('sampling rate mismatch');
end;
%returns indexes of lfp time bins corresponding to position time bins
idx = zeros(1,length(lfpt));
for it = 1:length(lfpt)
    
    d = abs(post-lfpt(it));
    
    idx(it) = find(d == min(d));% get the closest value
 
end;

chck = any(sign(diff(idx)))==-1;
if chck
    error('index values must be monotonically increasing');
end;

return
%code by F.Roux, Sept 2015