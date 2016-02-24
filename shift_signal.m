function [shif_sig] = shift_signal(sig,tw,Fs,N)

shif_sig = zeros(N,length(sig));
for it = N
        
    k = randperm(length(tw));
    k = k(end);
    
    iv = randperm(length(sig));
    iv = iv(1);
    shfix = iv:iv+tw(k)*Fs;
    
    while shfix(end) > length(sig)
        iv = randperm(length(sig));
        iv = iv(1);
        shfix = iv:iv+tw(k)*Fs;
    end;
    
    c = randperm(2);
    c = c(1);
    if c==1
        s = [shfix(end)+1:length(sig) shfix 1:shfix(1)-1];
        shif_sig(it,:) = sig(s);
    else
        s = [shfix shfix(end)+1:length(sig)  1:shfix(1)-1];
        shif_sig(it,:) = sig(s);
    end;
end;

return;
%code by F.Roux,Sept 2015
