function [CFC] = compute_CFC(data,method)

dim = size(data);

if strcmp(method,'pow-pow')
    
    CFC = zeros(dim(2),dim(2));
    for it = 1:dim(2)
        parfor jt = 1:dim(2)
            data2 = data;
            CFC(it,jt) = corr(data2(:,it),data2(:,jt));
            
        end;
    end;
end;

if strcmp(method,'PAC')
    
    pf = 4:2:20;
    af = 30:4:120;
    sig = data.sig;
    n  = length(sig);
    Fs = data.Fs;
    pbin = 0:pi/9:2*pi;
    pad = [zeros(1,n) sig' zeros(1,n)];

    PAC = zeros(length(af),length(pf),length(pbin));       
    for it = 1:length(pf)
                        
        [phi_data,~,~] = band_pass_filter_LFP(pad,data.Fs,[pf(it)-2 pf(it)+2]);
        phi_data = phi_data(n+1:end-n);
        phi_data = angle(hilbert(phi_data));
        %irad = phi_data+pi;

        [~,irad] = interp_phi(phi_data);        
        irad = irad+pi;

        [irad,s_idx] = sort(irad);
        del_idx = find(isnan(irad));
        irad(del_idx) = [];
        
        for jt = 1:length(af)

            [amp_data,~,~] = band_pass_filter_LFP(pad,Fs,[af(jt)-3 af(jt)+3]);
            amp_data = amp_data(n+1:end-n);
            amp_data = abs(amp_data).^2;
            amp_data = amp_data(s_idx);
            amp_data(del_idx) = [];

            parfor kt = 1:length(pbin)-1
                pbin2 = pbin;
                amp_data2 = amp_data;
                idx = find(irad >= pbin2(kt) & irad < pbin2(kt+1));
                PAC(jt,it,kt) = mean(amp_data2(idx));
            end;                        
        end;
        PAC(:,:,end) = PAC(:,:,1);
        
    end;
    N = size(PAC,3);
    P = PAC./repmat(sum(PAC,3),[1 1 size(PAC,3)]);
    H = -sum(log(P).*P,3);
    CFC.mi = (log(N)-H)./log(N);
    CFC.pac = PAC;
    CFC.af = af;
    CFC.pf = pf;
    CFC.pbin = pbin;
end;

return;

%code by F.Roux, Sept.2015