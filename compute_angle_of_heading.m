function [dfr] = compute_bearing(pos_x,pos_y)
%% normalize spatial coordinates
[x] = normalize_data(pos_x);
[y] = normalize_data(pos_y);

R = zeros(1,length(x)-1);
TH = zeros(1,length(x)-1);
for it = 2:length(x)
    
   dx = (diff([x(it-1) x(it)]));%change along x-axis
   dy = (diff([y(it-1) y(it)]));%change along y-axis
   t = atan(dy/dx);% angle of hypothenuse to base
   r = sqrt(dx.^2+dy.^2);
   if sign(dx)==-1
       if sign(dy)==1
           t = t+pi;
       else
           t = t-pi;
       end;
   end;
   R(it-1) = r;
   TH(it-1) = t;
    
end;
TH = TH(sp_idx);
R = R(sp_idx);
TH = rad2deg(TH);
% %% 
% x = pos_x(rmv_idx);
% y = pos_y(rmv_idx);
% 
% THr = zeros(1,length(x)-1);
% Rr = zeros(1,length(x)-1);
% for it = 2:length(x)
%     dx = (diff([x(it-1) x(it)]));
%     dy = (diff([y(it-1) y(it)]));
%     t = atan(dy/dx)*180/pi;
%     Rr(it-1) = sqrt(dy.^2+dx.^2);
%     if sign(dy)==1
%         THr(it-1) = t+180;
%     else
%         THr(it-1) = t-180;
%     end;
% end;
% %%
% [bf] = linspace(-90,90,30);
% nf = zeros(length(bf),1);
% for it = 1:length(bf)-1
%     idx = find(THf >= bf(it) & THf < bf(it+1));
%     nf(it) = length(idx);
% end;
% nf = nf./length(nf);
% %%
% [br] = [linspace(-180,-90,15) linspace(90,180,15)];
% nr = zeros(length(br),1);
% for it = 1:length(br)-1
%     idx = find(THr >= br(it) & THr < br(it+1));
%     nr(it) = length(idx);
% end;
% nr = nr./length(nr);
% %%
% dfr(1).theta = THf.*pi/180;
% dfr(1).rho = Rf;
% dfr(1).b = bf.*pi/180;
% dfr(1).fr = nf';
% 
% dfr(2).theta = THr.*pi/180;
% dfr(2).rho = Rr;
% dfr(2).b = br.*pi/180;
% dfr(2).fr = nr';
%%

return;
%code by F.Roux, Sept 2015
