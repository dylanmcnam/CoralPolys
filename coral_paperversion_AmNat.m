% the last version of the model we messed with - dimensionless here.
%all the "Plot" and "ParamSpace" scipts are using the same version as this
%one.


clear all
rand('seed',0)
randn('seed',0)

n = 100;      %number of vertical cells
m = n;      %number of horizontal cells

cap = 10;


%PARAMS---------------------------------------

Dh = 0.5;          %spatial resolution
Dt = 0.01;      %time step length
time_steps = 10000;

w = 0; %1 is along axis 0 is normal diffusion

nh = 3; %must be odd #


%[rr,cc]=meshgrid(linspace(.1,6,m),linspace(2,.1,n));
alpha = 4;%DHQ+
beta = 1;%rQ+
C = 0.6;%c
gamma = 0.1;%DQ+




store_Htot = zeros(time_steps,1);
store_Qtot = zeros(time_steps,1);


store_H = zeros(n,m,time_steps/200);
store_Q = zeros(n,m,time_steps/200);

%store_logH = zeros(n,m,time_steps);
%store_logQ = zeros(n,m,time_steps);
%store_advection = zeros(n,m,time_steps);
%store_diffHQ = zeros(n,m,time_steps);
%store_diffQ = zeros(n,m,time_steps);

H = (1 - 0.001*rand(n,m));
Q = (1 - C.*H);% + 0.01*randn(n,m));

H(H<0) = 0;
Q(Q<0) = 0;


prevH=H;

for t =1:time_steps
    
    t
    
    prevH=H;
    
    if(abs(mean(mean(prevH))-mean(mean(H)))>1e-3&t>10)
        disp('cons problem')
        return
    end
    
    %calculuate H and Q neighbors
    Hr = circshift(H,[0 -1]);
    Hrr = circshift(H,[0 -2]);
    Hl = circshift(H,[0 1]);
    Hll = circshift(H,[0 2]);
    Hd = circshift(H,[-1 0]);
    Hdd = circshift(H,[-2 0]);
    Hu = circshift(H,[1 0]);
    Huu = circshift(H,[2 0]);
    
    Hur = circshift(H,[1 -1]);
    Hurur = circshift(H,[2 -2]);
    Hul = circshift(H,[1 1]);
    Hulul = circshift(H,[2 2]);
    Hdr = circshift(H,[-1 -1]);
    Hdrdr = circshift(H,[-2 -2]);
    Hdl = circshift(H,[-1 1]);
    Hdldl = circshift(H,[-2 2]);
    
    Qr = circshift(Q,[0 -1]);
    Ql = circshift(Q,[0 1]);
    Qd = circshift(Q,[-1 0]);
    Qu = circshift(Q,[1 0]);
    
    Qur = circshift(Q,[1 -1]);
    Qul = circshift(Q,[1 1]);
    Qdr = circshift(Q,[-1 -1]);
    Qdl = circshift(Q,[-1 1]);
    
    
    
    HLR   = (Hl  + Hr + Hll + Hrr);
    HUD  = (Hu + Hd + Huu + Hdd);
    HDiagUR  = (Hur  + Hdl + Hurur + Hdldl);
    HDiagLR  = (Hul  + Hdr + Hulul + Hdrdr);
    dum=ones(5,5);
    dum(1,2)=0;
    dum(1,4)=0;
    dum(2,1)=0;
    dum(2,5)=0;
    dum(4,1)=0;
    dum(4,5)=0;
    dum(5,2)=0;
    dum(5,4)=0;
    
    H_tot = convolve2(H,dum,'circular')-H;
    
    v_ud = (HUD)./H_tot;
    v_diagur = (HDiagUR)./H_tot;
    v_lr = (HLR)./H_tot;
    v_diaglr = (HDiagLR)./H_tot;
    
    v_ud(H_tot==0)=0;
    v_diagur(H_tot==0)=0;
    v_lr(H_tot==0)=0;
    v_diaglr(H_tot==0)=0;
    
    v_ud(v_ud>v_diagur&v_ud>v_lr&v_ud>v_diaglr)=1;
    v_lr(v_lr>v_diagur&v_lr>v_ud&v_lr>v_diaglr)=1;
    v_diagur(v_diagur>v_lr&v_diagur>v_ud&v_diagur>v_diaglr)=1;
    v_diaglr(v_diaglr>v_lr&v_diaglr>v_ud&v_diaglr>v_diagur)=1;
    
    v_ud(v_ud<1)=0;
    v_lr(v_lr<1)=0;
    v_diagur(v_diagur<1)=0;
    v_diaglr(v_diaglr<1)=0;
    
    
%     v_ud = ((1-w)*0.25*ones(n,m) + w*v_ud);
%     v_diagur = ((1-w)*0.25*ones(n,m) + w*v_diagur);
%     v_lr = ((1-w)*0.25*ones(n,m) + w*v_lr);
%     v_diaglr = ((1-w)*0.25*ones(n,m) + w*v_diaglr);
%     
    
    
    
    vxr=(v_lr+circshift(v_lr,[0 -1]))./2;
    vxl=(v_lr+circshift(v_lr,[0 1]))./2;
    flux_in=zeros(n,m);
    flux_out=zeros(n,m);
    rt_bigger = find(Hr>H);
    rt_smaller = find(Hr<H);
    lt_bigger = find(Hl>H);
    lt_smaller = find(Hl<H);
    flux_in(rt_bigger)=flux_in(rt_bigger)+(Hr(rt_bigger)-H(rt_bigger)).*vxr(rt_bigger);
    flux_out(rt_smaller)=flux_out(rt_smaller)+(H(rt_smaller)-Hr(rt_smaller)).*vxr(rt_smaller);
    flux_in(lt_bigger)=flux_in(lt_bigger)+(Hl(lt_bigger)-H(lt_bigger)).*vxl(lt_bigger);
    flux_out(lt_smaller)=flux_out(lt_smaller)+(H(lt_smaller)-Hl(lt_smaller)).*vxl(lt_smaller);
    
    x_flux_grad = (flux_in-flux_out)./Dh^2;
    
    vyu=(v_ud+circshift(v_ud,[1 0]))./2;
    vyd=(v_ud+circshift(v_ud,[-1 0]))./2;
    flux_in=zeros(n,m);
    flux_out=zeros(n,m);
    up_bigger = find(Hu>H);
    up_smaller = find(Hu<H);
    dn_bigger = find(Hd>H);
    dn_smaller = find(Hd<H);
    flux_in(up_bigger)=flux_in(up_bigger)+(Hu(up_bigger)-H(up_bigger)).*vyu(up_bigger);
    flux_out(up_smaller)=flux_out(up_smaller)+(H(up_smaller)-Hu(up_smaller)).*vyu(up_smaller);
    flux_in(dn_bigger)=flux_in(dn_bigger)+(Hd(dn_bigger)-H(dn_bigger)).*vyd(dn_bigger);
    flux_out(dn_smaller)=flux_out(dn_smaller)+(H(dn_smaller)-Hd(dn_smaller)).*vyd(dn_smaller);
    
    y_flux_grad = (flux_in-flux_out)./Dh^2;
    
    vur=(v_diagur+circshift(v_diagur,[1 -1]))./2;
    vdl=(v_diagur+circshift(v_diagur,[-1 1]))./2;
    flux_in=zeros(n,m);
    flux_out=zeros(n,m);
    up_dbigger = find(Hur>H);
    up_dsmaller = find(Hur<H);
    dn_dbigger = find(Hdl>H);
    dn_dsmaller = find(Hdl<H);
    flux_in(up_dbigger)=flux_in(up_dbigger)+(Hur(up_dbigger)-H(up_dbigger)).*vur(up_dbigger);
    flux_out(up_dsmaller)=flux_out(up_dsmaller)+(H(up_dsmaller)-Hur(up_dsmaller)).*vur(up_dsmaller);
    flux_in(dn_dbigger)=flux_in(dn_dbigger)+(Hdl(dn_dbigger)-H(dn_dbigger)).*vdl(dn_dbigger);
    flux_out(dn_dsmaller)=flux_out(dn_dsmaller)+(H(dn_dsmaller)-Hdl(dn_dsmaller)).*vdl(dn_dsmaller);
    
    diagur_flux_grad = (flux_in-flux_out)./(sqrt(2)*Dh)^2;
    
    vlr=(v_diaglr+circshift(v_diaglr,[-1 -1]))./2;
    vul=(v_diaglr+circshift(v_diaglr,[1 1]))./2;
    flux_in=zeros(n,m);
    flux_out=zeros(n,m);
    upd_bigger = find(Hdr>H);
    upd_smaller = find(Hdr<H);
    dnd_bigger = find(Hul>H);
    dnd_smaller = find(Hul<H);
    flux_in(upd_bigger)=flux_in(upd_bigger)+(Hdr(upd_bigger)-H(upd_bigger)).*vlr(upd_bigger);
    flux_out(upd_smaller)=flux_out(upd_smaller)+(H(upd_smaller)-Hdr(upd_smaller)).*vlr(upd_smaller);
    flux_in(dnd_bigger)=flux_in(dnd_bigger)+(Hul(dnd_bigger)-H(dnd_bigger)).*vul(dnd_bigger);
    flux_out(dnd_smaller)=flux_out(dnd_smaller)+(H(dnd_smaller)-Hul(dnd_smaller)).*vul(dnd_smaller);
    
    diagdr_flux_grad = (flux_in-flux_out)./(sqrt(2)*Dh)^2;
    
    advection = zeros(n,n);
    advection = x_flux_grad+y_flux_grad+diagur_flux_grad+diagdr_flux_grad;
    
    
    if(nh==1)
        Neigh_sum = H;
    else
        Neighborhood = ones(nh,nh);
        Neigh_sum = convolve2(H,Neighborhood,'circular');
    end
    Cq = 1-C.*Neigh_sum/(nh^2);
    Cq(find(Cq<.0001))=.0001;
    
    
    Hlog = H.*(1-H);
    Qlog = beta.*Q.*(1-Q./Cq);
    
    LapQ = (Qr+Ql+Qu+Qd-4*Q)./Dh^2;
    LapH = (Hr+Hl+Hu+Hd-4*H)./Dh^2;
    diffH = w*advection+(1-w)*LapH+alpha.*LapQ;
    diffQ = gamma.*LapQ;
    
    %diffH(find(H>=1&diffH>0))=0;
    
    %nextH = H + Dt*(Hlog + advection + D_hq.*LapQ);
    nextH = H + Dt*(Hlog + diffH);
    nextQ = Q + Dt*(Qlog + diffQ);

    
    
    H = nextH;
    Q = nextQ;
    
    H(H<0)=0;
    Q(Q<0)=0;
    H(H>cap)=cap;
    Q(Q>cap)=cap;
    
    
    if(length(find(H<0)))
        disp('neg H')
        return
    end
    
    if round(nansum(diffQ,'all'),0)
        disp('diffQ not conserved')
        return
    end
    
    
    store_Htot(t) = sum(sum(H))/(n*m);
    store_Qtot(t) = sum(sum(Q))/(n*m);
    
    if(~mod(t,200))
        store_H(:,:,t/200) = H;
        store_Q(:,:,t/200) = Q;
    end
    
    if(~mod(t,200))
        if(1)
        
        %         figure(1)
        %         imagesc(advection)
        %         title('advection')
        %         colorbar
        %         figure(2)
        %         imagesc(diffH)
        %         colorbar
        %         title('diffH')
        %         colorbar
        %         figure(3)
        %         imagesc(Hlog)
        %         title('Hlog')
        %         colorbar
        %         figure(12)
        %         imagesc(diffQ)
        %         colorbar
        %         title('diffQ')
        %         colorbar
        %         figure(13)
        %         imagesc(Qlog)
        %         title('Qlog')
        %         colorbar
        
                 %store_H(:,:,t/1) = H;
                 %store_Q(:,:,t/1) = Q;
        %
        figure(6)       
        subplot(2,2,1);
        imagesc(H)
        colorbar     
        title('H')     
        subplot(2,2,2);       
        imagesc(Q)       
        colorbar
        title('Q')    
        subplot(2,2,[3,4]);        
        plot(store_Htot(1:t))        
        hold on       
        plot(store_Qtot(1:t))      
        hold off
        
        
        xlim([0,time_steps]);
        
        xlabel('time','FontSize',15)
        legend({'\langle H \rangle','\langle Q \rangle'},'FontSize',12)
        legend boxoff
        
        %pause(0.1)
        %         title(['n_{h} = ',num2str(nh),',  S = ',num2str(round(S,3)),',  D = ',num2str(D),'  R= ',num2str(R)])
        
    else
        figure(1)
        subplot(2,1,1)
        imagesc(H)
        axis equal
        axis tight
        colorbar
        subplot(2,1,2)
        imagesc(Q)
        axis equal
        axis tight
        colorbar
        pause(0.1)
        end
    end
    
    
end

