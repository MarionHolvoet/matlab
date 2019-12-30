function err = testIk(par)

global IHC IHCm stim titre t STIM t_IV window_IV MonBesselChoupiParam macolonne current;

f0=sigmoide(IHC.IKF.xinf, stim.VMAINTIEN*ones(size(stim.VSTEP)))
s0=sigmoide(IHC.IKS.xinf, stim.VMAINTIEN*ones(size(stim.VSTEP)))
n0=sigmoide(IHC.IKN.xinf, stim.VMAINTIEN*ones(size(stim.VSTEP)))
c0=sigmoide(IHC.ICA.xinf, stim.VMAINTIEN*ones(size(stim.VSTEP)))

%params = par;

[t,f] = ode15s('gateIKFnewmodel',t,f0) ;
[t,s] = ode15s('gateIKSnewmodel',t,s0) ;
[t,n] = ode15s('gateIKNnewmodel',t,n0) ;
[t,c] = ode15s('gateICanewmodel',t,c0) ;

Im = par(8)*(f.^4).*(STIM-par(9))+...  %IKF
    par(17)*(s.^4).*(STIM-par(18))-... %IKS
    par(26)*(c.^3).*(STIM-par(27)) + ...   %ICA
    par(35)*(n.^4).*(STIM-par(36))+ ...   %IKN
    par(37).*(STIM-par(38));              % Ifuite

if MonBesselChoupiParam.is==1;
     [b,a]=MonBesselChoupi(MonBesselChoupiParam.ordre1, 2*pi*MonBesselChoupiParam.fc1, stim.Fe)
     Im=filter(b,a,Im);
    [b,a]=MonBesselChoupi(MonBesselChoupiParam.ordre2, 2*pi*MonBesselChoupiParam.fc2, stim.Fe)
    Im=filter(b,a,Im);
end

indices=find(t<=t_IV+window_IV/2  & t>=t_IV-window_IV/2);

IVexp=mean(current(indices, :));
IVmod=mean(Im(indices, :));

% t_0 = 0.030
% t_f = 0.033
%whos current I IVexp IVmod
residu=current-Im;

err = var(residu(:)); 
%


Vm=[-70:5:30]*1e-3;
Vm=Vm(:,macolonne);
%% fig generale :

figure(2) ;%'Position',[20 60 900 1200], 'PaperPositionMode', 'auto');

subplot(4,4,[1 2 5 6]);
plot(t,current,'k','linewidth',0.4) ; hold on ;
plot(t,Im,'r','linewidth',1.2) ;
%ylim([-0.6 1.2]);
%xlim([0.240 0.280]);
set(gca, 'tickdir', 'out');
title(titre)

plot(t_IV*[1 1],[min(min(current)) max(max(current))],'b'); hold off ;
 
drawnow ;

subplot(4,4,[3 7]);
plot(t,current,'k','linewidth',1.2);
%ylim([-0.5 0.5]);
%xlim([0.05 0.3]);

subplot(4,4,[4 8]);
plot(t,Im,'r','linewidth',1.2) ;
%ylim([-0.5 0.5]);
%xlim([0.05 0.3]);

hold off ;
V = stim.VSTEP ;   % ?????????????????????? = -50 tt le temps !


subplot(4,4,9);
%indx=find(t>=0.10&t<=0.105);
%IV=mean(current(indx,:))
% plot(Vm,IVexp,'o-k')
% title('IV - curve')
% hold on;
% %IVm=mean(Im(indx,:))
% plot(Vm,IVmod,'o-r')
% hold off;
% % 
% subplot(4,4,10);
% G=IVexp/max(IVexp);
% plot(Vm,G,'k')
% hold on;
% Gm=IVmod/max(IVmod);
% plot(Vm,Gm,'r')
% 
% subplot(4,4,11);
% x_kf=sigmoide(IHC.IKF.xinf,Vm);
% x_ks=sigmoide(IHC.IKS.xinf,Vm);
% x_kn=sigmoide(IHC.IKN.xinf,Vm);
% x_ca=sigmoide(IHC.ICA.xinf,Vm);

IHCm.IKF.xinf.demie=par(1); %mV
IHCm.IKF.xinf.k=par(2);  %mV
IHCm.IKF.taux.amp=par(3); %sec
IHCm.IKF.taux.offset=par(4); %sec
IHCm.IKF.taux.V=par(5); %mV
IHCm.IKF.taux.k1=par(6); %mV
IHCm.IKF.taux.k2=par(7); %mV
IHCm.IKF.g=par(8); %S
IHCm.IKF.E=par(9); %mV

IHCm.IKS.xinf.demie=par(10); %mV
IHCm.IKS.xinf.k=par(11);  %mV
IHCm.IKS.taux.amp=par(12); %sec
IHCm.IKS.taux.offset=par(13); %sec
IHCm.IKS.taux.V=par(14); %mV %mV
IHCm.IKS.taux.k1=par(15); %mV
IHCm.IKS.taux.k2=par(16); %mV
IHCm.IKS.g=par(17); %S
IHCm.IKS.E=par(18); %mV

IHCm.ICA.xinf.demie=par(19); %mV
IHCm.ICA.xinf.k=par(20);  %mV
IHCm.ICA.taux.amp=par(21); %sec
IHCm.ICA.taux.offset=par(22); %sec
IHCm.ICA.taux.V=par(23); %mV
IHCm.ICA.taux.k1=par(24); %mV
IHCm.ICA.taux.k2=par(25); %mV
IHCm.ICA.g=par(26); %S
IHCm.ICA.E=par(27); %mV

IHCm.IKN.xinf.demie=par(28); %mV
IHCm.IKN.xinf.k=par(29);  %mV
IHCm.IKN.taux.amp=par(30); %sec
IHCm.IKN.taux.offset=par(31); %sec
IHCm.IKN.taux.V=par(32); %mV
IHCm.IKN.taux.k1=par(33); %mV
IHCm.IKN.taux.k2=par(34); %mV
IHCm.IKN.g=par(35); %S
IHCm.IKN.E=par(36); %mV

IHCm.Ileak.g=par(37);
IHCm.Ileak.E=par(38);

idx=find(t>=0.007&t<=0.01502);

IHCm.Im=Im*1e9;
IHCm.current=current*1e9;
IHCm.t=t;
IHCm.STIM=STIM*1e3;
IHCm.stim=stim;


%now=[IHCm.IKF.g*1e9, IHCm.IKF.E, IHCm.IKS.g*1e9, IHCm.IKS.E, IHCm.IKS.xinf.demie, IHCm.IKS.taux.V]; 


%% Les Gcourb

% xm_kn=sigmoide(IHCm.IKN.xinf,Vm);
% plot(Vm,x_kf,'r');
% hold on;
% plot(Vm,x_ks,'b');
% plot(Vm,x_kn,'g');

% moy=load('tails.mat')
% plot(Vm(1:10),xm_kn(1:10), 'r');
% hold on;
%plot(Vm(1:10),moy.courb(1:10),'k');
%plot(V,x_ca,'k');


% %% les TauX
% subplot(4,4,12);
%  taum_kf=cinetique2(IHCm.IKF.taux,Vm)
% plot(Vm,taum_kf,'r')
% % hold on;
% % plot(Vm,moy.tauM(1:12),'k')
% % hold off;
% 
% %% Les G individuels :
% gN=(IHCm.IKN.g)*n.^4;
% gF=(IHCm.IKF.g)*f.^4;
% gS=(IHCm.IKS.g)*s.^4;
% gCa=(IHCm.ICA.g)*c.^3;
% 
% subplot(4,4,13);
% plot(t,gF*1e9,'r')
% %pause
% subplot(4,4,14);
% plot(t,gS*1e9,'b')
% subplot(4,4,15);
% plot(t,gN*1e9,'k')
% subplot(4,4,16);
% plot(t,gCa*1e9,'Color',[0.9 0.1 0])
% 
% 

%%
save(['IHCm_mod_STIM_m201216_2','.mat'], '-struct', 'IHCm');
%print('-dpng', '-r300', [list(1:9),list(end-6:end-4),'.png']);

%%
%if maxfunevals==MaxFun
    


%end

