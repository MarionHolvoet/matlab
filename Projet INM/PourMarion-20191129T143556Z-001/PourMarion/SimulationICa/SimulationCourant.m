clc; close all;
clear all;

global IHC stim Vstep STIM
addpath('C:\Users\holvo\OneDrive\Documents\Projet INM\PourMarion-20191129T143556Z-001\PourMarion\MaLibrairie2019F')

police=14;
epaisseurTrait=1.2;

%% commande
t0=evalin('base','temps0'); t1=evalin('base','temps1'); t2=evalin('base','temps2'); t3=evalin('base','temps3'); 
stim.T=[t0 t1 t2 t3];
% stim.T=[0 0.2 0.6 0.8]
stim.VREPOS=evalin('base','Vrepos'*1e-3);
stim.VSTEP=evalin('base','Vstim'*1e-3);
stim.Fe=evalin('base','fe');
% stim.VREPOS=[-60]*1e-3;
% stim.VSTEP=[10]*1e-3';
% stim.Fe=200000;
Te=1/stim.Fe;
t=stim.T(1):Te:stim.T(end);
F0=1/stim.T(4);
phase=(stim.T(2)-stim.T(1))*2*pi/stim.T(4);
dutyCycle=100*(stim.T(3)-stim.T(2))/stim.T(4);

%% paramètres du modèle de CCI
IHC.ICA.E=0.0429; % en Volt
IHC.ICA.g=2.6200e-009; % en Siemens
IHC.ICA.xinf.demie=-0.0371; % en Volt
IHC.ICA.xinf.k=0.0140; %  en Volt
IHC.ICA.taux.amp= 9.6000e-004; % en secondes
IHC.ICA.taux.offset= 1.1000e-004; % en secondes
IHC.ICA.taux.V= -0.0381; % en Volt
IHC.ICA.taux.k1= 0.0155; % en Volt
IHC.ICA.taux.k2= 0.0870; % en Volt


%% courbes d'activation et de cinétique
Vm=[-100:5:50]*1e-3;
cinf=sigmoide(IHC.ICA.xinf,Vm); % courbe d'activation
ctau=cinetique2(IHC.ICA.taux, Vm); % cinétique d'activation

%% Affiche des courbes de la courbe d'activation et de cinétique
figure('position',[600 80 1000 400],'paperpositionmode','auto');
axes('position',[0.1 0.17 0.35 0.75]);
set(gca,'fontsize',police,'linewidth',epaisseurTrait)
hold on;
plot(Vm*1000, cinf.^3,'b','linewidth',epaisseurTrait+0.3);
set(gca,'tickdir','out','ticklength',[0.02 0.02])
xlabel('V_m (mV)','fontsize',police+2)
ylabel('c_{oo} ','fontsize',police+2);
grid on;

axes('position',[0.6 0.17 0.35 0.75]);
set(gca,'fontsize',police,'linewidth',epaisseurTrait)

hold on;
plot(Vm*1000, 1000*ctau,'b','linewidth',1.5);
set(gca,'tickdir','out','ticklength',[0.02 0.02])
xlabel('V_m (mV)','fontsize',police+2);
ylabel('\tau (ms)','fontsize',police+2);
grid on;
%print('-dpng','-r300','CaracteristiqueDuModele');


positionV=[0.15 0.8 0.8 0.15];
macouleurV=[0 0 0];
positionI=[0.15 0.1 0.8 0.6];
for n=1:length(stim.VSTEP)
    STIM(:,n)=(stim.VSTEP(n)-stim.VREPOS)*(0.5*square(2*pi*F0*t-phase,dutyCycle)+0.5)+stim.VREPOS;
end


c0=sigmoide(IHC.ICA.xinf, stim.VREPOS*ones(size(stim.VSTEP)));
[t,c] = ode15s('gateICanewmodel',t,c0) ;

ICa=IHC.ICA.g*(c.^3).*(STIM-IHC.ICA.E);


figure('position',[600 80 600 800],'paperpositionmode','auto');
axes('position',[0.17 0.75 0.75 0.2]);
set(gca,'fontsize',police,'linewidth',epaisseurTrait);
hold on;
plot(t*1000,1000*STIM,'k','linewidth',epaisseurTrait+0.3);
set(gca,'tickdir','out','ticklength',[0.02 0.02])
xlabel('Time (ms)','fontsize',police+2);
ylabel('V_m (mV)','fontsize',police+2);
ylim(1000*[min(stim.VREPOS)-0.01 max(stim.VSTEP)+0.01]);
grid on

axes('position',[0.17 0.08 0.75 0.55]);
set(gca,'fontsize',police,'linewidth',epaisseurTrait);
hold on;
plot(t*1000, ICa*1e12,'b','linewidth',2);
xlim([stim.T(1) stim.T(4)]*1000)
ylabel(' Current (pA)','fontsize',police+2);
xlabel(' Time (ms)','fontsize',police+2);
set(gca,'tickdir','out','ticklength',[0.02 0.02])
grid on
box off
ylim([-300 50]);

%print('-dpng','-r300','VoltClamp');

