function dxdt = gateIKSnewmodel(t,x)

global IHC stim STIM

indice=floor(t*stim.Fe)+1;
%indice=indice-indice(1)+1

if length(x)>1
    V0=STIM(indice,:)';
else
    V0=STIM(indice);
end

xinf=sigmoide(IHC.IKS.xinf, V0);
tau=cinetique2(IHC.IKS.taux,V0);

dxdt = (xinf-x)./tau ;

