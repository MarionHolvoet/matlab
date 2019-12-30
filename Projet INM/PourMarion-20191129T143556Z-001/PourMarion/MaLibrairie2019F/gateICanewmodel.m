function dxdt = gateICanewmodel(t,x)

global IHC stim STIM 

indice=floor(t*stim.Fe)+1;

if length(x)>1
    V0=STIM(indice,:)';
else
    V0=STIM(indice);
end

xinf=sigmoide(IHC.ICA.xinf, V0);
tau=cinetique2(IHC.ICA.taux,V0);

dxdt = (xinf-x)./tau ;

