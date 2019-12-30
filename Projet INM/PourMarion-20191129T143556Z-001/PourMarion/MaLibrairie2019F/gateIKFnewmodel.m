function dxdt = gateIKFnewmodel(t,x)

global IHC stim STIM 

indice=floor(t*stim.Fe)+1; %-1073;
%indice=indice-indice(1)+1;     % 475   (1) 296

if length(x)>1
    V0=STIM(indice,:)';
    
else
    V0=STIM(indice);
    %disp(V0)
    
end

xinf=sigmoide(IHC.IKF.xinf, V0);
tau=cinetique2(IHC.IKF.taux,V0);

dxdt = (xinf-x)./tau ;

