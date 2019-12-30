function dXdt = CurrentClampIHCmodel(t,X)
global Te Stim IHC;
 
%Stimulation
I = Stim(round(t/Te)+1);      %pA
 
%Equation diff?rentiel
dXdt = zeros(5,1) ;
dXdt(1) = (1/IHC.Cm)*...                                            %Courant de membrane
        (- IHC.Ileak.g *(X(1)-IHC.Ileak.E) ...                        %Courant de fuite
         - IHC.IKF.g*(X(2)^4)*(X(1)-IHC.IKF.E)...                   %Courant potassique rapide
         -IHC.IKS.g*(X(3)^4)*(X(1)-IHC.IKS.E) ...                   %Courant potassique lent
         -IHC.IKN.g*(X(4)^4)*(X(1)-IHC.IKN.E) ...                   %Courant potassique n
         - IHC.ICA.g*(X(5)^3)*(X(1)-IHC.ICA.E) ...
         + I);                                                      %Courant ext?rieur

%X(2)                                                               %Valeur initial de la variable d'activation de Ikf       
fi = sigmoide(IHC.IKF.xinf,X(1)) ;
tf = cinetique2(IHC.IKF.taux,X(1)) ;
dXdt(2) = (fi-X(2))/tf ;
%X(3)                                                               %Valeur initial de la variable d'activation de Iks  
si = sigmoide(IHC.IKS.xinf,X(1)) ;
ts = cinetique2(IHC.IKS.taux,X(1)) ;
dXdt(3) = (si-X(3))/ts ;
%X(4)                                                               %Valeur initial de la variable d'activation de Ikn  
ni = sigmoide(IHC.IKN.xinf,X(1)) ;
tn = cinetique2(IHC.IKN.taux,X(1)) ;
dXdt(4) = (ni-X(4))/tn ;
%X(5)                                                               %Valeur initial de la variable d'activation de ICa  
ci = sigmoide(IHC.ICA.xinf,X(1)) ;
tc = cinetique2(IHC.ICA.taux,X(1)) ;
dXdt(5) = (ci-X(5))/tc ;

return
 
end

