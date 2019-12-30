function x=sigmoide(parametres, Vm)
%USAGE : x=sigmoide(parametres, Vm)
% Vm : potentiel de membrane
% parametres : structure de param�tres
% x : la valeur de la sigmoide au(x) potentiel(s) Vm

x=1./(1+exp(-(Vm-parametres.demie)/parametres.k));
