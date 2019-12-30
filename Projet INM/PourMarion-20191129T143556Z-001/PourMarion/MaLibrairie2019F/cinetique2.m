function tau=cinetique(parametres,Vm)
    tau=(parametres.amp)./(exp((Vm-parametres.V)/parametres.k1) ...
        + exp(-(Vm-parametres.V)/parametres.k2))+parametres.offset;
    
    
    
    