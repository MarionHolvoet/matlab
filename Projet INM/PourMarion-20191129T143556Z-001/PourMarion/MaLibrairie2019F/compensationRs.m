function commandeComp=compensationRs(I, commande, Rs, compensation)

% Usage commandeComp=compensationRs(commande, Rs, compensation)

dV=I*Rs*(100-compensation)*0.01;
commandeComp=commande-dV;