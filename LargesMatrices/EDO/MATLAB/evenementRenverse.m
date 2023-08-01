function [valeur, fini, direction] = evenementRenverse(t,X)
% identifie le moment où l'on croise un bord
% et signale au solveur de s'arrêter
%valeur = min([abs(X(2)), abs(X(2)-1)]) ; 
valeur = X(2)*(X(2)-1);
fini   = 1 ; % donne au solveur l'information qu'il doit s'arrêter
direction = +1 ; % c'est exactement pareil que la fonction evenement
                % mais je la remets pour signaler que ce paramètre
                % prend la même valeur, puisque le champ de vecteurs
                % que l'on suit est l'opposé, donc on rencontre le 
                % bord comme avant
                
end