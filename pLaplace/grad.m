%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient (forward difference)
%% Chombolle's implementation 
%% See "An Algorithm for Total Variation Minimization and Applications". So help you God
function [fx,fy] = grad(P)
    fx = P(:,[2:end end])-P;
    fy = P([2:end end],:)-P;

