%%%%%%%%%%%%%%%%%%%%%%%%%%% Divergence (backward difference), adjoint BC
%% Chombolle's implementation 
%% See "An Algorithm for Total Variation Minimization and Applications". So help you God
function fd = div(Px,Py)

    fx = [Px(:,1:end-1) zeros(size(Px,1),1)] - [zeros(size(Px,1),1) Px(:,1:end-1)]; 
    fy = [Py(1:end-1,:); zeros(1,size(Py,2))] - [zeros(1,size(Py,2)); Py(1:end-1,:)];
    fd = fx+fy;


