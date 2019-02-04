%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% p-Laplace operator
%% according to "Stable Explicit p-Laplacian Flows based on Nonlinear Eigenvalue Analysis"

function Delta_pu = lapi(u,p)

[grad_x,grad_y] = grad(u);
theta = atan2(grad_y,grad_x);
cosTH = cos(theta);
sinTH = sin(theta);

grad_x2 = grad_x.^2;
grad_y2 = grad_y.^2;
aGrad = sqrt(grad_x2 + grad_y2);
aGradP = aGrad.^(p-1);

Delta_pu = div(cosTH.*aGradP,sinTH.*aGradP);
