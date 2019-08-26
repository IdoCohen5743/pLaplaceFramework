function [S,T,Phi,f_r] = specTV_evolve(f, Max_time, dt, Method)
% private function by Guy Gilboa (Jan 2015)
% Computes TV spectrum of an image, returns S, time interval T, Phi(T)
% and residual image f_r (to be added in the reconstruction).
% Method is an optional struct specifying numerical method and params
% Example: [S,T,Phi] = specTV_evolve(f, Max_time, dt)
% Based on: [1] G. Gilboa, "A total variation spectral framework for scale and texture analysis." SIAM Journal on Imaging Sciences 7.4 (2014): 1937-1961.


if exist('Method','var')
    Num_method = Method.Num_method;
    dt_proj = Method.dt_proj;
    iter_proj = Method.iter_proj;
else
    Num_method = 'proj'; % Projection algorithm, default (the only method implemented here)       
    dt_proj=0.2; iter_proj=20000; %500
end    

mu = 1/(2*dt);
NumIter = round(Max_time/dt);
scale = 5; % for numerical reasons

S = zeros(1,NumIter); 
Phi = zeros(size(f,1),size(f,2),NumIter);
T = (1:NumIter)*dt;

u0 = f*scale;
u1=proj_tvl2(u0,mu,iter_proj,dt_proj);
u2=proj_tvl2(u1,mu,iter_proj,dt_proj);

for i=1:NumIter,
    ddu = (u0+u2-2*u1)/(dt*dt);  % one/two more iter
    t = i*dt;
    phi = ddu*t;
    Phi(:,:,i) = phi;
    S(i) = sum(abs(phi(:)));
    if (i<NumIter) % not last iteration
        u0=u1;
        u1=u2;
        u2=proj_tvl2(u2,mu,iter_proj,dt_proj);
    end

end % for i

f_r = (NumIter+1)*u1-NumIter*u2;  % residual image

% rescale
S = S/scale; Phi = Phi/scale; f_r = f_r/scale;

end

