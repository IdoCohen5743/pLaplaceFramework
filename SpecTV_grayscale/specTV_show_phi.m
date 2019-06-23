function I = specTV_show_phi( Phi,dt, Ind, RbyC, fignum, Rescale, Mask)
% Private function by Guy Gilboa (June 2013)
% Shows Phi*dt channels
% Phi = Phi(x,y;t) function
% dt = evolution time step
% Ind = selected indices in time
% RbyC = Rob by Column matrix size
% fignum = figure number to show image
% Rescale = rescale image size (optional)
% Mask = Masking parts of Phi (optional)
% example: I = specTV_show_color_phi( Phi,Ind, RbyC, fignum)

if (~exist('Mask','var'))
    Mask=ones(size(Phi,1),size(Phi,2));
end
if (~exist('Rescale','var'))
    Rescale = 1;
end
% show combined phi's of color images on screen 
Nr = RbyC(1);
Nc = RbyC(2);
Rs = Rescale;

phi = imresize(Phi(:,:,1),Rs,'nearest');
Ny = size(phi,1);
Nx = size(phi,2);

Phi = Phi*dt; % rescale by dt to standard image scale

figure(fignum);

I = zeros(Nr*Ny+(Nr+1)*2,Nc*Nx+(Nc+1)*2);
cnt = 1;
for j=1:Nr,
    for i=1:Nc,
        phi = imresize(Mask.*Phi(:,:,Ind(cnt)),Rs,'nearest');
        phi = 1.5*phi+0.5;
        I(j*(Ny+2)-Ny+1:j*(Ny+2),i*(Nx+2)-Nx+1:i*(Nx+2))=phi;
        cnt = cnt+1;
    end % for i
end % for j
imshow(I);

end

