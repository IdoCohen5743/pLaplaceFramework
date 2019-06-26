% Decompose textures in image with spectral TV, gray-scale images
% Script by Guy Gilboa (Jan 2015).
% Based on: [1] G. Gilboa, "A total variation spectral framework for scale and texture analysis." SIAM Journal on Imaging Sciences 7.4 (2014): 1937-1961.

close all;
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Important params to change according to image / application %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Max_time = 4.2;       % Maximal scale to be processed (in evolution time)
%Max_time = 40;     % large scale
Num_of_bands = 30000;  % Number of bands phi(x,t);
%%

dt = Max_time/Num_of_bands; 

f = double(rgb2gray(imread('zebra_media_gmu.jpg'))); 
%f = f(50:100,50:100); % debug
% f = double(f);
% f = double(rgb2gray(f));
f = imresize(f,0.16);
% f = f(21:185,71:240);  % fruits, melon
f = f - mean(f(:));
f = double(f)/double(max(f(:)));
% f = f/255;  % pixels are in the range [0,1]
% f = f- mean(f(:));
figure(1); imshow(f,[]); title('f')

% Compute Phi bands and residual f_r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,T,Phi,f_r] = specTV_evolve(f, Max_time, dt);  % evolve image
% Define Filters
factor = 1;
hp_i = (Num_of_bands*1.5/100);
bp_low_i=(Num_of_bands*7.5/100);
bp_high_i = (Num_of_bands*20/100); 
H1 = zeros(size(T));
H1(1:hp_i)=1;  % high pass, melon
H2 = zeros(size(T));
H2(hp_i:bp_low_i)=1;  % band pass, melon
H3 = zeros(size(T));
H3(bp_low_i:bp_high_i)=1;  % band pass, melon
H4 = zeros(size(T));
H4(bp_high_i:end)=1;  % band pass, melon

%% Reconstruct filtered image given Phi, H and f_r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_H1 = specTV_filter( Phi, H1, f_r, dt );
f_H2 = specTV_filter( Phi, H2, f_r, dt );
f_H3 = specTV_filter( Phi, H3, f_r, dt );
f_H4 = specTV_filter( Phi, H4, f_r, dt );
% Plot S(t) and filter response
%%
h = figure(); plot(T,S,...
                T(hp_i:bp_low_i),S(hp_i:bp_low_i),...
                T(bp_low_i:bp_high_i),S(bp_low_i:bp_high_i),...
                T(bp_high_i:end),S(bp_high_i:end));hold off;
for iii=1:1:length(h.Children.Children)
    h.Children.Children(iii).LineWidth = 8;
end
grid on;
h.Children.FontSize = 70;
h.Children.TickLabelInterpreter = 'Latex';
h.Children.YLim = [0,1.1*max(S)];
% h.Children.YLim = [0,1000*max(S)];
% h.Children.YScale = 'log';
h.Children.XLim = [0 T(end)];
h.Children.XTick = [0:ceil((T(end)+0.1)/5):T(end)+0.1];
% h.Children.YTick = 10.^[floor(-log(max(S))-1):ceil((log(max(S)))/7):log(max(S))+1];
% h.Children.XTick = [0:ceil(4.2/5):4.2];
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
drawnow
ax_s=gca; outerpos = ax_s.OuterPosition;
ti = ax_s.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax_s.Position = [left bottom ax_width ax_height];
h = legend(['$0\%-1.5\%$'],['$1.5\%-7.5\%$'],['$7.5\%-20\%$'],['$20\%-100\%$']);
h.Interpreter = 'Latex';            
% xlabel('t'); 
% legend('Spectrum S(t)','High Pass Filter H1(t)', 'Band Pass Filter H2(t)')
% show a few Phi(t) instances
%phi_show = [4 20 30 50];
%%
h=figure(); 
imshow(f_H1,[]);
[row,col] = size(f_H1);
h.InnerPosition(3)=col;
    drawnow;
    %     h.InnerPosition(4)=col;
    ha = gca;
    h.InnerPosition(4)=floor(h.InnerPosition(3)*row/col);
    set(ha,'position',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'OuterPosition',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'Units','pixels');
    %     pos=get(ha,'position');
    h = gcf;
    set(h,'color','w');
h = figure(); 
imshow(f_H2,[]);
h.InnerPosition(3)=col;
    drawnow;
    %     h.InnerPosition(4)=col;
    ha = gca;
    h.InnerPosition(4)=floor(h.InnerPosition(3)*row/col);
    set(ha,'position',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'OuterPosition',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'Units','pixels');
    %     pos=get(ha,'position');
    h = gcf;
    set(h,'color','w');
h= figure(); 
imshow(f_H3,[]);
h.InnerPosition(3)=col;
    drawnow;
    %     h.InnerPosition(4)=col;
    ha = gca;
    h.InnerPosition(4)=floor(h.InnerPosition(3)*row/col);
    set(ha,'position',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'OuterPosition',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'Units','pixels');
    %     pos=get(ha,'position');
    h = gcf;
    set(h,'color','w');
h=figure(); 
imshow(f_H4,[]);
h.InnerPosition(3)=col;
    drawnow;
    %     h.InnerPosition(4)=col;
    ha = gca;
    h.InnerPosition(4)=floor(h.InnerPosition(3)*row/col);
    set(ha,'position',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'OuterPosition',[0 0 1 1]); %[left bottom width height]
    %     set(ha,'Units','pixels');
    %     pos=get(ha,'position');
    h = gcf;
    set(h,'color','w');
% figure(6); 
% imshow(f_r); title('Residual f_r (Low-pass)')

% Ind = 1:10;    % Index of phi bands
% RbyC = [2 5];  % Row by Column boxes
% fignum =10;    % Num of figure
% Rescale = 1;   % Possible image size rescaling (1 - no rescale, >1 larger images)
% contrast = 5;
% I = specTV_show_phi( Phi*contrast,dt, Ind, RbyC, fignum, Rescale);
% title(['A few instances of \phi(x;t), from bands ' num2str(Ind(1)) ' to ' num2str(Ind(end))])
