%% 
% gpuDevice([])

nLenslet_ = [20 40 64 84 150];
D_        = [3.6 8 5 42 30];
    
r0 = 15e-2;
L0 = 30;
atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
phase2nm = 1e9;%*lambda/2/pi;
    
D = 42;

%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

nLenslet = 84;
d = D/nLenslet;
nPxLenslet = 16;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;

ne = nLenslet+1;

compile = true;

%%
ceodir = '~/CEO';
if compile
cd([ceodir,'/include'])
unix(['sed -i ',...
    '-e ''s/#define N_SIDE_LENSLET [0-9]*/#define N_SIDE_LENSLET ',num2str(nLenslet),'/g'' ',...
    '-e ''s/#define _N_PX_PUPIL_ [0-9]*/#define _N_PX_PUPIL_ ',num2str(nPxLenslet),'/g'' ',...
    '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
    num2str((nLenslet*nPxLenslet)^2),'/g'' definitions.h']);
unix('cat definitions.h');
cd(ceodir)
unix('make clean all')
% unix('make imaging.mex')
% clear ceo_imaging
% mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu
end
%%
nIt = 100;
cd([ceodir,'/iterativeSolvers'])
unix('make iterativeSolvers.bin')
% fprintf(' ==>>> CG (N=%d)\n',nLenslet)
% tic
% unix(sprintf('./a.out %3.1f CG > CG_%03d_%03d.log',D,nIt,nLenslet));
% toc
%%
fprintf(' ==>>> MINRES (N=%d)\n',nLenslet)
tic
unix(sprintf('./a.out %3.1f MINRES > MINRES_%03d_%03d.log',D,nIt,nLenslet));
toc
%%
% ps = loadBin('phaseScreen',[nxy,nxy]);
% pup = logical( loadBin('A_mask','map','type','char') );
pup = tools.piston(ne,'type','logical','shape','square');
ps = loadBin('phaseScreenLowRes','map','scaling',1e9);
ps = pup.*(ps - mean(ps(pup(:))));
figure(102)
% subplot(2,3,[1,2])
% imagesc(ps)
% axis square
% xlabel(colorbar('location','northOutside'),'[nm]')

c = loadBin('centroids',[nLenslet,nLenslet*2]);
subplot(3,4,[9,12])
imagesc(c)
axis equal tight
colorbar('location','north')

%ps_e = phase2nm*loadBin('phaseScreenEst',[ne,ne]);
%ps_e = phase2nm*loadBin(sprintf('CG_phaseEst_%3d_%2d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = phase2nm*loadBin(sprintf('MINRES_phaseEst_%d',nIt),'map');
ps_e = pup.*bsxfun( @minus, ps_e, mean(ps_e(pup(:)),1) );
ps_e_k = reshape(ps_e,[ne,ne]);
subplot(3,4,[1,6])
imagesc([ps,ps_e_k])
axis equal tight
xlabel(colorbar('location','northOutside'),'[nm]')

ps_err = bsxfun( @minus, ps(:), ps_e(:));
rms_ps_err = std(ps_err);
ps_err_k = reshape(ps_err,[ne,ne]);
subplot(3,4,[3,8])
imagesc(ps_err_k)
axis equal tight
title(sprintf('wfe=%6.2fnm',rms_ps_err))
xlabel(colorbar('location','southOutside'),'[nm]')

