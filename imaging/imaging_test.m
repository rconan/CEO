%% 
% gpuDevice([])

 r0 = 15e-2;
 L0 = 30;
 atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

nLenslet = 4;
nPxLenslet = 16;
nxy = nLenslet*nPxLenslet;
ceodir = '~/CEO';
% cd([ceodir,'/imaging'])
% unix(['sed -i ',...
%     '-e ''s/#define _N_LAYER_ [0-9]*/#define _N_LAYER_ ',num2str(atm.nLayer),'/g'' ',...
%     '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
%     num2str(nxy^2),'/g'' definitions.h']);
% unix('cat definitions.h');
cd(ceodir)
unix('make clean lib imaging.mex')
cd([ceodir,'/imaging'])
clear ceo_imaging
mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu

D = 0.4;
u = single( 0.5*D*gpuArray.linspace(-1,1,nxy) );
[x,y] = meshgrid( u );
%%
[phs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
figure(1)
clf
ax1 = axes('pos',[0.1300    0.6149    0.6589    0.3412]);
imagesc(u,u,phs)
axis square
colorbar
ax2 = axes('pos',[0.1300    0.0634    0.6589    0.3412]);
imagesc(frame)
axis square
colorbar
ax3 = axes('pos',[0.1000    0.32    0.3643    0.4306]);
imagesc([reshape(cx-7.5,4,4),reshape(cy-7.5,4,4)])
axis equal tight
colorbar('location','NorthOutside')
ax4 = axes('pos',[ 0.5979    0.445   0.1879    0.1799]);
imagesc(reshape(flux/256,4,4))
axis equal tight
colorbar('location','NorthOutside')

%%
nSample = 2500;
Gcc = gpuArray.zeros(nLenslet^2*2,'single');
nP = nxy;%3*nLenslet+1;
Gpc = gpuArray.zeros( nP^2 , nLenslet^2*2, 'single' );
h = waitbar(0,'Building covariance matrices ...!');
ui = linspace(1,nxy,nP);
[xi,yi] = meshgrid( ui );
for kSammple=1:nSample
    [phs,~,cx,cy,~] = ceo_imaging(x,y,1,L0,0);
    cx = cx - 7.5;
    cy = cy - 7.5;
    c = [cx;cy];
    Gcc = Gcc + c*c';
%     phsi = interp2(gather(phs),xi,yi);
    Gpc = Gpc + phs(:)*c';
    waitbar(kSammple/nSample,h)
end
close(h)
Gcc = Gcc./nSample;
Gpc = Gpc./nSample;
figure(22)
ha1 = subplot(1,2,1);
imagesc(Gcc)
axis equal tight
colorbar('location','NorthOutside')
ha2 = subplot(1,2,2);
imagesc(Gpc)
colorbar('location','NorthOutside')

%%
[phs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
cx = cx - 7.5;
cy = cy - 7.5;
c = [cx;cy];
M = Gpc/Gcc;
phse = M*c;
phase2nm = 500/2/pi;
phs_zm = phase2nm*( phs-mean(phs(:)) );
phse_zm = phase2nm*( reshape(phse-mean(phse),nP,nP) );
phs_err = phs_zm - phse_zm;
figure(23)
subplot(2,1,1)
imagesc([ phs_zm, phse_zm])
title(sprintf('WF rms [nm] : %5.2f',std(phs_zm(:)) ) )
axis equal tight
colorbar
subplot(2,1,2)
imagesc(phs_err)
title(sprintf('WFE rms [nm] : %5.2f',std(phs_err(:)) ) )
axis square
colorbar