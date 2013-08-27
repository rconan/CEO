%% 
% gpuDevice([])

 r0 = 15e-2;
 L0 = 30;
 atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
D = 25;
phase2nm = 1e9*lambda/2/pi;

%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

nLenslet = 60;
d = D/nLenslet;
nPxLenslet = 8;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;
clear ceo_imaging
ceodir = '~/CEO';
cd([ceodir,'/imaging'])
unix(['sed -i ',...
    '-e ''s/#define N_SIDE_LENSLET [0-9]*/#define N_SIDE_LENSLET ',num2str(nLenslet),'/g'' ',...
    '-e ''s/#define _N_PX_PUPIL_ [0-9]*/#define _N_PX_PUPIL_ ',num2str(nPxLenslet),'/g'' ',...
    '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
    num2str((nLenslet*nPxLenslet)^2),'/g'' definitions.h']);
unix('cat definitions.h');
cd(ceodir)
unix('make clean lib imaging.mex')
cd([ceodir,'/imaging'])
clear ceo_imaging
mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu

u = single( 0.5*D*gpuArray.linspace(-1,1,nxy) );
[x,y] = meshgrid( u );

[phs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
figure(1)
subplot(7,3,[1,9])
imagesc(u,u,phs)
axis square
colorbar
subplot(7,3,[13,21])
imagesc(frame)
axis square
colorbar
subplot(7,3,[10,11])
imagesc([reshape(cx-cxy0,nLenslet,nLenslet),reshape(cy-cxy0,nLenslet,nLenslet)])
axis equal tight
colorbar
subplot(7,3,12)
imagesc(reshape(flux/nPxLenslet^2,nLenslet,nLenslet))
axis equal tight
colorbar
drawnow

%% slopes-to-slopes covariance matrix
nF = nLenslet*2*10;%128;
[fx,fy] = freqspace(nF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*fx;
fy = lf*fy;

r0 = 15e-2;
L0 = 30;
delta = 2*lf/nF;
spectrum = @(u,v) lambda.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
        delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
        (tools.sinc(d*fx).*tools.sinc(d*fy)).^2;
% spectrum = ...
%     fx.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
%     tools.sombrero(1,pi*d*hypot(fx,fy)).^2;
spectrum0 = ...
    phaseStats.spectrum(hypot(fx,fy),atm);

nm = ones(1,2)*nLenslet;
b0 = nF/2+1;
b  = ((1-nLenslet)*sf:sf:sf*(nLenslet-1)) + b0;

tic
covxx = real( fftshift( fft2( fftshift( spectrum([1,0],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covxx(b,b) );
CTBT{1,1} = T;
covyy = real( fftshift( fft2( fftshift( spectrum([0,1],[0,1]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covyy(b,b) );
CTBT{2,2} = T;
cov = real( fftshift( fft2( fftshift( spectrum([0,1],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, cov(b,b) );
CTBT{1,2} = T;
CTBT{2,1} = T';
elapsedTime = toc;
fprintf(' ==> slopes-to-slopes covariance matrix computed in %5.2fs\n',elapsedTime);

%% phase-to-slopes covariance matrix
alpha = 4;
nP = alpha*nLenslet+1;
nPF = nP*2*4;%32;
[fx,fy] = freqspace(nPF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*alpha*fx;
fy = lf*alpha*fy;
delta = 2*lf*alpha/nPF;
spectrum1 = @(u) -lambda.*1i*(fx.*u(1) + fy.*u(2)).*...
    delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
    tools.sinc(d*fx).*tools.sinc(d*fy);

nm = ones(1,2)*nP;
b0 = nPF/2+1;
b  = ((1-nP)*sf:sf:sf*(nP-1)) + b0;

tic
covx  = fftshift(real( fft2( fftshift( spectrum1([1,0]) ) ) ) );
covy  = fftshift(real( fft2( fftshift( spectrum1([0,1]) ) ) ) );
STx = toeplitzBlockToeplitz( nm, nm, covx(b,b) );
STy = toeplitzBlockToeplitz( nm, nm, covy(b,b) );
elapsedTime = toc;
fprintf(' ==> phase-to-slopes covariance matrix computed in %5.2fs\n',elapsedTime);

%%
w = (alpha-1):alpha:nP;
ww = w'*ones(1,nLenslet);
idx = sub2ind( ones(1,2)*nP , ww ,  ww');
mask = tools.piston(nP,'type','logical');
mask_c = tools.piston(nP-4,nP,'type','logical');
mask_c = mask_c(idx);
mask_c_c = repmat( mask_c(:), 2 ,1);
[ix,iy] = meshgrid(0:nP-1);
figure(21)
plot(ix(mask),iy(mask),'.')
ix_c = ix(idx);
iy_c = iy(idx);
line(ix_c(mask_c),iy_c(mask_c),'LineStyle','none','marker','o','color','r')
axis equal tight
xytick = 0:alpha:nP;
set(gca,'xtick',xytick,'ytick',xytick)
grid
legend('Phase','Slopes','location','EastOutside')

%%
[gphs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
cx = cx - cxy0;
cy = cy - cxy0;

ui = linspace(1,nxy,nP);
[xi,yi] = meshgrid( ui );
phs = interp2(gather(gphs),xi,yi);
phs_zm = mask.*phase2nm.*( phs-mean(phs(mask)) );

slopes2Angle = (lambda/2/d);
c = slopes2Angle*[cx.*mask_c(:);cy.*mask_c(:)];
fun = @(x) mtimes4squareBlocks(CTBT,x);

cpx = zeros(nP^2,1);
cpy = zeros(nP^2,1);

tic
[yy,flag,relres,iter,resvec] = my_minres(fun,gather(c),1e-3,50,[],[],[],mask_c_c);
cpx(idx) = yy(1:end/2);
cpy(idx) = yy(1+end/2:end);
phse_2 = STx*cpx + STy*cpy;
phse_2_zm = mask.*phase2nm.*( reshape(phse_2-mean(phse_2(mask)),nP,nP) );
elapsedTime = toc;
fprintf(' ==> phase estimate computed in %5.2fms\n',elapsedTime*1e3);

phse_2_err = phs_zm - phse_2_zm;

figure(23)
subplot(2,3,[1,4])
imagesc([ phs_zm; phse_2_zm])
title(sprintf('Orig.(WF rms [nm] : %5.2f) / Est.Theo.Iter',std(phs_zm(:)) ) )
axis equal tight
colorbar('location','south')
subplot(2,3,[2,6])
imagesc(phse_2_err)
title(sprintf('Est.Theo.Iter wfe rms [nm] : %5.2f',std(phse_2_err(:)) ) )
axis equal tight
colorbar

%%
%{
nIt = 50;
[gphs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
cx = cx - cxy0;
cy = cy - cxy0;
c = slopes2Angle*[cx.*mask_c(:);cy.*mask_c(:)];

for kiT=1:nIt
    
    [yy,flag,relres,iter,resvec] = my_minres(fun,gather(c),1e-3,50,[],[],[],mask_c_c);
    cpx(idx) = yy(1:end/2);
    cpy(idx) = yy(1+end/2:end);
    phse_2 = STx*cpx + STy*cpy;
    phse_2_zm = mask.*phase2nm.*( reshape(phse_2-mean(phse_2(mask)),nP,nP) );

    
end
%}