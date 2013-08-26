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
alpha = 4;
nP = alpha*nLenslet+1;
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
%     [phsi,~,~,~,~] = ceo_imaging(xi,yi,0,L0,0);
    phsi = interp2(gather(phs),xi,yi);
    Gpc = Gpc + reshape(phsi,1,[])'*c';
    waitbar(kSammple/nSample,h)
end
close(h)
slopes2Angle = (lambda/2/d);
Gcc = slopes2Angle^2*Gcc./nSample;
Gpc = slopes2Angle*Gpc./nSample;
figure
ha1 = subplot(1,2,1);
imagesc(Gcc)
axis equal tight
colorbar('location','NorthOutside')
ha2 = subplot(1,2,2);
imagesc(Gpc)
colorbar('location','NorthOutside')

%%
d = D/nLenslet;
nF = nLenslet*2*128;
[fx,fy] = freqspace(nF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*fx;
fy = lf*fy;

r0 = 15e-2;
L0 = 30;
atm = atmosphere(photometry.V,r0,L0);
lambda = atm.wavelength;
delta = 2*lf/nF;
spectrum = @(u,v) lambda.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
        delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
        (tools.sinc(d*fx).*tools.sinc(d*fy)).^2;
% spectrum = ...
%     fx.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
%     tools.sombrero(1,pi*d*hypot(fx,fy)).^2;
spectrum0 = ...
    phaseStats.spectrum(hypot(fx,fy),atm);

%%
nm = ones(1,2)*nLenslet;
b0 = nF/2+1;
b  = ((1-nLenslet)*sf:sf:sf*(nLenslet-1)) + b0;
covxx = real( fftshift( fft2( fftshift( spectrum([1,0],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covxx(b,b) );
CC{1,1} = full(T);
covyy = real( fftshift( fft2( fftshift( spectrum([0,1],[0,1]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covyy(b,b) );
CC{2,2} = full(T);
cov = real( fftshift( fft2( fftshift( spectrum([0,1],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, cov(b,b) );
CC{1,2} = full(T);
CC{2,1} = CC{1,2}';
C = cell2mat(CC);
figure(103)
imagesc([ C;Gcc])
axis equal tight
colorbar('location','NorthOutside')

% figure
% imagesc([covxx,covyy])
%%
Gcc_cell = mat2cell( gather(Gcc) , ones(1,nLenslet*2)*nLenslet, ones(1,nLenslet*2)*nLenslet );
Gcxcx_cell = gather( Gcc_cell(1:nLenslet,1:nLenslet) );
idx = [nLenslet:-1:1 nLenslet*(1:3)+1];
cov_cxcx = cellfun( @(x) x(idx)' , Gcxcx_cell(idx) , 'UniformOutput', false);
Gcycy_cell = gather( Gcc_cell((1:nLenslet)+nLenslet,(1:nLenslet)+nLenslet) );
idx = [nLenslet:-1:1 nLenslet*(1:3)+1];
cov_cycy = cellfun( @(x) x(idx)' , Gcycy_cell(idx) , 'UniformOutput', false);
figure
imagesc([cell2mat(cov_cxcx),cell2mat(cov_cycy)])
%%
nPF = nP*2*32;
[fx,fy] = freqspace(nPF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*alpha*fx;
fy = lf*alpha*fy;
delta = 2*lf*alpha/nPF;
spectrum1 = @(u) -lambda.*1i*(fx.*u(1) + fy.*u(2)).*...
        delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
        tools.sinc(d*fx).*tools.sinc(d*fy);
covx  = fftshift(real( fft2( fftshift( spectrum1([1,0]) ) ) ) );
covy  = fftshift(real( fft2( fftshift( spectrum1([0,1]) ) ) ) );
% c = reshape( cov(1:nLenslet,1:nLenslet) , 1 , [] );
% S = myToeplitz( c,c );
% figure
% imagesc(imag(S))

%%
nm = ones(1,2)*nP;
b0 = nPF/2+1;
b  = ((1-nP)*sf:sf:sf*(nP-1)) + b0;
T = toeplitzBlockToeplitz( nm, nm, covx(b,b) );
Sx = full(T);
T = toeplitzBlockToeplitz( nm, nm, covy(b,b) );
Sy = full(T);
% for k=0:nP-1
%     S = [S ; reshape( cov((nP+1:end)-k,(nP+1:end)-k) , 1 ,[] )];
% end
w = (alpha-1):alpha:nP;
nm = ones(1,nP)*nP;
Sx = mat2cell( Sx , nm, nm );
Sx = cellfun( @(x) x(:,w) , Sx(:,w), 'UniformOutput', false );
Sy = mat2cell( Sy , nm, nm );
Sy = cellfun( @(x) x(:,w) , Sy(:,w), 'UniformOutput', false );
S = [ cell2mat(Sx) cell2mat(Sy) ];
figure(102)
imagesc([S,Gpc])
colorbar('location','NorthOutside')

%%
iC = calibrationVault(C);
iC.cond = 100;
% iC.nThresholded = 1;
%%
[gphs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
    phs = interp2(gather(gphs),xi,yi);
cx = cx - 7.5;
cy = cy - 7.5;
c = slopes2Angle*[cx;cy];
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

%%
scale = 1.2;
phse_1 = S*iC.M*c;
phse_1_zm = scale*phase2nm*( reshape(phse_1-mean(phse_1),nP,nP) );
figure(7)
imagesc(phse_1_zm)
axis square
colorbar

[yy,flag,relres] = lsqr(C,c,1e-2);
phse_2 = S*yy;
phse_2_zm = scale*phase2nm*( reshape(phse_2-mean(phse_2),nP,nP) );
figure(8)
imagesc(phse_2_zm)
axis square
colorbar

phse_1_2_err = phse_2_zm - phse_1_zm;
std(phse_1_2_err(:))
phse_1_err = phse_zm - phse_1_zm;
std(phse_1_err(:))
phse_2_err = phse_zm - phse_2_zm;
std(phse_2_err(:))

