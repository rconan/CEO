%% 
% gpuDevice([])

 r0 = 15e-2;
 L0 = 30;
 atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
D = 8;
%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

nLenslet = 20;
d = D/nLenslet;
nPxLenslet = 8;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;
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
imagesc([reshape(cx-cxy0,nLenslet,nLenslet),reshape(cy-cxy0,nLenslet,nLenslet)])
axis equal tight
colorbar('location','NorthOutside')
ax4 = axes('pos',[ 0.5979    0.445   0.1879    0.1799]);
imagesc(reshape(flux/256,nLenslet,nLenslet))
axis equal tight
colorbar('location','NorthOutside')

%%
nSample = 500;
Gcc = gpuArray.zeros(nLenslet^2*2,'single');
alpha = 4;
nP = alpha*nLenslet+1;
Gpc = gpuArray.zeros( nP^2 , nLenslet^2*2, 'single' );
h = waitbar(0,'Building covariance matrices ...!');
ui = linspace(1,nxy,nP);
[xi,yi] = meshgrid( ui );
for kSammple=1:nSample
    [phs,~,cx,cy,~] = ceo_imaging(x,y,1,L0,0);
    cx = cx - cxy0;
    cy = cy - cxy0;
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
CTBT{1,1} = T;
CC{1,1} = full(T);
covyy = real( fftshift( fft2( fftshift( spectrum([0,1],[0,1]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covyy(b,b) );
CTBT{2,2} = T;
CC{2,2} = full(T);
cov = real( fftshift( fft2( fftshift( spectrum([0,1],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, cov(b,b) );
CTBT{1,2} = T;
CTBT{2,1} = T';
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
idx = [nLenslet:-1:1 nLenslet*(1:nLenslet-1)+1];
cov_cxcx = cellfun( @(x) x(idx)' , Gcxcx_cell(idx) , 'UniformOutput', false);
Gcycy_cell = gather( Gcc_cell((1:nLenslet)+nLenslet,(1:nLenslet)+nLenslet) );
idx = [nLenslet:-1:1 nLenslet*(1:nLenslet-1)+1];
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
STx = toeplitzBlockToeplitz( nm, nm, covx(b,b) );
Sx = full(STx);
STy = toeplitzBlockToeplitz( nm, nm, covy(b,b) );
Sy = full(STy);
Sfull = [Sx Sy];
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
ww = w'*ones(1,nLenslet);
idx = sub2ind( ones(1,2)*nP , ww ,  ww');
mask = tools.piston(nP,'type','logical');
mask_c = tools.piston(nP-4,nP,'type','logical');
mask_c = mask_c(idx);
[ix,iy] = meshgrid(0:nP-1);
figure(21)
plot(ix(mask),iy(mask),'.')
ix_c = ix(idx);
iy_c = iy(idx);
line(ix_c(mask_c),iy_c(mask_c),'LineStyle','none','marker','o','color','r')
axis equal tight
xytick = [0:alpha:nP];
set(gca,'xtick',xytick,'ytick',xytick)
grid
legend('Phase','Slopes','location','EastOutside')
%%
mask_c_c = repmat( mask_c(:), 2 ,1);
iC = calibrationVault(C(mask_c_c,mask_c_c));
iC.cond = 1000;
% iC.nThresholded = 1;
%%
iGcc = calibrationVault(Gcc(mask_c_c,mask_c_c));
iGcc.cond = 1000;
%%
[gphs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
    phs = interp2(gather(gphs),xi,yi);
cx = cx - cxy0;
cy = cy - cxy0;
c = slopes2Angle*[cx(mask_c);cy(mask_c)];
M = Gpc(mask,mask_c_c)*iGcc.M;
phse = zeros(nP^2,1);
phse(mask) = gather( M*c );
phase2nm = 500/2/pi;
phs_zm = mask.*phase2nm.*( phs-mean(phs(mask)) );
phse_zm = mask.*phase2nm.*( reshape(phse-mean(phse(mask)),nP,nP) );
phs_err = phs_zm - phse_zm;

scale = 1;
cpx = zeros(nP^2,1);
cpy = zeros(nP^2,1);

yy = iC.M*gather(c);
% phse_1 = S*yy;

cpx(idx(mask_c)) = yy(1:end/2);
cpy(idx(mask_c)) = yy(1+end/2:end);
% cp = [cpx ; cpy];
% phse_1 = Sfull*cp;
phse_1 = STx*cpx + STy*cpy;
phse_1_zm = mask.*scale.*phase2nm.*( reshape(phse_1-mean(phse_1(mask)),nP,nP) );

% c = slopes2Angle*[cx;cy];
% [yy0,flag,relres,iter,resvec] = minres(C,gather(c),1e-3,50);
c = slopes2Angle*[cx.*mask_c(:);cy.*mask_c(:)];
fun = @(x) mtimes4squareBlocks(CTBT,x);
tic;
[yy,flag,relres,iter,resvec] = my_minres(fun,gather(c),1e-3,50,[],[],[],mask_c_c);
toc
%yy = yy.*mask_c_c;
cpx(idx) = yy(1:end/2);
cpy(idx) = yy(1+end/2:end);
% phse_2 = S*yy;
phse_2 = STx*cpx + STy*cpy;
phse_2_zm = mask.*scale.*phase2nm.*( reshape(phse_2-mean(phse_2(mask)),nP,nP) );

phse_1_2_err = phse_2_zm - phse_1_zm;
% std(phse_1_2_err(:))
phse_1_err = phs_zm - phse_1_zm;
% std(phse_1_err(:))
phse_2_err = phs_zm - phse_2_zm;
% std(phse_2_err(:))

figure(23)
subplot(2,1,1)
imagesc([ phs_zm, phse_zm, phse_1_zm, phse_2_zm])
title(sprintf('Orig.(WF rms [nm] : %5.2f) / Est.Num.Inv / Est.Theo.Inv / Est.Theo.Iter',std(phs_zm(:)) ) )
axis equal tight
colorbar
subplot(2,1,2)
imagesc([phs_err, phse_1_err, phse_2_err])
title(sprintf(...
    'Est.Num.Inv wfe rms [nm] : %5.2f / Est.Theo.Inv wfe rms [nm] : %5.2f / Est.Theo.Iter wfe rms [nm] : %5.2f',...
    std(phs_err(:)),std(phse_1_err(:)),std(phse_2_err(:)) ) )
axis equal tight
colorbar
