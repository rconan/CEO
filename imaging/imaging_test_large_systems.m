%% 
% gpuDevice([])

 r0 = 15e-2;
 L0 = 30;
 atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
D = 0.4;
phase2nm = 1e9*lambda/2/pi;

%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

nLenslet = 4;
d = D/nLenslet;
nPxLenslet = 16;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;

%%
clear ceo_imaging
ceodir = '~/CEO';
cd([ceodir,'/include'])
unix(['sed -i ',...
    '-e ''s/#define N_SIDE_LENSLET [0-9]*/#define N_SIDE_LENSLET ',num2str(nLenslet),'/g'' ',...
    '-e ''s/#define _N_PX_PUPIL_ [0-9]*/#define _N_PX_PUPIL_ ',num2str(nPxLenslet),'/g'' ',...
    '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
    num2str((nLenslet*nPxLenslet)^2),'/g'' definitions.h']);
unix('cat definitions.h');
cd(ceodir)
unix('make clean all')
cd([ceodir,'/imaging'])
unix('make imaging.mex')
clear ceo_imaging
mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu

%%

u = single( 0.5*D*gpuArray.linspace(0,2,nxy) );
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
nF = 16;%1024;%2^nextpow2(nLenslet*10);%nLenslet*2*10;%128;
[fx,fy] = freqspace(nF,'meshgrid');
sf = 1;%4;
lf = sf/(d*2);
fx = lf*fx;
fy = lf*fy;

r0 = 15e-2;
L0 = 30;
delta = 2*lf/nF;
spectrum = @(fx,fy,u,v) lambda.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
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
covxx = real( fftshift( fft2( fftshift( spectrum(fx,fy,[1,0],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covxx(b,b) );
CTBT{1,1} = T;
covyy = real( fftshift( fft2( fftshift( spectrum(fx,fy,[0,1],[0,1]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, covyy(b,b) );
CTBT{2,2} = T;
cov = real( fftshift( fft2( fftshift( spectrum(fx,fy,[0,1],[1,0]) ) ) ) );
T = toeplitzBlockToeplitz( nm, nm, cov(b,b) );
CTBT{1,2} = T;
CTBT{2,1} = T';
elapsedTime = toc;
fprintf(' ==> slopes-to-slopes covariance matrix computed in %5.2fs\n',elapsedTime);

%% phase-to-slopes covariance matrix
alpha = 2;
nP = alpha*nLenslet+1;
nPF = 2^nextpow2(nP*8);%nP*2*4;%32;
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
w = (alpha/2+1):alpha:nP;
ww = w'*ones(1,nLenslet);
idx = sub2ind( ones(1,2)*nP , ww ,  ww');
mask = tools.piston(nP,'type','logical');
mask_c = tools.piston(nP-4,nP,'type','logical');
mask_c = mask_c(idx);
mask_c_c = repmat( mask_c(:), 2 ,1);
[ix,iy] = meshgrid(0:nP-1);
if nLenslet<41
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
end
%% Deformable Mirror
bifLR = influenceFunction('monotonic',0.5);
dmLR = deformableMirror(nLenslet+1,'modes',bifLR,'resolution',nP);
F = bifLR.modes;
bif = influenceFunction('monotonic',0.5);
dm = deformableMirror(nLenslet+1,'modes',bif,'resolution',nxy);
pupil = tools.piston(nxy-nPxLenslet*0,nxy,'type','logical');
%%
[gphs,frame,cx,cy,flux] = ceo_imaging(x,y,1,L0,0);
cx = cx - cxy0;
cy = cy - cxy0;

phs = gather(gphs);
phs_zm = pupil.*phase2nm.*( phs-mean(phs(pupil)) );

slopes2Angle = (lambda/2/d);
c = slopes2Angle*[cx.*mask_c(:);cy.*mask_c(:)];
fun = @(x) mtimes4squareBlocks(CTBT,x,mask_c(:));

cpx = zeros(nP^2,1);
cpy = zeros(nP^2,1);

tic
% [yy,flag,relres,iter,resvec] = my_minres(fun,gather(c),1e-3,50,[],[],[],mask_c_c);
%[yy,flag,relres,iter,resvec] = minres(fun,gather(c),1e-3,50);
[yy,flag,relres,iter,resvec] = pcg(fun,gather(c),1e-3,50);
cpx(idx) = yy(1:end/2);
cpy(idx) = yy(1+end/2:end);
phse_2 = STx*cpx + STy*cpy;
phse_2_zm = mask.*phase2nm.*( reshape(phse_2-mean(phse_2(mask)),nP,nP) );
dm.coefs = lsqr(F,phse_2_zm(:),1e-3,50);
phs_dm = dm.surface;
phs_dm_zm = pupil.*( phs_dm-mean(phs_dm(pupil)) );
elapsedTime = toc;
fprintf(' ==> phase estimate computed in %5.2fms\n',elapsedTime*1e3);

% phse_2_err = phs_zm - phse_2_zm;
wfe = phs_zm - phs_dm_zm;
wfe_rms0 = std(wfe(pupil));
marechal_strehl = exp(-(1e-9*wfe_rms0*2*pi/2.2e-6).^2);
fprintf(' ==> Marechal Strehl: %5.2f%%\n',marechal_strehl*1e2);

figure(23)
subplot(2,3,[1,4])
imagesc([ phs_zm; phs_dm_zm])
title(sprintf('Orig.(WF rms [nm] : %5.2f) / Est.Theo.Iter',std(phs_zm(pupil)) ) )
axis equal tight
colorbar('location','south')
subplot(2,3,[2,6])
imagesc(wfe)
title(sprintf('Est.Theo.Iter. wfe rms [nm] : %5.2f', wfe_rms0) )
axis equal tight
colorbar
drawnow

%% OPEN-LOOP CONTROL
nIt = 500;

c_dm = zeros(nLenslet^2*2,1);
c_est = zeros(nLenslet^2*2,1);
phs_dm_zm = zeros(nxy);
gain = 1;
tau = 1e-3;
elt_in_loop = zeros(1,nIt);
wf_rms = zeros(1,nIt);
wfe_rms = zeros(1,nIt);

fprintf(' ==> %s: Open-loop control started!\n',datestr(now,'dd mmm. yyyy @ HH:MM:SS'));

t_loop = tic;
for kIt=1:nIt
            
    [gphs,frame,cx,cy,flux] = ceo_imaging(x,y,0,L0,kIt*tau);
    cx = cx - cxy0;
    cy = cy - cxy0;
    c = gather( slopes2Angle*[cx.*mask_c(:);cy.*mask_c(:)] );
    
    phs = gather(gphs);
    phs_zm = pupil.*phase2nm.*( phs-mean(phs(pupil)) );
    
    t_in_loop = tic;

    [c_dm,flag_minres,relres_minres,iter_minres,resvec_minres] = minres(fun,c,5e-2,5,[],[],c_dm);    
    %c_dm = minres(fun,c,5e-2,10,[],[],c_dm);    
    
    cpx(idx) = c_dm(1:end/2);
    cpy(idx) = c_dm(1+end/2:end);
    phse_2 = STx*cpx + STy*cpy;
    phse_2_zm = mask.*phase2nm.*( reshape(phse_2-mean(phse_2(mask)),nP,nP) );
    [dm.coefs,flag_lsqr,relres_lsqr,iter_lsqr,resvec_lsqr] = lsqr(F,phse_2_zm(:),5e-2,5,[],[],dm.coefs);
    %dm.coefs = lsqr(F,phse_2_zm(:),5e-2,5,[],[],dm.coefs);
    
    elt_in_loop(kIt) = toc(t_in_loop);
    
    phs_dm = dm.surface;
    phs_dm_zm = pupil.*( phs_dm-mean(phs_dm(pupil)) );
   
    wfe = phs_zm - phs_dm_zm;
    wf_rms(kIt) = std(phs_zm(pupil));
    wfe_rms(kIt) = std(wfe(pupil));
    
%     fprintf(' -- It#:%d - Est.Theo.Iter. wfe rms [nm] : %5.2f [ %5.2f ] \n',...
%         kIt, wfe_rms(kIt),wfe_rms0)
    
%     figure(23)
%     subplot(2,3,[1,4])
%     imagesc([ phs_zm; phs_dm_zm])
%     title(sprintf('Orig.(WF rms [nm] : %5.2f) / Est.Theo.Iter', wf_rms(kIt)) )
%     axis equal tight
%     colorbar('location','south')
%     subplot(2,3,[2,6])
%     imagesc(wfe)
%     title(sprintf('It#:%d - Est.Theo.Iter. wfe rms [nm] : %5.2f', kIt, wfe_rms(kIt)) )
%     axis equal tight
%     colorbar
% 
%     drawnow

end
elt_loop = toc(t_loop);

fprintf(' ==> %s: Open-loop control ended (%.0fs)!\n',datestr(now,'dd mmm. yyyy @ HH:MM:SS'),elt_loop);

fprintf(' ==> wavefront reconstruction computing time: %5.0fms +/- %2.0fms\n',...
    median(elt_in_loop)*1e3,std(elt_in_loop)*1e3);
fprintf(' ==> average number of iterations per second: time: %5.2f\n',nIt/elt_loop);

