%% 
% gpuDevice([])

nLenslet_ = [20 40 64 84 150];
D_        = [3.6 8 5 42 30];

for kRun = 1:length(D_)
    
D = D_(kRun);

%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
 atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

 nLenslet = nLenslet_(kRun);
d = D/nLenslet;
nPxLenslet = 16;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;

r0 = d;
L0 = 30;
%atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
phase2nm = 1e9*lambda/2/pi;

nIt = 200;
ne = 2*nLenslet+1;

compile = true;

if compile
%%
ceodir = '~/CEO';
cd([ceodir,'/include'])
unix(['sed -i ',...
    '-e ''s/#define N_SIDE_LENSLET [0-9]*/#define N_SIDE_LENSLET ',num2str(nLenslet),'/g'' ',...
    '-e ''s/#define _N_PX_PUPIL_ [0-9]*/#define _N_PX_PUPIL_ ',num2str(nPxLenslet),'/g'' ',...
    '-e ''s/#define _N_LAYER_ [0-9]*/#define _N_LAYER_ ',num2str(atm.nLayer),'/g'' ',...
    '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
    num2str((nLenslet*nPxLenslet)^2),'/g'' definitions.h']);
unix('cat definitions.h');
cd(ceodir)
unix('make clean all')
cd([ceodir,'/iterativeSolvers'])
unix('make iterativeSolvers.bin')
% unix('make imaging.mex')
% clear ceo_imaging
% mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu
%%
% fprintf(' ==>>> CG (N=%d)\n',nLenslet)
% tic
%unix(sprintf('./a.out %3.1f CG > CVGCE_CG_%03d_%03d.log',D,nIt,nLenslet));
%toc
fprintf(' ==>>> MINRES (N=%d)\n',nLenslet)
tic
unix(sprintf('./a.out %3.1f MINRES > CVGCE_MINRES_%03d_%03d.log',D,nIt,nLenslet));
toc

end
%%
pup = tools.piston(ne,'type','logical');
pup = pup(:);
% ps = loadBin('phaseScreen',[nxy,nxy]);
ps = phase2nm*loadBin(sprintf('CVGCE_phaseScreenLowRes_%03d',nLenslet),[ne*ne,nIt]);
ps = bsxfun( @minus, ps, mean(ps(pup,:),1) );
ps = bsxfun( @times, pup, ps);
ps_k = reshape(ps(:,nIt),[ne,ne]);
figure(102)
% subplot(2,3,[1,2])
% imagesc(ps)
% axis square
% xlabel(colorbar('location','northOutside'),'[nm]')

%c = loadBin('centroids',[nLenslet,nLenslet*2]);
%subplot(3,4,[9,12])
%imagesc(c)
%axis equal tight
%colorbar('location','north')

%ps_e = phase2nm*loadBin('phaseScreenEst',[ne,ne]);
%ps_e = phase2nm*loadBin(sprintf('CG_phaseEst_%3d_%2d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = phase2nm*loadBin(sprintf('CVGCE_MINRES_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = bsxfun( @minus, ps_e, mean(ps_e(pup,:),1) );
ps_e = bsxfun( @times, pup, ps_e);
ps_e_k = reshape(ps_e(:,nIt),[ne,ne]);
subplot(3,4,[1,6])
imagesc([ps_k,ps_e_k])
axis equal tight
xlabel(colorbar('location','northOutside'),'[nm]')

ps_err = bsxfun( @minus, ps, ps_e);
rms_ps_err = std(ps_err);
ps_err_k = reshape(ps_err(:,nIt),[ne,ne]);
subplot(3,4,[3,8])
imagesc(ps_err_k)
axis equal tight
title(sprintf('wfe=%6.2fnm',rms_ps_err(nIt)))
xlabel(colorbar('location','southOutside'),'[nm]')

%%
rms_ps = std(ps);
rms_ps_err_minres = std(ps_err);

logdata    = fileread(sprintf('CVGCE_MINRES_%03d_%03d.log',nIt,nLenslet));
logexcerpt = regexp(logdata,'[^\n]*WAVEFRONT ESTIMATION: Elapsed time[^\n]*','match');
est_ET = cellfun( @(x) str2double(x(42:49)) , logexcerpt);

u = 1:nIt;
figure(3141)   
subplot(2,3,kRun)
hold all
loglog(u,rms_ps,'k',u,rms_ps_err_minres,'.-')
grid
title(sprintf('N=%d - recon. time: %4.2fms+/-%2.0f\\mus',...
    nLenslet,mean(est_ET),std(est_ET)*1e3))
xlabel('Iteration #')
ylabel('WFE [nm]')
%legend('MINRES','CG',0)

drawnow

end

%%
%{
logdata    = fileread(sprintf('CVGCE_MINRES_%03d_%03d.log',nIt,nLenslet));
logexcerpt = regexp(logdata,'[^\n]*Solver residue norm[^\n]*','match');
res_minres = cellfun( @(x) str2double(x(33:end)) , logexcerpt);
logexcerpt = regexp(logdata,'[^\n]*r norm=[^\n]*','match');
res_it_minres = cellfun( @(x) str2double(x(8:end)) , logexcerpt);
w = 1:length(res_it_minres);
nStep = length(res_it_minres)/nIt;
figure(1)
hold all
plot(w(nStep:nStep:end),res_minres,'--.',w,res_it_minres,':.')
%}