%% 
% gpuDevice([])

nLenslet_ = [20 40 64 84 150];
D_        = [3.6 8 5 42 30];
    
nIt = 200;
rms_ps_err_minres = zeros(nIt,length(D_));

for kRun = 1:length(D_)
    
D = D_(kRun);

%  atm = atmosphere(photometry.V,r0,L0,...
%      'altitude',[0, 500, 1000, 2000, 5000, 8000. 13000],...
%      'fractionnalR0',[0.2, 0.1, 0.1, 0.3, 0.2, 0.05, 0.05],...
%      'windSpeed',[10, 5, 7.5, 5, 10, 12, 15],...
%      'windDirection',[0, 0.25, 0.5, 1, 1.5, 1.75, 2]);
% atm = gmtAtmosphere(1);
% r0 = atm.r0;
% L0 = atm.L0;

 nLenslet = nLenslet_(kRun);
d = D/nLenslet;
nPxLenslet = 16;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;

r0 = d;
L0 = 30;
atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
phase2nm = 1e9*lambda/2/pi;

ne = 2*nLenslet+1;

compile = true;

if compile
%%
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
cd([ceodir,'/iterativeSolvers'])
unix('make iterativeSolvers.bin')
% unix('make imaging.mex')
% clear ceo_imaging
% mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu
%%
% fprintf(' ==>>> CG (N=%d)\n',nLenslet)
% tic
% unix(sprintf('./a.out CG > CG_%03d_%03d.log',nIt,nLenslet));
% toc
fprintf(' ==>>> STATS_MINRES (N=%d)\n',nLenslet)
tic
unix(sprintf('./a.out %3.1f MINRES > STATS_MINRES_%03d_%03d.log',D,nIt,nLenslet));
toc
end
%%
ps = phase2nm*loadBin(sprintf('STATS_phaseScreenLowRes_%03d',nLenslet),[ne*ne,nIt]);
ps = bsxfun( @minus, ps, mean(ps,1) );
ps_e = phase2nm*loadBin(sprintf('STATS_MINRES_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = bsxfun( @minus, ps_e, mean(ps_e,1) );
ps_err = ps - ps_e;
rms_ps_err_minres(:,kRun) = std(ps_err);
% ps_e = phase2nm*loadBin(sprintf('CG_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
% ps_e = bsxfun( @minus, ps_e, mean(ps_e,1) );
% ps_err = bsxfun( @minus, ps(:), ps_e);
% rms_ps_err_cg = std(ps_err);

% figure(kRun+5)
% loglog(1:nIt,rms_ps_err_minres,'.-',1:nIt,rms_ps_err_cg,'.-')
% grid
% title(sprintf('N=%d',nLenslet))
% xlabel('Iteration #')
% ylabel('WFE [nm]')
% legend('MINRES','CG',0)
% 
% drawnow

end
%%
figure
boxplot(rms_ps_err_minres,nLenslet_)
grid
xlabel('Lenslet linear size')
ylabel('WFE [nm]')

%%
med_rms_ps_err_minres = median(rms_ps_err_minres)
h = @(fm) phase2nm*sqrt( ...
    phaseStats.variance(atm)-integral(@(x)2.*pi*x.*phaseStats.spectrum(x,atm),0,fm));
d = D_./nLenslet_;
sfit = arrayfun( h, 1./(2*d) )
((med_rms_ps_err_minres./phase2nm)).^2
((sfit./phase2nm)).^2