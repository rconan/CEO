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
r0 = atm.r0;
L0 = atm.L0;

 nLenslet = nLenslet_(kRun);
d = D/nLenslet;
nPxLenslet = 16;
cxy0 = 0.5*(nPxLenslet-1);
nxy = nLenslet*nPxLenslet;

r0 = d;
% L0 = 30;
% atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
phase2nm = 1e9;%*lambda/2/pi;

nIt = 200;
ne = nLenslet+1;

compile = true;
figure(102)
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
% unix('make imaging.mex')
% clear ceo_imaging
% mex -largeArrayDims -I../include -L../lib -lceo -o ceo_imaging imaging.mex.cu
%%
cd([ceodir,'/iterativeSolvers'])
unix('make iterativeSolvers.bin')
 fprintf(' ==>>> CG (N=%d)\n',nLenslet)
 tic
unix(sprintf('./a.out %3.1f CG > CVGCE_CG_%03d_%03d.log',D,nIt,nLenslet));
toc
fprintf(' ==>>> MINRES (N=%d)\n',nLenslet)
tic
unix(sprintf('./a.out %3.1f MINRES > CVGCE_MINRES_%03d_%03d.log',D,nIt,nLenslet));
toc

c = loadBin('centroids',[nLenslet,nLenslet*2]);
subplot(3,4,[9,12])
imagesc(c)
axis equal tight
colorbar('location','north')

drawnow

end
%%
pup = tools.piston(ne,'type','logical','shape','square');
% ps = loadBin('phaseScreen',[nxy,nxy]);
ps = phase2nm*loadBin(sprintf('CVGCE_phaseScreenLowRes_%03d',nLenslet),[ne,ne]);
ps = pup.*(ps - mean(ps(pup)));
% subplot(2,3,[1,2])
% imagesc(ps)
% axis square
% xlabel(colorbar('location','northOutside'),'[nm]')

%ps_e = phase2nm*loadBin('phaseScreenEst',[ne,ne]);
%ps_e = phase2nm*loadBin(sprintf('CG_phaseEst_%3d_%2d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = phase2nm*loadBin(sprintf('CVGCE_MINRES_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = bsxfun( @minus, ps_e, mean(ps_e(pup,:),1) );
ps_e = bsxfun( @times, pup(:), ps_e);
ps_e_k = reshape(ps_e(:,nIt),[ne,ne]);
subplot(3,4,[1,6])
imagesc([ps,ps_e_k])
axis equal tight
xlabel(colorbar('location','northOutside'),'[nm]')

ps_err = bsxfun( @minus, ps(:), ps_e);
rms_ps_err = std(ps_err(pup(:),:));
ps_err_k = reshape(ps_err(:,nIt),[ne,ne]);
subplot(3,4,[3,8])
imagesc(ps_err_k)
axis equal tight
title(sprintf('wfe=%6.2fnm',rms_ps_err(nIt)))
xlabel(colorbar('location','southOutside'),'[nm]')

%%
ps = phase2nm*loadBin(sprintf('CVGCE_phaseScreenLowRes_%03d',nLenslet),[ne,ne]);
ps = pup.*(ps - mean(ps(pup)));
ps_e = phase2nm*loadBin(sprintf('CVGCE_MINRES_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = bsxfun( @minus, ps_e, mean(ps_e(pup,:),1) );
ps_e = bsxfun( @times, pup(:), ps_e);
ps_err = bsxfun( @minus, ps(:), ps_e);
rms_ps_err_minres = std(ps_err(pup(:),:));
ps_e = phase2nm*loadBin(sprintf('CVGCE_CG_phaseEst_%03d_%03d',nIt,nLenslet),[ne*ne,nIt]);
ps_e = bsxfun( @minus, ps_e, mean(ps_e(pup,:),1) );
ps_e = bsxfun( @times, pup(:), ps_e);
ps_err = bsxfun( @minus, ps(:), ps_e);
rms_ps_err_cg = std(ps_err(pup(:),:));

figure(314)   
subplot(3,2,kRun)
loglog(1:nIt,rms_ps_err_minres,'.-',1:nIt,rms_ps_err_cg,'.-')
grid
title(sprintf('N=%d',nLenslet))
xlabel('Iteration #')
ylabel('WFE [nm]')
legend('MINRES','CG',0)

drawnow

end

%% MINRES CONVERGENCE
%{
figure(101)
for kRun = 1:length(D_)
    
D = D_(kRun);
nLenslet = nLenslet_(kRun);

logdata    = fileread(sprintf('CVGCE_CG_%03d_%03d.log',nIt,nLenslet));
logexcerpt = regexp(logdata,'[^\n]*Solver residue norm[^\n]*','match');
res_cg     = cellfun( @(x) str2double(x(33:end)) , logexcerpt);

logdata    = fileread(sprintf('CVGCE_MINRES_%03d_%03d.log',nIt,nLenslet));
logexcerpt = regexp(logdata,'[^\n]*Solver residue norm[^\n]*','match');
res_minres = cellfun( @(x) str2double(x(33:end)) , logexcerpt);
subplot(2,3,kRun)
loglog(1:nIt,res_minres,'.-',1:nIt,res_cg,'.-')
grid
title(sprintf('N=%d',nLenslet))
xlabel('Iteration #')
ylabel('WFE [nm]')
legend('MINRES','CG',0)

end
%}

%% 
%{
nu = 2*nLenslet-1;
AA = loadBin('aaCovariance',[nu,nu*4]);
figure(301),subplot(2,1,1),imagesc(AA),axis equal tight

nm = ones(1,2)*nLenslet;
u = 1:nu;
T = toeplitzBlockToeplitz( nm, nm, AA(:,u) );
CTBT{1,1} = T;
u = u + nu;
T = toeplitzBlockToeplitz( nm, nm, AA(:,u) );
CTBT{1,2} = T;
u = u + nu;
T = toeplitzBlockToeplitz( nm, nm, AA(:,u) );
CTBT{2,1} = T;
u = u + nu;
T = toeplitzBlockToeplitz( nm, nm, AA(:,u) );
CTBT{2,2} = T;

AA_full = cell2mat(cellfun( @(x)full(x) , CTBT , 'UniformOutput', false));
% figure(302),imagesc(AA_full)
% M = spdiags(diag(AA_full), 0 , nLenslet^2*2, nLenslet^2*2);

nu = 2*nLenslet;
PA = loadBin('paCovariance',[nu,nu*2]);
figure(301),subplot(2,1,2),imagesc(PA),axis equal tight

nm = [ne,nLenslet];
u = 1:nu;
ST{1} = toeplitzBlockToeplitz( nm, nm, PA(:,u) );
u = u + nu;
ST{2} = toeplitzBlockToeplitz( nm, nm, PA(:,u) );

figure(303),imagesc(cell2mat(cellfun( @(x)full(x) , ST , 'UniformOutput', false)))

c = loadBin('centroids');
ps = phase2nm*loadBin(sprintf('CVGCE_phaseScreenLowRes_%03d',nLenslet));
ps = ps - mean(ps);

fun = @(x) mtimes4squareBlocks(CTBT,x);
% [yy,flag,relres,iter,resvec] = minres(fun,c,1e-3,50);

yy = minres(fun,c,5e-2,200);
ps_e = -phase2nm*(ST{1}*yy(1:end/2) + ST{2}*yy(1+end/2:end));
fprintf(' ==>> WFE=%5.2fnm\n',std(ps-ps_e))

% M = diag(diag(AA_full));
L = chol(AA_full,'lower');
yy = minres(fun,c,5e-2,200,L,L');
ps_e = -phase2nm*(ST{1}*yy(1:end/2) + ST{2}*yy(1+end/2:end));
fprintf(' ==>> WFE=%5.2fnm\n',std(ps-ps_e))

A = cellfun( @(x)full(x) , CTBT , 'UniformOutput', false);
A{1,2} = zeros(nLenslet^2,'single');
A{2,1} = zeros(nLenslet^2,'single');
A = cell2mat(A);
iL = ichol(sparse(double((A))));
iL = full(iL);
yy = minres(fun,c,5e-2,200,iL,iL');
ps_e = -phase2nm*(ST{1}*yy(1:end/2) + ST{2}*yy(1+end/2:end));
fprintf(' ==>> WFE=%5.2fnm\n',std(ps-ps_e))

% ps = reshape( ps, ne, ne);
% ps_e = reshape( ps_e, ne, ne);
wfs = shackHartmann(nLenslet,nLenslet*16);
G = wfs.sparseGradientMatrix;
D = delsq( numgrid('S', ne ) );
%}