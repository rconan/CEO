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
wavenumber = 2*pi/atm.wavelength;

nxy = 512;
ceodir = '~/CEO';
cd([ceodir,'/atmosphere'])
clear ceo_atmosphere
unix('make atmosphere.mex')
mex -largeArrayDims -I../include -L../lib -lceo -lcurl -ljsmn -o ceo_atmosphere atmosphere.mex.cu

u = single( L0*gpuArray.linspace(-1,1,nxy) );
[x,y] = meshgrid( u );
phs = ceo_atmosphere(x,y,0,L0,0);
figure(1)
imagesc(u,u,phs)
axis square
colorbar

%% Variance test
fprintf('__ Variance Test __\n')
clear x y ceo_atmosphere
tic
nxy = 2500;
x   = gpuArray.rand(1,nxy,'single');
y   = gpuArray.rand(1,nxy,'single');
L = 100;
x = (2*x-1)*L/2;
y = (2*y-1)*L/2;
phs_var =  gpuArray.zeros(1,nxy,'single');
h = waitbar(0,'Variance Test');
for kxy = 1:nxy
    phs_var(kxy) =  ceo_atmosphere(x(kxy),y(kxy),1,L0,0);
    waitbar(kxy/nxy,h)
end
close(h)
var_num = var(phs_var)*wavenumber^2;
fprintf(' . Theoretical variance: %8.2frd^2\n',phaseStats.variance(atm))
fprintf(' . Numerical variance:   %8.2frd^2\n',var_num)
fprintf(' . Variance ratio: %6.5f\n',var_num/phaseStats.variance(atm))
toc
%% Structure function test I
fprintf('__ Structure Function Test I __\n')
n_sample = 1000;
rho = 0:0.25:5;
rho(1) = 0.1;
nRho = length(rho);
mean_sf = zeros(1,nRho);
std_sf = zeros(1,nRho);
n_plps = 1000;
d_phs = gpuArray.zeros(n_plps,n_sample,'single');
hwb = waitbar(0,'Computing SF ...');
for kRho=1:nRho
    
    phi = gpuArray.rand(1,n_sample,'single')*2*pi;
    zRho = rho(kRho).*exp(1i*phi);
    zxy = (gpuArray.rand(1,n_sample,'single')*2-1)*0.5*L0 + ...
        1i*(gpuArray.rand(1,n_sample,'single')*2-1)*0.5*L0;
    zxy_rho = zxy + zRho;
    tic
    for k_plps = 1:n_plps
        
        phs_xy = ceo_atmosphere(real(zxy),imag(zxy),1,L0,0);
        phs_xy_rho =  ceo_atmosphere(real(zxy_rho),imag(zxy_rho),0,L0,0);
        d_phs(k_plps,:) = phs_xy - phs_xy_rho;
        
    end
    toc
    
    sf = var(d_phs);
    mean_sf(kRho) =  gather( mean(sf) );
    std_sf(kRho)  = gather( std(sf) );
    
    waitbar(kRho/nRho)
    
end
close(hwb)

figure(25)
heb = errorbar(rho,mean_sf*wavenumber^2, std_sf*wavenumber^2);
set(heb','Marker','o','MarkerSize',8,...
    'MarkerFaceColor','r','MarkerEdgeColor','k',...
    'Linewidth',2,'LineStyle','none')
hold all
plot(rho,phaseStats.structureFunction(rho,atm),'Linewidth',2)
hold off
grid
xlabel('Separation [m]')
ylabel('Structure function [rd^2]')

%% Structure function test II
fprintf('__ Structure Function Test II __\n')
L0_ = [1 5 25 300]; 
nL0 = length(L0_);

n_plps = 1000;

nxy = n_sample;
phs_xy = gpuArray.zeros(1,nxy,'single');
phs_xy_rho = gpuArray.zeros(1,nxy,'single');

d_phs = gpuArray.zeros(n_plps,n_sample,'single');

rho = logspace(-2,2,10)';
nRho = length(rho);
mean_sf = zeros(nRho,nL0);
std_sf = zeros(nRho,nL0);
th_sf = zeros(nRho,nL0);

for kL0 = 1:nL0
    
    L0 = L0_(kL0);
    atm.L0 = L0;
    clear ceo_atmosphere

    hwb = waitbar(0,sprintf('Computing SF for L0=%3.0fm ...',L0));
    for kRho=1:nRho
        
        phi = gpuArray.rand(1,n_sample,'single')*2*pi;
        zRho = rho(kRho).*exp(1i*phi);
        zxy = (gpuArray.rand(1,n_sample,'single')*2-1)*0.5*L0 + ...
            1i*(gpuArray.rand(1,n_sample,'single')*2-1)*0.5*L0;
        zxy_rho = zxy + zRho;
        
        tic
        for k_plps = 1:n_plps
            
            phs_xy = ceo_atmosphere(real(zxy),imag(zxy),1,L0,0);
            phs_xy_rho =  ceo_atmosphere(real(zxy_rho),imag(zxy_rho),0,L0,0);
            d_phs(k_plps,:) = phs_xy - phs_xy_rho;
            
        end
        toc
        
        sf = var(d_phs);
        mean_sf(kRho,kL0) =  gather( mean(sf) );
        std_sf(kRho,kL0)  = gather( std(sf) );
        
        waitbar(kRho/nRho)
        
    end
    close(hwb)
    
    th_sf(:,kL0) = phaseStats.structureFunction(rho,atm);
    
    figure(26)
    heb = errorbar(repmat(rho,1,kL0),...
        mean_sf(:,1:kL0)*wavenumber^2,...
        std_sf(:,1:kL0)*wavenumber^2);
    set(heb,'Marker','o','MarkerSize',8,...
        'MarkerFaceColor','r','MarkerEdgeColor','k',...
        'Linewidth',2,'LineStyle','none','color','b')
    hold all
    plot(rho,th_sf(:,1:kL0),'color','k','Linewidth',2)
    hold off
    grid
    xlabel('Separation [m]')
    ylabel('Structure function [rd^2]')
    set(gca,'xscale','log','yscale','log')  
    drawnow
end
for kL0=1:nL0
    text(rho(end),mean_sf(end,kL0)*.7,sprintf('L0=%3.0fm',L0_(kL0)),...
        'VerticalAlignment','top','BackgroundColor','w')
end

%% Zernike test
fprintf('__ Zernike Test __\n')

L0 = 30;
atm.L0 = L0;
nxy = 128;
clear ceo_atmosphere
cd([ceodir,'/atmosphere'])
unix(['sed -i ',...
    '-e ''s/#define _N_LAYER_ [0-9]*/#define _N_LAYER_ ',num2str(atm.nLayer),'/g'' ',...
    '-e ''s/#define _N_PIXEL_ [0-9]*/#define _N_PIXEL_ ',...
    num2str(nxy^2),'/g'' definitions.h']);
unix('cat definitions.h');
cd(ceodir)
unix('make clean lib atmosphere.mex')
cd([ceodir,'/atmosphere'])
mex -largeArrayDims -I../include -L../lib -lceo -o ceo_atmosphere atmosphere.mex.cu

u = 12.5*single(gpuArray.linspace(-1,1,nxy));
[x,y] = meshgrid(u);

ngs = source('zenith',0,'azimuth',0,'height',90e3);
zern = zernike(1:66,25,'resolution',nxy);

nIt = 4000;
zernCoefs = gpuArray.zeros(zern.nMode,nIt,'single');
h = waitbar(0,'Zernike Test !');
for kTau=1:nIt
    phs = ceo_atmosphere(x,y,1,L0,0);
    zern = zern.\phs;
    zernCoefs(:,kTau) = zern.c;
    waitbar(kTau/nIt,h)
end
close(h)
figure(29)
h = semilogy(zern.j,var(zernCoefs,0,2),'ko',...
    zern.j,zernikeStats.variance(zern,atm,ngs),'.-',...
    zern.j,zernikeStats.variance(zern,atm),'ko--');
set(h(1),'MarkerFaceColor','r')
grid
xlabel('Zernike mode')
ylabel('Zernike coef. variance [rd^2]')

%for kTau=1:nTau;phs = ceo_atmosphere(x,y,0,L0,(kTau-1)*tau);imagesc(phs);axis square;colorbar;drawnow;end

%% Taylor (frozen flow) hypothesis test
fprintf('__ Taylor (frozen flow) Hypothesis Test __\n')

tic
phs = ceo_atmosphere(x,y,0,L0,0);
toc

figure(27)
imagesc(u,u,phs)
axis square
colorbar

nIt = 1000;
tau = 1/10;
duration = 5;
nTau = duration/tau;
wind = 10;%.*exp(1i*pi/3);
% wind = 10.*exp(1i*sin(2*pi*(0:nIt-1)*tau*1));

zern = zernike(1:22,25,'resolution',nxy);
zernCoefs = gpuArray.zeros(zern.nMode,nTau,nIt,'single');

h = waitbar(0,'Taylor (frozen flow) hypothesis test!');
for kIt=1:nIt
    
    phs = ceo_atmosphere(x,y,1,L0,0);
    
    for kTau=1:nTau
        phs = ceo_atmosphere(x,y,0,L0,(kTau-1)*tau);
        zern = zern.\phs;
        zernCoefs(:,kTau,kIt) = zern.c;
        %     set(h,'Cdata',h_phs)
        %     drawnow
    end
    
    waitbar(kIt/nIt,h)
    
end
close(h)

tau_ = (0:nTau-1)*tau;
ngs = source;
zcov = zeros(zern.nMode,zern.nMode,nTau);
if matlabpool('size')==0
    matlabpool open
end
tic
parfor kTau=1:nTau
    zcov(:,:,kTau) = ...
        zernikeStats.temporalAngularCovariance(zern,atm,tau_(kTau),ngs,ngs);
end
toc
zcov_diag =cell2mat( ...
    arrayfun( @(x) squeeze( zcov(x,x,:) ) , 1:22, 'UniformOutput', false) );
figure(30)
h_th = plot(tau_,zcov_diag(:,2:8),'LineWidth',2);
grid
xlabel('Time [s]')
ylabel('Zernike coef. covariance [rd^2]')
legend(num2str((2:8)'),0)
hold off

C = mean( bsxfun( @times , zernCoefs(:,1,:) , zernCoefs ) , 3);
hold all
h_num = plot(tau_,C(2:8,:)','.','MarkerSize',15);
hold off

