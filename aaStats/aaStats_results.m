fileID = fopen('aaStats_results.txt');
C = textscan(fileID,'%f %f %f %f','headerLines',1);
fclose(fileID);
C = cell2mat(C);
C(C==0) = NaN;
N = C(:,1);
%%
figure(50)
subplot(1,2,1)
ht = loglog(N,C(:,2:end),'o--');
grid
xlabel('Lenslet Array Size')
ylabel('Runtime [ms]')
legend('MVM Full','MVM Compressed','Compressed cov.',0)
set(gca,'xtick',N)
set(ht,'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2)
arrayfun( @(k) set(ht(k),'MarkerFaceColor',get(ht(k),'color')), 1:length(ht) )
%%
0.5*(N(1:5)/N(3)).^2
0.1*(N(:).^1.5.*log(N(:)))./(N(3).^1.5.*log(N(3)))
%%
b_full = N.^4*4*4/2^20;
b_comp = (2*N-1).^2*4*4*2/2^20;
subplot(1,2,2)
hb = loglog(N,[b_full,b_comp],'o--');
grid
xlabel('Lenslet Array Size')
ylabel('Matrix Memory Rqt. [MB]')
legend('Full','Compressed',0)
set(gca,'xtick',N)
set(hb,'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2)
arrayfun( @(k) set(hb(k),'MarkerFaceColor',get(hb(k),'color')), 1:length(hb) )
%%
nLenslet = 40;
nSlope   = 2*nLenslet^2;
fid = fopen(sprintf('mvm%03d.bin',nLenslet), 'r');
y = fread(fid, nSlope, 'float');
fclose(fid);
% slopes-to-slopes covariance matrix
addpath('../imaging')
r0 = 15e-2;
L0 = 30;
atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
d = 0.1;
nF = 1024;%2^nextpow2(nLenslet*10);%nLenslet*2*10;%128;
[fx,fy] = freqspace(nF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*fx;
fy = lf*fy;

delta = 2*lf/nF;
spectrum = @(fx,fy,u,v) lambda.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
        delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
        (tools.sinc(d*fx).*tools.sinc(d*fy)).^2;

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

y0 = mtimes4squareBlocks(CTBT,(0:nSlope-1)',true(nSlope/2,1));

figure(33)
subplot(2,1,1)
plot([y0,y],'.')
subplot(2,1,2)
plot(y0-y)