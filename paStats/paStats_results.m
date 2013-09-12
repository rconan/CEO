fileID = fopen('paStats_results.txt');
C = textscan(fileID,'%f %f %f %f','headerLines',1);
fclose(fileID);
C = cell2mat(C);
C(C==0) = NaN;
N = C(:,1);
%%
figure(50)
subplot(1,2,1)
if ishandle(50)
    hold on
end
ht = loglog(N,C(:,2:end),'s--');
if ishandle(50)
    hold off
else
    grid
end
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
nP = 2*N + 1;
b_full = nP.^4*4*2/2^20;
b_comp = (2*nP-1).^2*2*4*2/2^20;
subplot(1,2,2)
if ishandle(50)
    hold on
end
hb = loglog(N,[b_full,b_comp],'s--');
if ishandle(50)
    hold off
else
    grid
end
xlabel('Lenslet Array Size')
ylabel('Matrix Memory Rqt. [MB]')
legend('Full','Compressed',0)
set(gca,'xtick',N)
set(hb,'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2)
arrayfun( @(k) set(hb(k),'MarkerFaceColor',get(hb(k),'color')), 1:length(hb) )
%%
nLenslet = 20;
nSlope   = 2*nLenslet^2;
nP = 2*nLenslet + 1;
fid = fopen(sprintf('mvm%03d.bin',nP), 'r');
y = fread(fid, nP^2, 'float');
fclose(fid);
% slopes-to-slopes covariance matrix
addpath('../imaging')
r0 = 15e-2;
L0 = 30;
atm = atmosphere(photometry.V,r0,L0,'windSpeed',10,'windDirection',0);
lambda = atm.wavelength;
d = 0.1;
% phase-to-slopes covariance matrix
alpha = 2;
nP = alpha*nLenslet+1;
nPF = 2048;%2^nextpow2(nP*8);%nP*2*4;%32;
[fx,fy] = freqspace(nPF,'meshgrid');
sf = 4;
lf = sf/(d*2);
fx = lf*alpha*fx;
fy = lf*alpha*fy;
delta = 2*lf*alpha/nPF;
spectrum1 = @(u) -lambda.*1i*(fx.*u(1) + fy.*u(2)).*...
    delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*...
    tools.sinc(d*fx).*tools.sinc(d*fy);
spectrum2 = @(u) -lambda.*1i*(fx.*u(1) + fy.*u(2)).*...
    delta.^2.*phaseStats.spectrum(hypot(fx,fy),atm);

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

s = (0:nP^2*2-1)';
y0 = STx*s(1:end/2) + STy*s(1+end/2:end);

figure(33)
subplot(2,1,1)
plot([y0,y],'.')
subplot(2,1,2)
plot(y0-y)