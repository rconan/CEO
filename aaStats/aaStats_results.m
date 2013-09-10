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
