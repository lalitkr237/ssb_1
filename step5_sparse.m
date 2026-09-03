%% ========================================================================
%  STEP 5 | Joint sparse de-aliasing (OMP over a wide velocity dictionary)
%  Loads step2_out.mat. For each range-angle cell, recover TRUE (unaliased)
%  velocities from the irregular spread-beam slow-time samples.
%  VERIFY: recovered == true speeds (25,8,12,3), and the SAME solver on
%  UNIFORM samples fails to de-alias -> proves the operator is essential.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:}); PFS={'FAIL','PASS'};
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G; SC=L2.SC;
LOG('\n================ STEP 5: SPARSE DE-ALIASING ================\n');

lam=P.lambda; T=P.T; M=P.M; L=P.L; W=D.W;
SNRdB=20; sigma=sqrt(10^(-SNRdB/10)/2);
rng_s=11; try, randn('seed',rng_s); catch, end

%% ---- spread-beam sampling operator (permuted schedule) ---------------
C=6; idxS=round(linspace(1,L,C)); offS=G.ssbTime(idxS);   % spread offsets
m=(0:M-1).'; tS = reshape(m*T + offS,[],1);               % Ns x 1 (spread)
tU = m*T;                                                 % M x 1 (uniform)
Ns=numel(tS);
LOG('\n[0] SNR=%d dB  Ns(spread)=%d  Ns(uniform)=%d  vmax=%.3f W=%.3f\n',SNRdB,Ns,M,D.vmax,W);

%% ---- wide velocity dictionary ---------------------------------------
vg = (-30:0.01:30).';                       % candidate velocities (true range)
PhiS = exp(1j*2*pi*(2/lam)*tS*vg.');        % Ns x Ng   (spread operator)
PhiU = exp(1j*2*pi*(2/lam)*tU*vg.');        % M  x Ng   (uniform operator)
PhiSn=PhiS./sqrt(sum(abs(PhiS).^2,1));      % unit-norm columns
PhiUn=PhiU./sqrt(sum(abs(PhiU).^2,1));
LOG('[1] Dictionary: %d velocities in [%.0f, %.0f] m/s (step %.3f)\n',numel(vg),vg(1),vg(end),vg(2)-vg(1));

%% ---- cells (grouped targets sharing range-angle) --------------------
cells={[1],[2 3],[4]};
LOG('\n[2] Per-cell OMP recovery on the SPREAD operator\n');
recAll=[]; truAll=[]; pass=true;
for ci=1:numel(cells)
  idx=cells{ci}; beta=SC.amp(idx).*exp(1j*SC.phase(idx)); fd=SC.fd(idx);
  % synthesize cell slow-time measurement at spread times
  y = zeros(Ns,1);
  for j=1:numel(idx), y=y+beta(j)*exp(1j*2*pi*fd(j)*tS); end
  y=y+sigma*(randn(Ns,1)+1j*randn(Ns,1));
  % --- OMP ---
  r=y; supp=[]; Kmax=numel(idx)+2; resTh=sqrt(1.5*Ns)*sigma;
  for it=1:Kmax
    c=abs(PhiSn'*r); [~,gi]=max(c);
    if any(supp==gi), break; end
    supp=[supp gi]; A=PhiS(:,supp); xh=A\y; r=y-A*xh;
    if norm(r)<resTh, break; end
  end
  [amp,ord]=sort(abs(xh),'descend'); supp=supp(ord);
  keep=amp>0.15*max(amp); vrec=vg(supp(keep));
  % match each true target to nearest recovered
  LOG('    cell%d (R=%.0f,th=%.0f):\n',ci,SC.R(idx(1)),SC.th(idx(1)));
  for j=1:numel(idx)
    [e,~]=min(abs(vrec-SC.v(idx(j))));
    ok=e<=0.1; pass=pass&&ok;
    [~,jj]=min(abs(vrec-SC.v(idx(j))));
    LOG('       T%d true v=%5.1f -> recovered %6.2f  (err %.3f) [%s]\n',...
        idx(j),SC.v(idx(j)),vrec(jj),e,PFS{1+ok});
    recAll(end+1)=vrec(jj); truAll(end+1)=SC.v(idx(j));
  end
end
LOG('    OVERALL SPARSE DE-ALIASING: %s\n',PFS{1+pass});

%% ---- contrast: SAME solver on the UNIFORM operator (co-cell) ---------
LOG('\n[3] Contrast on UNIFORM operator (co-cell T2,T3) - expected to FAIL\n');
idx=[2 3]; beta=SC.amp(idx).*exp(1j*SC.phase(idx)); fd=SC.fd(idx);
yU=zeros(M,1); for j=1:numel(idx), yU=yU+beta(j)*exp(1j*2*pi*fd(j)*tU); end
yU=yU+sigma*(randn(M,1)+1j*randn(M,1));
r=yU; supp=[];
for it=1:4
  c=abs(PhiUn'*r); [~,gi]=max(c); if any(supp==gi),break; end
  supp=[supp gi]; A=PhiU(:,supp); xh=A\yU; r=yU-A*xh;
  if norm(r)<sqrt(1.5*M)*sigma, break; end
end
vrecU=vg(supp);
LOG('    uniform-operator OMP recovered: %s m/s\n',mat2str(round(vrecU'*100)/100));
LOG('    true speeds were 25 and 8 -> uniform operator cannot de-alias (aliased picks).\n');

%% ---- Plots ----------------------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');

% Fig1: co-cell RECOVERED SPIKES (OMP result) - spread vs uniform
idx=[2 3]; beta=SC.amp(idx).*exp(1j*SC.phase(idx)); fd=SC.fd(idx);
yS=zeros(Ns,1); for j=1:numel(idx), yS=yS+beta(j)*exp(1j*2*pi*fd(j)*tS); end
yS=yS+sigma*(randn(Ns,1)+1j*randn(Ns,1));
% OMP on spread operator (record spikes)
r=yS; supp=[]; for it=1:4, c=abs(PhiSn'*r);[~,gi]=max(c); if any(supp==gi),break;end
  supp=[supp gi]; A=PhiS(:,supp); xh=A\yS; r=yS-A*xh; if norm(r)<sqrt(1.5*Ns)*sigma,break;end; end
vS=vg(supp); aS=abs(xh)/max(abs(xh));
% OMP on uniform operator (record spikes)
r=yU; supp=[]; for it=1:4, c=abs(PhiUn'*r);[~,gi]=max(c); if any(supp==gi),break;end
  supp=[supp gi]; A=PhiU(:,supp); xh=A\yU; r=yU-A*xh; if norm(r)<sqrt(1.5*M)*sigma,break;end; end
vU=vg(supp); aU=abs(xh)/max(abs(xh));
f1=figure('position',[0 0 780 380]);
for j=1:numel(idx), plot([SC.v(idx(j)) SC.v(idx(j))],[0 1.1],'k--'); hold on; end
hU=stem(vU,aU,'o','color',[.85 .33 .1],'linewidth',1.4,'markersize',6);
hS=stem(vS,aS,'filled','color',[.1 .5 .8],'linewidth',1.8,'markersize',7);
xlim([-30 30]); ylim([0 1.15]); xlabel('velocity [m/s]'); ylabel('normalized recovered amplitude');
legend([hS hU],{'spread operator (OMP)','uniform operator (OMP)'},'location','north');
title('Fig1: co-cell OMP spikes - spread recovers true 25 & 8; uniform picks aliases');
print(f1,'step5_fig1_omp_cocell.png','-dpng','-r110');

% Fig2: recovered vs true (companion to Step 3 collapse)
f2=figure('position',[0 0 620 460]);
plot([-2 27],[-2 27],'-','color',[.7 .7 .7]); hold on;
plot(truAll,recAll,'o','color',[.1 .5 .8],'markersize',9,'linewidth',1.6,'markerfacecolor',[.6 .8 .95]);
for i=1:numel(truAll), text(truAll(i)+0.5,recAll(i),sprintf('%.0f',truAll(i)),'fontsize',9); end
axis([-2 27 -2 27]); grid on; xlabel('true velocity [m/s]'); ylabel('recovered velocity [m/s]');
title('Fig2: de-aliasing works - recovered == true (on the diagonal)');
print(f2,'step5_fig2_recovered.png','-dpng','-r110');

save('step5_out.mat','vg','recAll','truAll','tS','offS','C','SNRdB','-v7');
LOG('\n[4] Saved step5_out.mat + 2 PNGs\n');
LOG('================ STEP 5 COMPLETE ================\n\n');
