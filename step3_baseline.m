%% ========================================================================
%  STEP 3 | Baseline range-Doppler periodogram (the aliasing failure)
%  Loads step2_out.mat. Standard uniform inter-burst slow-time FFT.
%  VERIFY: detected velocity == folded vmeas (NOT true v) -> proves aliasing.
%  Produces the "Figure 1" of the Letter.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:});
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G; A=L2.A; SC=L2.SC; Y=L2.Y;
LOG('\n================ STEP 3: BASELINE PERIODOGRAM ================\n');

Nsc=P.Nsc; M=P.M; L=P.L;

%% ---- 1. Range compression (IFFT over subcarriers) --------------------
Nr=1024; rax=(0:Nr-1)/Nr*D.Rmax;
Rcube=zeros(Nr,M,L);
for mm=1:M, for ll=1:L, Rcube(:,mm,ll)=ifft(Y(:,mm,ll),Nr); end, end
LOG('\n[1] Range-compressed cube: %dx%dx%d (range x m x beam)\n',Nr,M,L);

%% ---- 2. Doppler / velocity axis (uniform inter-burst, spacing T) ------
Nd=256; fax=(-Nd/2:Nd/2-1)/(Nd*P.T); vax=P.lambda*fax/2;   % spans +/- vmax
LOG('[2] Velocity axis: %d bins over [%.4f, %.4f] m/s (vmax=%.4f)\n', ...
     Nd, vax(1), vax(end), D.vmax);

%% ---- 3. Detect each target: peak velocity at its (range,beam) --------
LOG('\n[3] Per-target baseline detection (peak-matched; handles co-cell)\n');
LOG('    %-3s %6s %10s %10s %10s | verdict\n','id','v_true','v_meas*','v_detect','|det-meas|');
pass=true;
for k=1:SC.K
  [~,pR]=min(abs(rax-SC.R(k)));            % target range bin
  [~,lb]=min(abs(G.thetaBeam-SC.th(k)));   % beam pointing at target
  s = squeeze(Rcube(pR,:,lb));             % slow-time vector (1 x M) - may hold >1 target
  Sd= abs(fftshift(fft(s,Nd)));
  % local maxima above 30% of peak (resolves both co-cell folds)
  ispk=[false, (Sd(2:end-1)>Sd(1:end-2)) & (Sd(2:end-1)>Sd(3:end)), false] & (Sd>0.3*max(Sd));
  vpk=vax(ispk);
  [err,j]=min(abs(vpk-SC.vmeas(k))); vdet=vpk(j);   % peak nearest this target's fold
  ok = err <= 2*D.dv; pass=pass&&ok;
  vd='FAIL'; if ok, vd='PASS'; end
  LOG('    T%-2d %6.1f %10.4f %10.4f %10.4f | %s (true err=%.2f)\n', ...
       k, SC.v(k), SC.vmeas(k), vdet, err, vd, abs(vdet-SC.v(k)));
end
sAll='FAIL'; if pass, sAll='PASS'; end
LOG('    OVERALL (detected == folded vmeas): %s\n', sAll);
LOG('    => baseline reports all targets within +/-%.3f m/s of ZERO,\n', D.vmax);
LOG('       i.e. true speeds 3..25 m/s are LOST to aliasing.\n');

%% ---- 4. Fig 1: range-velocity map at beam=+10deg (T2 & T3) -----------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');
[~,lb10]=min(abs(G.thetaBeam-10));
RV=zeros(Nr,Nd);
for p=1:Nr, RV(p,:)=abs(fftshift(fft(squeeze(Rcube(p,:,lb10)),Nd))); end
RV=RV/max(RV(:));
f1=figure('position',[0 0 720 420]);
imagesc(vax,rax,20*log10(RV+1e-6)); axis xy; caxis([-30 0]);
xlabel('velocity [m/s]'); ylabel('range [m]'); ylim([230 270]); colorbar; hold on;
plot(SC.vmeas(2),SC.R(2),'w+','markersize',11,'linewidth',1.6);
plot(SC.vmeas(3),SC.R(3),'wx','markersize',11,'linewidth',1.6);
title('Fig1: baseline RV @10deg - T2(+),T3(x) fold into +/-0.134 m/s window');
print(f1,'step3_fig1_rvmap.png','-dpng','-r110');

%% ---- 5. Fig 2: the collapse - true vs baseline-measured velocity ------
f2=figure('position',[0 0 720 360]);
for k=1:SC.K
  plot([SC.v(k) SC.vmeas(k)],[k k],'-','color',[.7 .7 .7]); hold on;
  plot(SC.v(k),k,'o','color',[.1 .45 .8],'markersize',8,'linewidth',1.5);
  plot(SC.vmeas(k),k,'s','color',[.85 .33 .1],'markersize',8,'linewidth',1.5,'markerfacecolor',[.85 .33 .1]);
end
plot([D.vmax D.vmax],[0 SC.K+1],'r:'); plot([-D.vmax -D.vmax],[0 SC.K+1],'r:');
xlim([-2 27]); ylim([0 SC.K+1]); set(gca,'ytick',1:SC.K,'yticklabel',{'T1','T2','T3','T4'});
xlabel('velocity [m/s]'); grid on;
title('Fig2: TRUE (o) vs baseline-MEASURED (square) - all collapse into red window');
print(f2,'step3_fig2_collapse.png','-dpng','-r110');

%% ---- 6. Save ---------------------------------------------------------
save('step3_out.mat','vax','fax','rax','Nr','Nd','-v7');
LOG('\n[6] Saved step3_out.mat + 2 PNGs\n');
LOG('================ STEP 3 COMPLETE ================\n\n');
