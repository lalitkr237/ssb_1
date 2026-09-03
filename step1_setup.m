%% ========================================================================
%  STEP 1  |  Parameters & SSB sampling geometry
%  SSB-based ISAC velocity de-aliasing study.  MATLAB or Octave. No toolboxes.
%
%    (a) define physical/numerology parameters
%    (b) compute derived sensing limits and SELF-CHECK vs verified derivation
%    (c) build SSB slow-time sampling geometry (the two timescales)
%    (d) plot diagnostics to eyeball the geometry for bugs
%    (e) save to step1_out.mat for later steps
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG = @(varargin) fprintf(varargin{:});
LOG('\n================ STEP 1: PARAMETERS & SSB GEOMETRY ================\n');

%% ---- 1. Physical constants & numerology (FR2) --------------------------
P.c=3e8; P.fc=28e9; P.lambda=P.c/P.fc;
P.mu=3; P.scs=15e3*2^P.mu; P.Nsc=240; P.Nsym=4; P.L=64;
P.T=20e-3; P.M=16; P.Tsym=1e-3/(14*2^P.mu);
LOG('\n[1] Numerology\n');
LOG('    fc   = %.3g GHz   lambda = %.4f mm\n', P.fc/1e9, P.lambda*1e3);
LOG('    SCS  = %g kHz     Tsym   = %.4f us\n', P.scs/1e3, P.Tsym*1e6);
LOG('    Nsc=%d  Nsym/SSB=%d  L=%d  T=%.0f ms  M=%d\n', P.Nsc,P.Nsym,P.L,P.T*1e3,P.M);

%% ---- 2. Derived sensing limits + SELF-CHECK ---------------------------
D.B=P.Nsc*P.scs; D.dR=P.c/(2*D.B); D.Rmax=P.c/(2*P.scs);
D.W=P.lambda/(2*P.T); D.vmax=P.lambda/(4*P.T); D.dv=P.lambda/(2*P.M*P.T);
D.dSSB=4*P.Tsym; D.vmaxIntra=P.lambda/(4*D.dSSB);
LOG('\n[2] Derived limits\n');
LOG('    B      = %6.2f MHz -> dR   = %7.4f m   Rmax = %6.1f m\n', D.B/1e6,D.dR,D.Rmax);
LOG('    W=l/2T = %7.5f m/s   vmax = +/-%7.5f m/s   dv = %7.5f m/s\n', D.W,D.vmax,D.dv);
LOG('    intra SSB spacing = %.2f us -> intra vmax = %.2f m/s\n', D.dSSB*1e6,D.vmaxIntra);

names={'dR','Rmax','W','vmax','dv','vmaxIntra'};
vals =[D.dR,D.Rmax,D.W,D.vmax,D.dv,D.vmaxIntra];
refs =[5.2083,1250,0.26786,0.13393,0.016741,75.08];
tols =[1e-3,1e-3,1e-3,1e-3,1e-3,2e-3];
LOG('\n    --- self-check vs verified derivation ---\n');
allpass=true;
for i=1:numel(names)
  ok = abs(vals(i)-refs(i))<=tols(i)*max(1,abs(refs(i)));
  st='FAIL'; if ok, st='PASS'; else, allpass=false; end
  LOG('    check %-9s = %10.5f  [%s]\n', names{i}, vals(i), st);
end
LOG('    OVERALL: %s\n', mat2str(allpass));

%% ---- 3. SSB first-symbol indices (FR2 Case D, 120 kHz) ----------------
base=[4 8 16 20]; nset=[0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18];
symIdx=[]; for n=nset, symIdx=[symIdx, base+28*n]; end
symIdx=sort(symIdx);
assert(numel(symIdx)==P.L,'expected 64 SSB positions');
G.symIdx=symIdx; G.ssbTime=symIdx*P.Tsym;
fitstr='yes'; if G.ssbTime(end)>=5e-3, fitstr='NO'; end
LOG('\n[3] SSB timing (Case D)\n');
LOG('    #SSB positions       = %d\n', numel(symIdx));
LOG('    first / last start   = %.2f us / %.3f ms  (fits 5 ms: %s)\n', ...
     G.ssbTime(1)*1e6, G.ssbTime(end)*1e3, fitstr);
LOG('    early spacings       = %.1f / %.1f / %.1f us\n', ...
     diff(G.ssbTime(1:2))*1e6, diff(G.ssbTime(2:3))*1e6, diff(G.ssbTime(3:4))*1e6);

%% ---- 4. Beam sweep mapping -------------------------------------------
G.thetaBeam=linspace(-60,60,P.L); G.beamBW=120/P.L;
LOG('\n[4] Beam sweep: [%.0f,%.0f] deg over %d beams (%.2f deg/beam)\n', ...
     G.thetaBeam(1),G.thetaBeam(end),P.L,G.beamBW);

%% ---- 5. Full slow-time sample grid t(m,ell) --------------------------
m=(0:P.M-1).'; G.tGrid=m*P.T + G.ssbTime;
LOG('\n[5] Slow-time grid  %d x %d (MxL)   total span = %.2f ms\n', ...
     size(G.tGrid,1),size(G.tGrid,2),(G.tGrid(end,end)-G.tGrid(1,1))*1e3);

%% ---- 6. Example per-target sampling set ------------------------------
theta_k=6; C=6; [~,bc]=min(abs(G.thetaBeam-theta_k));
cover=max(1,bc-floor(C/2)):min(P.L,bc-floor(C/2)+C-1); G.coverIdx=cover;
tset=G.tGrid(:,cover); spanIntra=max(G.ssbTime(cover))-min(G.ssbTime(cover));
LOG('\n[6] Target @ %.0f deg: beams=%s\n', theta_k, mat2str(cover));
LOG('    intra span S=%.1f us -> coarse dv ~ %.1f m/s ; Ns=M*C=%d samples\n', ...
     spanIntra*1e6, P.lambda/(2*max(spanIntra,eps)), numel(tset));

%% ---- 7. SSB resource-grid schematic ----------------------------------
RG=zeros(P.Nsc,P.Nsym); ctr=56:182;
RG(ctr,1)=1; RG(:,2)=3; RG(1:4:end,2)=4;
RG(ctr,3)=2; RG([1:47,193:240],3)=3; RG(:,4)=3; RG(1:4:end,4)=4;
G.RG=RG;

%% ---- 8. Diagnostic plots ---------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');

f1=figure('position',[0 0 520 460]);
imagesc(1:P.Nsym,1:P.Nsc,RG); axis xy; colormap(jet(5)); caxis([0 4]); colorbar;
xlabel('OFDM symbol in SSB'); ylabel('subcarrier');
title('Fig 1: SSB grid  0 empty 1 PSS 2 SSS 3 PBCH 4 DMRS');
print(f1,'step1_fig1_grid.png','-dpng','-r110');

f2=figure('position',[0 0 720 340]);
for i=1:P.L, plot([G.ssbTime(i) G.ssbTime(i)]*1e3,[-62 G.thetaBeam(i)],'-','color',[.82 .82 .82]); hold on; end
scatter(G.ssbTime*1e3, G.thetaBeam, 26, G.thetaBeam,'filled'); colormap(jet);
xlabel('time within burst set [ms]'); ylabel('beam azimuth [deg]');
title('Fig 2: SSB sweep - non-uniform intra-burst spacing'); grid on;
print(f2,'step1_fig2_sweep.png','-dpng','-r110');

f3=figure('position',[0 0 760 340]);
plot(G.tGrid(:)*1e3, zeros(numel(G.tGrid),1),'.','color',[.7 .7 .7],'markersize',6); hold on;
plot(tset(:)*1e3, 0.05+zeros(numel(tset),1),'o','color',[.85 .33 .1],'markersize',5,'linewidth',1);
yl=[-0.3 0.4]; for mm=0:P.M-1, plot([mm*P.T mm*P.T]*1e3, yl,'-','color',[.92 .92 .92]); end
ylim(yl); set(gca,'ytick',[0 0.05],'yticklabel',{'all SSBs','target set'});
xlabel('absolute slow time [ms]');
title('Fig 3: two timescales - intra clusters repeated every T'); grid on;
print(f3,'step1_fig3_slowtime.png','-dpng','-r110');

%% ---- 9. Save ----------------------------------------------------------
save('step1_out.mat','P','D','G');
LOG('\n[9] Saved step1_out.mat + 3 PNGs\n');
LOG('================ STEP 1 COMPLETE ================\n\n');
