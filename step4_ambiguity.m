%% ========================================================================
%  STEP 4 | Irregular-sampling operator & two-timescale ambiguity function
%  Loads step2_out.mat. Builds the slow-time sample set T for a target and
%  compares the velocity ambiguity function chi(v) for:
%     (U) uniform inter-burst  (spacing T)        -> alias comb
%     (A) irregular, ADJACENT beams (clustered)   -> realistic, weak suppress
%     (S) irregular, SPREAD beams (permuted sweep)-> strong suppress
%  VERIFY: chi factorizes as chi_inter*chi_intra, and the measured alias
%  envelope equals the closed form (1/C)|sum exp(j2*pi*p*delta/T)|.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:});
PFS={'FAIL','PASS'};   % PFS{1+cond}
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G; SC=L2.SC;
LOG('\n================ STEP 4: AMBIGUITY FUNCTION ================\n');

lam=P.lambda; T=P.T; M=P.M; L=P.L; W=D.W;
ktar=2; vt=SC.v(ktar); thk=SC.th(ktar);       % study target T2 (v=25, +10deg)
C=6;                                          % illuminating beams
LOG('\n[0] Studying T%d: v=%.1f m/s, theta=%.0f deg, C=%d beams, W=%.4f m/s\n',ktar,vt,thk,C,W);

%% ---- 1. Build the three intra-burst offset sets ----------------------
[~,bc]=min(abs(G.thetaBeam-thk));
covA = max(1,bc-floor(C/2)) : min(L,bc-floor(C/2)+C-1);   % adjacent beams (monotonic sweep)
offA = G.ssbTime(covA);                                   % clustered offsets
idxS = round(linspace(1,L,C));                            % spread time-slots (permuted sweep)
offS = G.ssbTime(idxS);                                   % spread offsets
LOG('\n[1] Intra-burst offset spans\n');
LOG('    ADJACENT beams %s -> span %.1f us  (off/T max = %.4f)\n', mat2str(covA), (max(offA)-min(offA))*1e6, max(offA-min(offA))/T);
LOG('    SPREAD  slots  %s -> span %.1f ms  (off/T max = %.4f)\n', mat2str(idxS), (max(offS)-min(offS))*1e3, max(offS-min(offS))/T);

%% ---- 2. Sample sets T = { m*T + offset } -----------------------------
m=(0:M-1).';
tU = m*T;                          % uniform inter-burst only (M x 1)
tA = m*T + offA;                   % M x C
tS = m*T + offS;                   % M x C

%% ---- 3. Ambiguity function chi(v) ------------------------------------
chi=@(vg,tset) abs(mean(exp(1j*2*pi*(2/lam)*(vt-vg(:)).*reshape(tset,1,[])),2));
vg = linspace(vt-12*W, vt+12*W, 4001).';
cU = chi(vg,tU); cA = chi(vg,tA); cS = chi(vg,tS);

%% ---- 4. VERIFY factorization chi = chi_inter * chi_intra (spread) -----
chiInter = abs(mean(exp(1j*2*pi*(2/lam)*(vt-vg(:)).*reshape(tU,1,[])),2));   % uses inter grid
chiIntraS= abs(mean(exp(1j*2*pi*(2/lam)*(vt-vg(:)).*reshape(offS,1,[])),2));
facErr = max(abs(cS - chiInter.*chiIntraS));
LOG('\n[4] Factorization check (spread): max|chi - chi_inter*chi_intra| = %.2e  [%s]\n', ...
     facErr, PFS{1+(facErr<1e-9)});

%% ---- 5. VERIFY alias envelope vs closed form -------------------------
LOG('\n[5] Alias envelope at teeth v = vt - p*W   (spread case)\n');
LOG('    %3s | %10s %10s %10s | %s\n','p','closedform','measured','uniform','verdict');
allok=true;
for p=0:4
  cf = abs(mean(exp(1j*2*pi*p*offS(:)/T)));        % closed-form intra envelope
  vp = vt - p*W;
  meas = chi(vp,tS);                                % measured on spread set
  uni  = chi(vp,tU);                                % uniform (should stay ~1)
  ok = abs(meas-cf)<1e-6; allok=allok&&ok;
  LOG('    %3d | %10.4f %10.4f %10.4f | %s\n',p,cf,meas,uni,PFS{1+ok});
end
LOG('    ENVELOPE MATCH: %s\n', PFS{1+allok});

% peak-to-worst-alias (exclude +/-0.5W around true peak)
psl=@(cv) 20*log10(max(cv(abs(vg-vt)>0.5*W))/max(cv));
LOG('\n    worst alias level (dB below true peak):\n');
LOG('      uniform  = %+6.2f dB   adjacent = %+6.2f dB   spread = %+6.2f dB\n', psl(cU),psl(cA),psl(cS));
LOG('    => uniform & adjacent: aliases ~0 dB (fully ambiguous).\n');
LOG('       spread breaks the degeneracy so Step 6 sparse solver can pick p=0.\n');

%% ---- 6. Plots --------------------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');

f1=figure('position',[0 0 780 420]);
plot(vg,cU,'-','color',[.85 .33 .1],'linewidth',1.1); hold on;
plot(vg,cA,'-','color',[.4 .4 .4],'linewidth',1.0);
plot(vg,cS,'-','color',[.1 .5 .8],'linewidth',1.6);
plot([vt vt],[0 1.05],'k--');
xlabel('candidate velocity [m/s]'); ylabel('|\chi(v)|'); ylim([0 1.08]);
legend('uniform (comb)','irregular adjacent','irregular spread','true v','location','south');
title('Fig1: ambiguity function - spread sampling breaks the alias comb');
print(f1,'step4_fig1_ambiguity.png','-dpng','-r110');

f2=figure('position',[0 0 620 340]);
pp=0:4; cfv=arrayfun(@(p) abs(mean(exp(1j*2*pi*p*offS(:)/T))),pp);
mev=arrayfun(@(p) chi(vt-p*W,tS),pp);
bar(pp, [cfv(:) mev(:)]); grid on;
xlabel('alias index p'); ylabel('|\chi_{intra}|');
legend('closed form','measured','location','northeast');
title('Fig2: alias envelope - closed form vs measured (spread)');
print(f2,'step4_fig2_envelope.png','-dpng','-r110');

save('step4_out.mat','vg','cU','cA','cS','offA','offS','covA','idxS','C','ktar','-v7');
LOG('\n[7] Saved step4_out.mat + 2 PNGs\n');
LOG('================ STEP 4 COMPLETE ================\n\n');
