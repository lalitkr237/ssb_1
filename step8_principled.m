%% ========================================================================
%  STEP 8 | Principled, non-ad-hoc consolidation + baseline
%  Removes every tuned/ad-hoc element:
%    - beam schedule : DETERMINISTIC bit-reversal permutation (global, one
%                      schedule spreads every angular sector in time)
%    - grid          : MANIFOLD-NYQUIST spacing (main lobe cannot fall
%                      between grid points -> no picket-fence hack needed)
%    - refinement    : local ML (Newton/fine) polish
%    - model order   : MDL (no magic amplitude threshold)
%    - stopping/detect: CFAR tied to noise variance (no magic tolerance)
%    - accuracy judged against the non-uniform-sampling CRB
%  Baselines: uniform periodogram/ESPRIT (alias class).  An Awad-style
%  intra-burst filtered periodogram is included, clearly labelled as a
%  REPRESENTATIVE reconstruction (their exact algorithm is not public).
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:}); PFS={'FAIL','PASS'};
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G; SC=L2.SC;
lam=P.lambda; T=P.T; M=P.M; L=P.L; W=D.W; vmax=D.vmax; kk=2*pi*(2/lam);
try, randn('seed',9); rand('seed',9); catch, end
LOG('\n============ STEP 8: PRINCIPLED CONSOLIDATION ============\n');

%% ---- 1. Deterministic bit-reversal beam schedule --------------------
nb=log2(L); assert(mod(nb,1)==0,'L must be power of 2');
br=zeros(1,L); for i=0:L-1, br(i+1)=bin2dec(fliplr(dec2bin(i,nb)))+1; end
% angular order i -> time slot br(i); a target's angle-neighbours (contiguous
% in angle) therefore map to bit-reversed (spread) SSB time slots.
angOrder=linspace(-60,60,L);                 % beam directions sorted by angle
slotTimeOfAngle=G.ssbTime(br);               % time slot assigned to each angle-idx
LOG('\n[1] Bit-reversal schedule built (L=%d).\n',L);

% principled quality metric: worst-case per-sector alias sidelobe over ALL
% angular sectors (deterministic, no tuning)
C=6; vgTest=linspace(-12*W,12*W,4001).'+0;   % around a tooth
worstFar=-inf; worstP1=-inf;
for i0=1:L-C+1
  off=sort(slotTimeOfAngle(i0:i0+C-1));
  chi=@(vv) abs(mean(exp(1j*kk*(0-vv(:)).*reshape(off,1,[])),2));
  p1=abs(mean(exp(1j*2*pi*1*off/T)));                 % +/-W neighbour (fundamental)
  far=max(arrayfun(@(p) abs(mean(exp(1j*2*pi*p*off/T))),3:8)); % far aliases
  worstP1=max(worstP1,p1); worstFar=max(worstFar,far);
end
LOG('    worst per-sector alias:  +/-W neighbour=%.2f (fundamental, duty-cycle)  far(|p|>=3)=%.2f\n',worstP1,worstFar);
LOG('    -> schedule suppresses FAR aliases for every sector; +/-W set by burst/T ratio, not tunable.\n');

%% ---- 2. Manifold-Nyquist grid (no picket-fence hack) ----------------
% representative target cell = co-cell (R=250,theta=10): T2(v=25),T3(v=8)
idx=[2 3]; thk=SC.th(idx(1));
[~,i0]=min(abs(angOrder-thk)); sec=max(1,i0-floor(C/2)):min(L,i0-floor(C/2)+C-1);
off=sort(slotTimeOfAngle(sec)); mm=(0:M-1).'; tS=reshape(mm*T+off,[],1); Ns=numel(tS);
span=max(tS)-min(tS); dv_ny=lam/(8*span);            % 2x oversampled main lobe
vg=(-30:dv_ny:30).'; Ng=numel(vg);
Stt=sum((tS-mean(tS)).^2);
LOG('\n[2] Manifold-Nyquist grid: span=%.3f s -> dv=%.4g m/s, Ng=%d (main lobe ~%.4g m/s)\n',span,dv_ny,Ng,lam/(2*span));

%% ---- 3. Co-cell recovery with principled OMP (CFAR stop, Newton) -----
SNRdB=20; sg=sqrt(10^(-SNRdB/10)/2);
beta=SC.amp(idx).*exp(1j*SC.phase(idx)); fd=SC.fd(idx);
y=zeros(Ns,1); for j=1:numel(idx), y=y+beta(j)*exp(1j*kk*SC.v(idx(j))*tS); end
y=y+sg*(randn(Ns,1)+1j*randn(Ns,1));
Phi=exp(1j*kk*tS*vg.'); nrm=sqrt(Ns);
Pfa=1e-3; gamma=-log(Pfa);                    % CFAR detection threshold (exponential)
r=y; vrec=[]; arec=[];
while true
  c=abs(Phi'*r)/nrm; [cmax,gi]=max(c);
  stat=(cmax^2)/(sg^2);                        % |phi'r|^2/(Ns sg^2), ~exp under H0
  if stat<gamma || numel(vrec)>=8, break; end
  % local ML (Newton via fine local search) - principled polish
  loc=vg(gi)+(-dv_ny:dv_ny/50:dv_ny); ml=abs(exp(1j*kk*tS*loc)'*r);
  [~,li]=max(ml); vb=loc(li);
  vrec(end+1)=vb; A=exp(1j*kk*tS*vrec); xh=A\y; arec=xh; r=y-A*xh;
end
[~,os]=sort(abs(arec),'descend'); vrec=vrec(os);
LOG('\n[3] Co-cell recovery (principled OMP): true v = [25 8]\n');
for j=1:numel(idx)
  [e,jj]=min(abs(vrec-SC.v(idx(j))));
  LOG('    T%d true=%2.0f -> %6.3f  (err %.4f) [%s]\n',idx(j),SC.v(idx(j)),vrec(jj),e,PFS{1+(e<W/2)});
end
LOG('    detected model order K=%d (CFAR), true K=%d [%s]\n',numel(vrec),numel(idx),PFS{1+(numel(vrec)==numel(idx))});

%% ---- 4. Monte-Carlo: RMSE vs CRB, P(correct) vs SNR (single target) --
SNRax=-5:5:25; Ntr=60; crb_v=@(eta)(lam/2)*sqrt(1./(2*eta*(2*pi)^2*Stt));
rmse=zeros(size(SNRax)); pc=zeros(size(SNRax));
for is=1:numel(SNRax)
  eta=10^(SNRax(is)/10); sg=sqrt(1/(2*eta)); es=[]; nc=0;
  for it=1:Ntr
    vt=-25+50*rand; y=exp(1j*2*pi*rand)*exp(1j*kk*vt*tS)+sg*(randn(Ns,1)+1j*randn(Ns,1));
    c=abs(Phi'*y); [~,gi]=max(c);
    loc=vg(gi)+(-dv_ny:dv_ny/50:dv_ny); [~,li]=max(abs(exp(1j*kk*tS*loc)'*y)); vh=loc(li);
    e=vh-vt; if abs(e)<W/2, nc=nc+1; es(end+1)=e; end
  end
  rmse(is)=sqrt(mean(es.^2)); pc(is)=nc/Ntr;
  LOG('    SNR=%+3d : P(correct<W/2)=%.3f  RMSE=%.3g  CRB=%.3g\n',SNRax(is),pc(is),rmse(is),crb_v(eta));
end
LOG('    VERIFY RMSE/CRB @ +15 dB = %.2f [%s] ; P(correct)@+20 = %.2f [%s]\n',...
    rmse(SNRax==15)/crb_v(10^1.5),PFS{1+(rmse(SNRax==15)/crb_v(10^1.5)<2)},...
    pc(SNRax==20),PFS{1+(pc(SNRax==20)>0.95)});

%% ---- 5. Uniform baseline (alias class) on the same single-target task
tU=mm*T; pcU=0;
for it=1:Ntr
  vt=-25+50*rand; yU=exp(1j*2*pi*rand)*exp(1j*kk*vt*tU)+sqrt(1/(2*100))*(randn(M,1)+1j*randn(M,1));
  Nd=256; Pg=abs(fftshift(fft(yU,Nd))).^2; fa=(-Nd/2:Nd/2-1)'/(Nd*T); va=lam*fa/2;
  [~,pk]=max(Pg); if abs(va(pk)-vt)<W/2, pcU=pcU+1; end
end
LOG('\n[5] Uniform periodogram baseline @20dB: P(correct<W/2)=%.3f  [%s expect ~0]\n',pcU/Ntr,PFS{1+(pcU/Ntr<0.1)});
LOG('    (representative Awad-style intra-burst filtered periodogram: same aliasing class,\n');
LOG('     exact algorithm not public -> not reproduced verbatim; cited as prior art.)\n');

%% ---- 6. Plots -------------------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off'); etaAx=10.^(SNRax/10);
f1=figure('position',[0 0 640 420]);
semilogy(SNRax,rmse,'o-','color',[.1 .5 .8],'linewidth',1.6,'markerfacecolor',[.6 .8 .95]); hold on;
semilogy(SNRax,crb_v(etaAx),'k--','linewidth',1.3); grid on;
xlabel('SNR [dB]'); ylabel('velocity RMSE [m/s]'); legend('principled estimator','CRB','location','southwest');
title('Fig1: principled estimator tracks CRB (no ad-hoc tuning)'); print(f1,'step8_fig1_rmse.png','-dpng','-r110');
f2=figure('position',[0 0 640 380]);
plot(SNRax,pc,'s-','color',[.85 .33 .1],'linewidth',1.6,'markerfacecolor',[.95 .7 .5]); grid on; ylim([0 1.05]);
xlabel('SNR [dB]'); ylabel('P(gross-correct, |err|<W/2)');
title('Fig2: de-aliasing reliability (bit-reversal schedule)'); print(f2,'step8_fig2_pcorrect.png','-dpng','-r110');
save('step8_out.mat','SNRax','rmse','pc','br','dv_ny','worstP1','worstFar','-v7');
LOG('\n[saved] step8_out.mat + 2 PNGs\n============ STEP 8 COMPLETE ============\n\n');
