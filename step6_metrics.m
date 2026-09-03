%% ========================================================================
%  STEP 6 | Performance metrics (Monte-Carlo) - refined estimator
%  A) velocity RMSE vs SNR  (verify approaches CRB)
%  B) P(correct de-aliasing) vs SNR
%  C) success vs number of co-cell targets K  (sparsity limit)
%  Estimator: top-candidate REFINED matching pursuit (defeats picket fence).
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:}); PFS={'FAIL','PASS'};
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G;
lam=P.lambda; T=P.T; M=P.M; L=P.L; W=D.W;
try, randn('seed',3); rand('seed',3); catch, end
LOG('\n================ STEP 6: PERFORMANCE METRICS ================\n');

C=6; idxS=round(linspace(1,L,C)); offS=G.ssbTime(idxS);
m=(0:M-1).'; tS=reshape(m*T+offS,[],1); Ns=numel(tS);
Stt=sum((tS-mean(tS)).^2); crb_v=@(eta)(lam/2)*sqrt(1./(2*eta*(2*pi)^2*Stt));
kk=2*pi*(2/lam);                                  % phase constant
vgc=(-30:0.01:30); PhiC=exp(1j*kk*tS*vgc);        % coarse dictionary
LOG('\n[0] Ns=%d  sum(t-tbar)^2=%.3e   grid=%.3f  refine=1e-4\n',Ns,Stt,vgc(2)-vgc(1));

%% ===== A/B: RMSE & P(correct) vs SNR  (K=1) ===========================
SNRdB=-10:5:30; Ntr=150; vspan=25;
rmse=zeros(size(SNRdB)); pcor=zeros(size(SNRdB));
LOG('\n[A/B] RMSE & P(correct) vs SNR  (K=1, %d trials/pt)\n',Ntr);
for is=1:numel(SNRdB)
  eta=10^(SNRdB(is)/10); sg=sqrt(1/(2*eta)); errs=[]; nc=0;
  for it=1:Ntr
    vt=-vspan+2*vspan*rand; ph=2*pi*rand;
    y=exp(1j*ph)*exp(1j*kk*vt*tS)+sg*(randn(Ns,1)+1j*randn(Ns,1));
    % ---- refined single-atom estimate ----
    mf=abs(PhiC'*y);
    lm=find(mf(2:end-1)>=mf(1:end-2)&mf(2:end-1)>=mf(3:end))+1;
    [~,o]=sort(mf(lm),'descend'); cand=lm(o(1:min(12,numel(o))));
    bb=-inf; vhat=vgc(cand(1));
    for cc=cand(:).'
      loc=vgc(cc)+(-0.02:1e-4:0.02);
      ml=abs(exp(1j*kk*tS*loc)'*y); [mv,li]=max(ml);
      if mv>bb, bb=mv; vhat=loc(li); end
    end
    e=vhat-vt; if abs(e)<0.15, nc=nc+1; errs(end+1)=e; end
  end
  rmse(is)=sqrt(mean(errs.^2)); pcor(is)=nc/Ntr;
  LOG('    SNR=%+3d : P(correct)=%.3f  RMSE=%.3g  CRB=%.3g\n',SNRdB(is),pcor(is),rmse(is),crb_v(eta));
end
i10=find(SNRdB==10); ratio=rmse(i10)/crb_v(10^(1));
LOG('    VERIFY RMSE/CRB at +10 dB = %.2f  [%s]\n',ratio,PFS{1+(ratio<2)});
LOG('    VERIFY P(correct) at +20 dB = %.3f  [%s]\n',pcor(SNRdB==20),PFS{1+(pcor(SNRdB==20)>0.95)});

%% ===== C: success vs K (refined OMP) ==================================
SNRc=20; eta=10^(SNRc/10); sg=sqrt(1/(2*eta)); Ks=1:8; NtrK=80; succ=zeros(size(Ks));
LOG('\n[C] Success vs co-cell K  (SNR=%d dB, %d trials/pt, refined OMP)\n',SNRc,NtrK);
for ik=1:numel(Ks)
  K=Ks(ik); ok=0;
  for it=1:NtrK
    vt=[]; while numel(vt)<K, cd=-vspan+2*vspan*rand; if isempty(vt)||min(abs(vt-cd))>0.5, vt(end+1)=cd; end, end
    b=(0.5+rand(K,1)).*exp(1j*2*pi*rand(K,1));
    y=zeros(Ns,1); for j=1:K, y=y+b(j)*exp(1j*kk*vt(j)*tS); end
    y=y+sg*(randn(Ns,1)+1j*randn(Ns,1));
    % refined matching pursuit, K atoms
    r=y; vrec=[];
    for q=1:K
      mf=abs(PhiC'*r);
      lm=find(mf(2:end-1)>=mf(1:end-2)&mf(2:end-1)>=mf(3:end))+1;
      if isempty(lm),[~,lm]=max(mf);end
      [~,o]=sort(mf(lm),'descend'); cand=lm(o(1:min(12,numel(o))));
      bb=-inf; vb=vgc(cand(1));
      for cc=cand(:).'
        loc=vgc(cc)+(-0.02:1e-4:0.02); ml=abs(exp(1j*kk*tS*loc)'*r);
        [mv,li]=max(ml); if mv>bb, bb=mv; vb=loc(li); end
      end
      vrec(end+1)=vb; A=exp(1j*kk*tS*vrec); xh=A\y; r=y-A*xh;
    end
    good=true; for j=1:K, if min(abs(vrec-vt(j)))>0.15, good=false; break; end, end
    ok=ok+good;
  end
  succ(ik)=ok/NtrK; LOG('    K=%d : success=%.3f\n',K,succ(ik));
end

%% ---- Plots ----------------------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off'); etaAx=10.^(SNRdB/10);
f1=figure('position',[0 0 620 420]);
semilogy(SNRdB,rmse,'o-','color',[.1 .5 .8],'linewidth',1.6,'markerfacecolor',[.6 .8 .95]); hold on;
semilogy(SNRdB,crb_v(etaAx),'k--','linewidth',1.3); grid on;
xlabel('SNR [dB]'); ylabel('velocity RMSE [m/s]'); legend('sparse estimator','CRB','location','southwest');
title('Fig1: velocity RMSE vs SNR (tracks CRB)'); print(f1,'step6_fig1_rmse.png','-dpng','-r110');
f2=figure('position',[0 0 620 380]);
plot(SNRdB,pcor,'s-','color',[.85 .33 .1],'linewidth',1.6,'markerfacecolor',[.95 .7 .5]); grid on;
xlabel('SNR [dB]'); ylabel('P(correct de-aliasing)'); ylim([0 1.05]);
title('Fig2: probability of correct fold vs SNR'); print(f2,'step6_fig2_pfold.png','-dpng','-r110');
f3=figure('position',[0 0 620 380]);
plot(Ks,succ,'^-','color',[.2 .5 .2],'linewidth',1.6,'markerfacecolor',[.6 .85 .6]); grid on;
xlabel('number of co-cell targets K'); ylabel('recovery success rate'); ylim([0 1.05]);
title('Fig3: success vs K - the sparsity limit'); print(f3,'step6_fig3_sparsity.png','-dpng','-r110');
save('step6_out.mat','SNRdB','rmse','pcor','Ks','succ','Stt','-v7');
LOG('\n[D] Saved step6_out.mat + 3 PNGs\n================ STEP 6 COMPLETE ================\n\n');
