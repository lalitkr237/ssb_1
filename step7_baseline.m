%% ========================================================================
%  STEP 7 | Baseline + complexity comparison
%  Methods on the SAME single-target scenes vs SNR (Nsnap snapshots):
%    (1) Periodogram  (uniform inter-burst)  - ALIASES
%    (2) ESPRIT       (uniform inter-burst)  - ALIASES
%    (3) 2D-MUSIC     (spread operator)      - de-aliases
%    (4) Proposed     refined MP (spread)    - de-aliases
%  Outputs: RMSE vs SNR (+CRB), P(correct) vs SNR, Big-O table + runtimes.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:}); PFS={'FAIL','PASS'};
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G;
lam=P.lambda; T=P.T; M=P.M; L=P.L; W=D.W; vmax=D.vmax;
try, randn('seed',5); rand('seed',5); catch, end
LOG('\n================ STEP 7: BASELINE + COMPLEXITY ================\n');

kk=2*pi*(2/lam);
C=6; idxS=round(linspace(1,L,C)); offS=G.ssbTime(idxS);
mm=(0:M-1).'; tS=reshape(mm*T+offS,[],1); Ns=numel(tS); tU=mm*T;
Stt=sum((tS-mean(tS)).^2); Nsnap=32;
crb_v=@(eta)(lam/2)*sqrt(1./(2*Nsnap*eta*(2*pi)^2*Stt));
vgc=(-30:0.01:30); Ng=numel(vgc);
PhiC=exp(1j*kk*tS*vgc);          % Ns x Ng  (spread manifold)
Nd=512;                          % periodogram FFT size
LOG('\n[0] Ns=%d M=%d Nsnap=%d Ng=%d  vmax=%.3f W=%.3f\n',Ns,M,Nsnap,Ng,vmax,W);

SNRdB=-10:5:25; Ntr=100; vspan=25;
mList={'Periodogram-U','ESPRIT-U','MUSIC-S','Proposed-S'};
RMSE=zeros(4,numel(SNRdB)); PC=zeros(4,numel(SNRdB));

for is=1:numel(SNRdB)
  eta=10^(SNRdB(is)/10); sg=sqrt(1/(2*eta));
  esum=zeros(4,1); nc=zeros(4,1);
  for it=1:Ntr
    vt=-vspan+2*vspan*rand;
    aU=exp(1j*kk*vt*tU); aS=exp(1j*kk*vt*tS);
    bset=exp(1j*2*pi*rand(1,Nsnap));                 % random snapshot gains
    YU=aU*bset + sg*(randn(M,Nsnap)+1j*randn(M,Nsnap));
    YS=aS*bset + sg*(randn(Ns,Nsnap)+1j*randn(Ns,Nsnap));

    % (1) Periodogram uniform
    Pg=mean(abs(fftshift(fft(YU,Nd,1),1)).^2,2);
    fa=(-Nd/2:Nd/2-1)'/(Nd*T); va=lam*fa/2; [~,pi1]=max(Pg); v1=va(pi1);

    % (2) ESPRIT uniform (K=1)
    Ru=YU*YU'/Nsnap; [Eu,Du]=eig(Ru); [~,oi]=max(real(diag(Du))); es=Eu(:,oi);
    psi=es(1:end-1)\es(2:end); v2=lam*(angle(psi)/(2*pi*T))/2;

    % (3) MUSIC spread (K=1): top-candidate refined search (fair vs proposed)
    Rs=YS*YS'/Nsnap; [Es,Ds]=eig(Rs); [~,si]=sort(real(diag(Ds)),'descend');
    En=Es(:,si(2:end));
    music=1./sum(abs(En'*PhiC).^2,1);
    lm=find(music(2:end-1)>=music(1:end-2)&music(2:end-1)>=music(3:end))+1;
    [~,o]=sort(music(lm),'descend'); cand=lm(o(1:min(10,numel(o))));
    bb=-inf; v3=vgc(cand(1));
    for cc=cand(:).'
      loc=vgc(cc)+(-0.005:1e-5:0.005); Al=exp(1j*kk*tS*loc);
      msp=1./sum(abs(En'*Al).^2,1); [mv,li]=max(msp);
      if mv>bb, bb=mv; v3=loc(li); end
    end

    % (4) Proposed: refined incoherent MP on spread (no EVD)
    mf=sum(abs(PhiC'*YS).^2,2);
    lm=find(mf(2:end-1)>=mf(1:end-2)&mf(2:end-1)>=mf(3:end))+1;
    [~,o]=sort(mf(lm),'descend'); cand=lm(o(1:min(10,numel(o))));
    bb=-inf; v4=vgc(cand(1));
    for cc=cand(:).'
      loc=vgc(cc)+(-0.005:1e-5:0.005); ml=sum(abs(exp(1j*kk*tS*loc)'*YS).^2,2);
      [mv,li]=max(ml); if mv>bb, bb=mv; v4=loc(li); end
    end

    vh=[v1 v2 v3 v4];
    for q=1:4, e=vh(q)-vt; esum(q)=esum(q)+e^2; if abs(e)<0.15, nc(q)=nc(q)+1; end, end
  end
  RMSE(:,is)=sqrt(esum/Ntr); PC(:,is)=nc/Ntr;
  LOG('  SNR=%+3d | Pgram P=%.2f ESPRIT P=%.2f | MUSIC P=%.2f RMSE=%.2g | Prop P=%.2f RMSE=%.2g (CRB %.2g)\n',...
      SNRdB(is),PC(1,is),PC(2,is),PC(3,is),RMSE(3,is),PC(4,is),RMSE(4,is),crb_v(eta));
end

LOG('\n[verify] high-SNR:\n');
LOG('   uniform methods de-alias?  Pgram P=%.2f ESPRIT P=%.2f  [%s expect ~0]\n',PC(1,end),PC(2,end),PFS{1+(PC(1,end)<0.1&&PC(2,end)<0.1)});
LOG('   spread methods de-alias?   MUSIC P=%.2f Prop P=%.2f     [%s]\n',PC(3,end),PC(4,end),PFS{1+(PC(3,end)>0.95&&PC(4,end)>0.95)});
LOG('   proposed tracks CRB @20dB: RMSE/CRB=%.2f [%s]\n',RMSE(4,SNRdB==20)/crb_v(10^2),PFS{1+(RMSE(4,SNRdB==20)/crb_v(10^2)<2)});

%% ---- Complexity: theory + measured runtime --------------------------
LOG('\n[complexity] per-estimate Big-O:\n');
LOG('   Periodogram-U : O(Nsnap*Nd*log Nd)            no de-alias\n');
LOG('   ESPRIT-U      : O(M^3 + Nsnap*M^2)            no de-alias\n');
LOG('   MUSIC-S       : O(Ns^3 + Nsnap*Ns^2 + Ng*Ns*(Ns-K))  de-alias\n');
LOG('   Proposed-S    : O(Nsnap*Ng*Ns + Kcand*R*Ns)  de-alias (no EVD)\n');
% measured runtime at SNR=10
eta=10^1; sg=sqrt(1/(2*eta)); Nrep=60; tmr=zeros(4,1);
for q=1:4
 tic;
 for it=1:Nrep
   vt=-vspan+2*vspan*rand; aU=exp(1j*kk*vt*tU); aS=exp(1j*kk*vt*tS);
   bset=exp(1j*2*pi*rand(1,Nsnap));
   YU=aU*bset+sg*(randn(M,Nsnap)+1j*randn(M,Nsnap));
   YS=aS*bset+sg*(randn(Ns,Nsnap)+1j*randn(Ns,Nsnap));
   if q==1, Pg=mean(abs(fft(YU,Nd,1)).^2,2); [~,~]=max(Pg);
   elseif q==2, Ru=YU*YU'/Nsnap; [Eu,Du]=eig(Ru); [~,oi]=max(real(diag(Du))); es=Eu(:,oi); psi=es(1:end-1)\es(2:end);
   elseif q==3, Rs=YS*YS'/Nsnap; [Es,Ds]=eig(Rs); [~,si]=sort(real(diag(Ds)),'descend'); En=Es(:,si(2:end)); music=1./sum(abs(En'*PhiC).^2,1); [~,~]=max(music);
   else, mf=sum(abs(PhiC'*YS).^2,2); [~,~]=max(mf);
   end
 end
 tmr(q)=toc/Nrep*1e3;
end
LOG('\n   measured runtime [ms/estimate] (Ng=%d, Ns=%d, Nsnap=%d):\n',Ng,Ns,Nsnap);
for q=1:4, LOG('     %-14s : %8.3f ms\n',mList{q},tmr(q)); end
LOG('   speedup proposed vs MUSIC = %.1fx\n',tmr(3)/tmr(4));

%% ---- Plots ----------------------------------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off'); etaAx=10.^(SNRdB/10);
cols=[.6 .6 .6; .85 .33 .1; .2 .5 .2; .1 .5 .8]; mk='osd^';
f1=figure('position',[0 0 640 440]);
for q=1:4, semilogy(SNRdB,RMSE(q,:),['-' mk(q)],'color',cols(q,:),'linewidth',1.5,'markerfacecolor',cols(q,:)); hold on; end
semilogy(SNRdB,crb_v(etaAx),'k--','linewidth',1.3); grid on;
xlabel('SNR [dB]'); ylabel('velocity RMSE [m/s]'); legend([mList,{'CRB'}],'location','southwest');
title('Fig1: RMSE vs SNR - uniform methods alias (flat, high)'); print(f1,'step7_fig1_rmse.png','-dpng','-r110');
f2=figure('position',[0 0 640 400]);
for q=1:4, plot(SNRdB,PC(q,:),['-' mk(q)],'color',cols(q,:),'linewidth',1.5,'markerfacecolor',cols(q,:)); hold on; end
grid on; ylim([0 1.05]); xlabel('SNR [dB]'); ylabel('P(correct de-aliasing)'); legend(mList,'location','east');
title('Fig2: P(correct) - only spread methods de-alias'); print(f2,'step7_fig2_pfold.png','-dpng','-r110');
f3=figure('position',[0 0 640 400]);
barh(1:4,tmr); set(gca,'ytick',1:4,'yticklabel',mList); grid on; xlabel('runtime [ms/estimate]');
title('Fig3: complexity - measured runtime'); print(f3,'step7_fig3_runtime.png','-dpng','-r110');
save('step7_out.mat','SNRdB','RMSE','PC','tmr','mList','-v7');
LOG('\n[saved] step7_out.mat + 3 PNGs\n================ STEP 7 COMPLETE ================\n\n');
