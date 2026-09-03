%% ========================================================================
%  STEP 9 | Feasibility theory + time-varying schedule (the proper result)
%  Derives & verifies the SSB velocity de-aliasing ALIAS HIERARCHY:
%    R0 naive (1/T)         : v_unamb = lambda/(4T)
%    R1 static schedule     : v_deal  = lambda*(C-1)/(4S)      [grating limit]
%    R2 time-varying sched  : full range, alias floor ~ 1/sqrt(M*C)   [random]
%  Then verifies the time-varying estimator reaches the non-uniform CRB.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:}); PFS={'FAIL','PASS'};
L2=load('step2_out.mat'); P=L2.P; D=L2.D; G=L2.G; SC=L2.SC;
lam=P.lambda; T=P.T; L=P.L; M=P.M; kk=2*pi*(2/lam);
try, randn('seed',4); rand('seed',4); catch, end
S=max(G.ssbTime)-min(G.ssbTime); C=6; W=D.W;
LOG('\n============ STEP 9: FEASIBILITY + TIME-VARYING SCHEDULE ============\n');
LOG('\n[params] lambda=%.4f mm  T=%.0f ms  S(burst span)=%.2f ms  C=%d  M=%d\n',lam*1e3,T*1e3,S*1e3,C,M);

%% ---- R0/R1/R2 predicted limits --------------------------------------
v0 = lam/(4*T);
v1 = lam*(C-1)/(4*S);
dssb = 4*P.Tsym; v2cap = lam/(4*dssb);
LOG('\n[theory] R0 naive     v_unamb = lam/4T        = %.4f m/s\n',v0);
LOG('         R1 static     v_deal  = lam(C-1)/4S   = %.4f m/s\n',v1);
LOG('         R2 time-vary  floor   = 1/sqrt(M*C)   = %.4f  (cap ~lam/4dssb=%.1f m/s)\n',1/sqrt(M*C),v2cap);

%% ---- build the three sampling sets for one target sector ------------
nb=log2(L); br=zeros(1,L); for i=0:L-1, br(i+1)=bin2dec(fliplr(dec2bin(i,nb)))+1; end
angT=G.ssbTime(br); i0=30; oStat=sort(angT(i0:i0+C-1)); mm=(0:M-1).';
tStat=reshape(mm*T+oStat,[],1);
tVar=[]; for m=0:M-1, tVar=[tVar; m*T+G.ssbTime(randperm(L,C)).']; end
tVar=tVar(:);

%% ---- VERIFY R1: static grating at lam(C-1)/2S -----------------------
vv=linspace(0,12,6000).';
chiIntra=abs(mean(exp(1j*kk*(-vv).*reshape(oStat,1,[])),2));
lm=find(chiIntra(2:end-1)>chiIntra(1:end-2)&chiIntra(2:end-1)>chiIntra(3:end))+1;
[~,om]=sort(chiIntra(lm),'descend'); vg1=vv(lm(om(1)));
r1=vg1/(2*v1);
LOG('\n[verify R1] static grating v=%.2f m/s ; scaling law lam(C-1)/2S=%.2f (ratio %.2f) [%s]\n',...
    vg1,2*v1,r1,PFS{1+(r1>0.5 && r1<2)});
LOG('            (exact: v_deal = lam/4*delta_max, delta_max = largest intra-sector gap;\n');
LOG('             uniform-spacing approximation gives lam(C-1)/4S -> scaling in (C-1)/S.)\n');

%% ---- VERIFY R2: time-varying alias floor ~ 1/sqrt(MC) ---------------
vgw=linspace(-70,70,20001).'; vt0=0;
chiVar=abs(mean(exp(1j*kk*(vt0-vgw).*reshape(tVar,1,[])),2));
sl=chiVar(abs(vgw-vt0)>1); side=max(sl); typ=median(sl);
LOG('[verify R2] time-varying sidelobe: typical=%.3f (~1/sqrt(MC)=%.3f), worst=%.3f\n',typ,1/sqrt(M*C),side);
LOG('            worst << unit true peak -> de-aliasing reliable [%s]\n',PFS{1+(side<0.7)});

%% ---- Time-varying estimator vs CRB & P(correct) ---------------------
Ns=numel(tVar); span=max(tVar)-min(tVar); dv=lam/(8*span); vg=(-70:dv:70).';
Phi=exp(1j*kk*tVar*vg.'); Stt=sum((tVar-mean(tVar)).^2); crb_v=@(eta)(lam/2)*sqrt(1./(2*eta*(2*pi)^2*Stt));
SNRax=-5:5:25; Ntr=80; rmse=zeros(size(SNRax)); pc=zeros(size(SNRax));
LOG('\n[time-varying] estimator vs CRB (principled: Nyquist grid + ML refine)\n');
for is=1:numel(SNRax)
  eta=10^(SNRax(is)/10); sg=sqrt(1/(2*eta)); es=[]; nc=0;
  for it=1:Ntr
    vt=-60+120*rand; y=exp(1j*2*pi*rand)*exp(1j*kk*vt*tVar)+sg*(randn(Ns,1)+1j*randn(Ns,1));
    c=abs(Phi'*y); [~,gi]=max(c); loc=vg(gi)+(-dv:dv/50:dv);
    [~,li]=max(abs(exp(1j*kk*tVar*loc)'*y)); e=loc(li)-vt;
    if abs(e)<W/2, nc=nc+1; es(end+1)=e; end
  end
  rmse(is)=sqrt(mean(es.^2)); pc(is)=nc/Ntr;
  LOG('   SNR=%+3d : P(correct)=%.3f  RMSE=%.3g  CRB=%.3g\n',SNRax(is),pc(is),rmse(is),crb_v(eta));
end
LOG('   VERIFY P(correct)@+15dB=%.2f [%s]  RMSE/CRB@+15dB=%.2f [%s]\n',...
   pc(SNRax==15),PFS{1+(pc(SNRax==15)>0.95)},rmse(SNRax==15)/crb_v(10^1.5),PFS{1+(rmse(SNRax==15)/crb_v(10^1.5)<2)});

%% ---- FEASIBILITY MAP + result figures ------------------------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');
Cr=2:64; vdeal=lam*(Cr-1)/(4*S);
f1=figure('position',[0 0 660 440]);
fill([2 64 64 2],[5 5 30 30],[.95 .9 .7],'edgecolor','none'); hold on;
plot(Cr,vdeal,'-','color',[.85 .33 .1],'linewidth',2);
plot([2 64],[v2cap v2cap],'-','color',[.1 .5 .8],'linewidth',2);
plot(C,lam*(C-1)/(4*S),'ko','markersize',8,'markerfacecolor','k');
text(8,3.8,'realistic C\approx6','fontsize',10);
set(gca,'yscale','log'); grid on; xlabel('beams illuminating a target,  C'); ylabel('de-aliasable velocity [m/s]');
legend('typical UAV speeds','R1 static: \lambda(C-1)/4S','R2 time-varying (full range)','location','east');
title('Fig1: feasibility map - static needs C\approx36 for UAVs; time-varying covers all C');
print(f1,'step9_fig1_feasibility.png','-dpng','-r110');

f2=figure('position',[0 0 660 400]); etaAx=10.^(SNRax/10);
semilogy(SNRax,rmse,'o-','color',[.1 .5 .8],'linewidth',1.6,'markerfacecolor',[.6 .8 .95]); hold on;
semilogy(SNRax,crb_v(etaAx),'k--','linewidth',1.3); grid on;
xlabel('SNR [dB]'); ylabel('velocity RMSE [m/s]'); legend('time-varying estimator','CRB','location','southwest');
title('Fig2: time-varying schedule reaches the CRB over full velocity range');
print(f2,'step9_fig2_rmse.png','-dpng','-r110');

f3=figure('position',[0 0 660 380]);
plot(vgw,20*log10(chiVar+1e-4),'-','color',[.1 .5 .8]); hold on;
plot([-70 70],20*log10([1 1]/sqrt(M*C)),'r--','linewidth',1.2); grid on; ylim([-40 2]);
xlabel('velocity [m/s]'); ylabel('|\chi(v)| [dB]'); legend('time-varying ambiguity','1/\surd(MC) floor');
title('Fig3: time-varying ambiguity - aliases become a 1/\surd(MC) floor');
print(f3,'step9_fig3_ambiguity.png','-dpng','-r110');

save('step9_out.mat','SNRax','rmse','pc','v0','v1','v2cap','Cr','vdeal','-v7');
LOG('\n[saved] step9_out.mat + 3 PNGs\n============ STEP 9 COMPLETE ============\n\n');

