%% ========================================================================
%  STEP 2 | Scene definition & channel synthesis (+ verification)
%  Loads step1_out.mat, builds received cube Y[n,m,ell], unit-tests the
%  injected delay/Doppler, and plots range + range-angle to confirm the
%  channel BEFORE any velocity/aliasing processing.
%% ========================================================================
clear; clc; close all;
try, pkg load signal; catch, end
LOG=@(varargin) fprintf(varargin{:});
S1=load('step1_out.mat'); P=S1.P; D=S1.D; G=S1.G;
LOG('\n================ STEP 2: SCENE & CHANNEL SYNTHESIS ================\n');

%% ---- 1. Array & scene -------------------------------------------------
A.Na = 16;                       % gNB ULA elements (half-wavelength)
% targets: [R(m)  v(m/s)  theta(deg)  |alpha|]
SC.tab = [100  12  -30  1.0 ;      % T1  fast, aliases
          250  25   10  1.0 ;      % T2  co-cell with T3 (same R,theta)
          250   8   10  0.7 ;      % T3  co-cell: separable only in Doppler
          400   3   40  0.8 ];     % T4  slow
SC.K = size(SC.tab,1);
SC.R=SC.tab(:,1); SC.v=SC.tab(:,2); SC.th=SC.tab(:,3); SC.amp=SC.tab(:,4);
rng_seed=7; try, rand('seed',rng_seed); randn('seed',rng_seed); catch, end
SC.phase = 2*pi*rand(SC.K,1);
SC.alpha = SC.amp.*exp(1j*SC.phase);
SC.SNRdB = 20;                    % nominal per-RE SNR at beam peak

SC.tau = 2*SC.R/P.c;              % delays
SC.fd  = 2*SC.v/P.lambda;         % Doppler (true, unaliased)
SC.kfold = round(SC.v./D.W);      % alias fold index
SC.vmeas = SC.v - SC.kfold*D.W;   % aliased (folded) velocity

LOG('\n[1] Scene (K=%d targets)\n', SC.K);
LOG('    %-3s %6s %6s %6s | %8s %9s %6s %9s\n','id','R[m]','v','th','tau[us]','fd[Hz]','fold','vmeas');
for k=1:SC.K
  LOG('    T%-2d %6.0f %6.1f %6.0f | %8.3f %9.1f %6d %9.4f\n', ...
      k,SC.R(k),SC.v(k),SC.th(k),SC.tau(k)*1e6,SC.fd(k),SC.kfold(k),SC.vmeas(k));
end
LOG('    NOTE: all |v| >> vmax=%.3f m/s  -> every target aliases\n', D.vmax);
LOG('    NOTE: T2 & T3 share (R=250,theta=10) -> unresolvable except in Doppler\n');

%% ---- 2. Beam gain (ULA array factor, magnitude) ----------------------
beamGain=@(th_t,th0) abs( sin(A.Na*pi/2*(sind(th_t)-sind(th0))) ...
                        ./ (A.Na*sin(pi/2*(sind(th_t)-sind(th0))+eps)) );

%% ---- 3. Synthesize received cube  Y[n,m,ell] -------------------------
Nsc=P.Nsc; M=P.M; L=P.L; n=(0:Nsc-1).';
Y=zeros(Nsc,M,L);
for k=1:SC.K
  rp  = exp(-1j*2*pi*n*P.scs*SC.tau(k));         % Nsc x 1  (delay)
  dop = exp( 1j*2*pi*SC.fd(k)*G.tGrid);          % M x L    (Doppler)
  g   = beamGain(SC.th(k),G.thetaBeam);          % 1 x L    (illumination)
  db  = dop .* g;                                % M x L
  Y   = Y + SC.alpha(k)*reshape(rp,[Nsc 1 1]).*reshape(db,[1 M L]);
end
% AWGN at nominal SNR (referenced to unit beam-peak amplitude)
sigma = sqrt(10^(-SC.SNRdB/10)/2);
Y = Y + sigma*(randn(size(Y))+1j*randn(size(Y)));
LOG('\n[3] Cube Y synthesized: size %dx%dx%d (n x m x ell), sigma=%.4g\n',Nsc,M,L,sigma);

%% ---- 4. UNIT TESTS on a clean single-target cube ---------------------
LOG('\n[4] Unit tests (noise-free single-target)\n');
Nfft=1024; rax=(0:Nfft-1)/Nfft*D.Rmax;   % range axis
pass=true;
for k=[1 2]
  rp=exp(-1j*2*pi*n*P.scs*SC.tau(k)); dop=exp(1j*2*pi*SC.fd(k)*G.tGrid);
  g=beamGain(SC.th(k),G.thetaBeam);
  Yk=SC.alpha(k)*reshape(rp,[Nsc 1 1]).*reshape(dop.*g,[1 M L]);
  % (a) delay: range-FFT at the peak beam, m=1
  [~,lb]=min(abs(G.thetaBeam-SC.th(k)));
  prof=abs(ifft(Yk(:,1,lb),Nfft));
  [~,pk]=max(prof); Rhat=rax(pk);
  okR = abs(Rhat-SC.R(k)) <= D.dR;             % within one range bin
  % (b) Doppler: intra-burst phase slope between two SSBs in burst m=1
  cand=[lb-1 lb+1]; cand=cand(cand>=1 & cand<=L);
  [~,ix]=min(abs(G.tGrid(1,cand)-G.tGrid(1,lb)));  % temporally-closest SSB
  la=cand(ix); dt=G.tGrid(1,lb)-G.tGrid(1,la);
  fdhat=angle(Yk(1,1,lb)/Yk(1,1,la))/(2*pi*dt);
  vhat=P.lambda*fdhat/2;
  okD = abs(vhat-SC.v(k)) <= 0.5;              % m/s tol
  pass = pass && okR && okD;
  sR='FAIL'; if okR, sR='PASS'; end
  sD='FAIL'; if okD, sD='PASS'; end
  LOG('    T%d: Rhat=%7.2f m (true %3.0f) [%s]   vhat=%6.2f m/s (true %2.0f) [%s]\n', ...
       k,Rhat,SC.R(k),sR,vhat,SC.v(k),sD);
end
sAll='FAIL'; if pass, sAll='PASS'; end
LOG('    OVERALL UNIT TESTS: %s\n', sAll);

%% ---- 5. Range profile (verify all peaks land correctly) --------------
try, graphics_toolkit('gnuplot'); catch, end
set(0,'defaultfigurevisible','off');
% incoherent range profile summed over all beams & burst sets
RP=zeros(Nfft,1);
for mm=1:M, for ll=1:L, RP=RP+abs(ifft(Y(:,mm,ll),Nfft)); end, end
RP=RP/max(RP);
f1=figure('position',[0 0 760 320]);
plot(rax,20*log10(RP+1e-6),'b'); hold on; grid on;
for k=1:SC.K, xl=SC.R(k); plot([xl xl],[-60 2],'r--'); end
xlim([0 500]); ylim([-50 2]); xlabel('range [m]'); ylabel('norm. power [dB]');
title('Fig 1: range profile - dashed = true target ranges');
print(f1,'step2_fig1_range.png','-dpng','-r110');

%% ---- 6. Range-angle map (verify range AND angle; show co-cell overlap)
RAmap=zeros(Nfft,L);
for ll=1:L
  acc=zeros(Nfft,1);
  for mm=1:M, acc=acc+abs(ifft(Y(:,mm,ll),Nfft)); end
  RAmap(:,ll)=acc;
end
RAmap=RAmap/max(RAmap(:));
f2=figure('position',[0 0 720 420]);
imagesc(G.thetaBeam, rax, 20*log10(RAmap+1e-6)); axis xy; caxis([-30 0]);
xlabel('beam azimuth [deg]'); ylabel('range [m]'); ylim([0 500]); colorbar;
title('Fig 2: range-angle map (T2&T3 overlap at 250 m, 10 deg)');
hold on; for k=1:SC.K, plot(SC.th(k),SC.R(k),'w+','markersize',10,'linewidth',1.5); end
print(f2,'step2_fig2_rangeangle.png','-dpng','-r110');

%% ---- 7. Save ---------------------------------------------------------
save('step2_out.mat','P','D','G','A','SC','Y','-v7');
LOG('\n[7] Saved step2_out.mat + 2 PNGs\n');
LOG('================ STEP 2 COMPLETE ================\n\n');

