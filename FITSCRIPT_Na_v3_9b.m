%%
%%  Global fit: 278, 283, 295K 
%%             T1, T2, Diffusion data (only at 293K)
%%
mu0=4*pi*1e-7;
hbar  = 1.05457e-34;
Na=6.022e23;
const.hbar=hbar;
gammaNa=70.761e6;
%%gamma1H=267.513e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read experiments    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load EXP_23Na_T1_hf.dat
load EXP_23Na_T2_hf.dat
load EXP_23Na_Diff.dat
const.omega0=EXP_23Na_T1_hf(1,1)*gammaNa;
R1_hf278 = 1e3./EXP_23Na_T1_hf(1,[2:6:25]);
R2_hf278 = 1e3./EXP_23Na_T2_hf(1,[2:6:25]);
R1_hf283 = 1e3./EXP_23Na_T1_hf(1,[4:6:25]);
R2_hf283 = 1e3./EXP_23Na_T2_hf(1,[4:6:25]);
R1_hf295 = 1e3./EXP_23Na_T1_hf(1,[6:6:25]);
R2_hf295 = 1e3./EXP_23Na_T2_hf(1,[6:6:25]);
Diff_295 = EXP_23Na_Diff(1,[6:6:25]);
%%
const.R1b=[R1_hf278(1,4),R1_hf283(1,4),R1_hf295(1,4)]';
const.R2b=[R2_hf278(1,4),R2_hf283(1,4),R2_hf295(1,4)]';
const.Db =[0.0,0.0,Diff_295(1,4)]';

%%
errT1_278 = EXP_23Na_T1_hf(1,[3:6:25]);
errT2_278 = EXP_23Na_T2_hf(1,[3:6:25]);
errT1_283 = EXP_23Na_T1_hf(1,[5:6:25]);
errT2_283 = EXP_23Na_T2_hf(1,[5:6:25]);
errT1_295 = EXP_23Na_T1_hf(1,[7:6:25]);
errT2_295 = EXP_23Na_T2_hf(1,[7:6:25]);
errDiff_295 = EXP_23Na_Diff(1,[7:6:25]);
%%
errR1_278=( 2e3./( EXP_23Na_T1_hf(1,[2:6:25]) - errT1_278) - 2e3./( EXP_23Na_T1_hf(1,[2:6:25]) + errT1_278) );
errR2_278=( 2e3./( EXP_23Na_T2_hf(1,[2:6:25]) - errT2_278) - 2e3./( EXP_23Na_T2_hf(1,[2:6:25]) + errT2_278) );
errR1_283=( 2e3./( EXP_23Na_T1_hf(1,[4:6:25]) - errT1_283) - 2e3./( EXP_23Na_T1_hf(1,[4:6:25]) + errT1_283) );
errR2_283=( 2e3./( EXP_23Na_T2_hf(1,[4:6:25]) - errT2_283) - 2e3./( EXP_23Na_T2_hf(1,[4:6:25]) + errT2_283) );
errR1_295=( 2e3./( EXP_23Na_T1_hf(1,[6:6:25]) - errT1_295) - 2e3./( EXP_23Na_T1_hf(1,[6:6:25]) + errT1_295) );
errR2_295=( 2e3./( EXP_23Na_T2_hf(1,[6:6:25]) - errT2_295) - 2e3./( EXP_23Na_T2_hf(1,[6:6:25]) + errT2_295) );
%%
const.Conc=[50,300,700,10];
kk=1:3;
%% ydata(NN x NNtypes) , NN-number of fields and temps, number of rates (R1,R2) and concentrations
ydata(1,:)=[R1_hf278(kk),R2_hf278(kk),0.0,0.0,0.0];
ydata(2,:)=[R1_hf283(kk),R2_hf283(kk),0.0,0.0,0.0];
ydata(3,:)=[R1_hf295(kk),R2_hf295(kk),Diff_295(kk)];
yerr(1,:)=[errR1_278(kk),errR2_278(kk),1.0,1.0,1.0];
yerr(2,:)=[errR1_283(kk),errR2_283(kk),1.0,1.0,1.0];
yerr(3,:)=[errR1_295(kk),errR2_295(kk),errDiff_295(kk)];
const.NN=3; %% #temps
const.NNtypes=3*3; %% konc x rates (r1&r2)
const.m_cmc  = 100.0;   %% mM
const.tauc_s = [1.6436e-09, 1.3302e-09,1.1459e-09];  %% 1.6436e-09,(275K)
const.T=[275,283,295];
%%const.Dmic = [0.0, 0.0, 0.78e-10]; %% micelle MD diffusion (R=27.8 Angstrom) T=293K
%%const.Dcl5 = [0.0, 0.0, 3e-10]; %% pentamer MD difffusion at 293K
%%% surfactant proton study used obstruction factor for diffusion constant, using same here
Mw=172.3; %% surfactant molar weigh
xx_mic(1:3) =  (const.Conc(1:3) - const.m_cmc)./const.Conc(1:3); xx_mic(1)=0;
Phi=(1e-6*Mw.*const.Conc(1:3)/(0.9)).*xx_mic(1:3);
Obstr=1.0./(1.0 + 2.0.*Phi);
const.kB=1.38064852e-23;
R=27.4;
%% Viscosity of D2O at +25C and 0.1 MPa (http://www.iapws.org/relguide/TransD2O-2007.pdf )
const.visc = 1.250e-3; % 
Dmic=Obstr*(const.kB*const.T(3)/(6*pi*const.visc*R*1e-10)); Dmic(1)=0;
const.Dmic=zeros(3,3); const.Dmic(:,3)=Dmic';
R5=6.3;
Obstr5=1.0./(1.0 + 0.5.*Phi);
Dcl5=Obstr5*(const.kB*const.T(3)/(6*pi*const.visc*R5*1e-10));
const.Dcl5=zeros(3,3); const.Dcl5(:,3)=Dcl5';
%%
%%const.tauc_s = [3.1e-9,2.57e-9];  %% 1.6436e-09,(275K)
%%const.tauc_f = [1.5*295/283,0.95*293/295].*0.12e-9;  %% s
%%const.tauc_cl5=[295/283, 293/295].*0.285e-9; %% s
const.tauc_cl5=[4.2,3.4,2.85].*1e-10; %% (seconds) Md@283K & 295K sim D=k*10/4, tauD=R^2/(6*D)
%const.tauc_c15=[2.77,2.93,1.80].*1e-10; %% sec Md@275K 283 & 295K sim D=k*10/6, tauD=R^2/(6*D)
const.RR=8.31446261815324; %% J/(K⋅mol)
%% const.tauc_cl5=[295/283,1.0].*0.15e-9; %% s
const.ChiQ=102.0e3; %%kHz, Quist, Blom, Halle JOURNAL OF MAGNETIC RESONANCE vol 100, 267-281 (1992)
MCFLAG=1;
MCSTART=1;
savefile='MCMC_LM_1_v9b_quat_1.mat';
%% Mole fraction free lipids
%%Xc=0.1.*ones(1,3); %%
rng('shuffle');
beta=[1.0];
%% A2
A2=3e-1;
pcl=0.9;
const.beta=beta(1);
Nloop=size(beta,2);
start = [0.05,0.8,0.99,0.88,const.ChiQ,2e10,30e3]; %% size 3 ,
opt=optimset('MaxIter',20000,'MaxFunEvals',25000,'TolX',1e-10,'TolFun',1e-10);
Estimates=fminsearch(@(x)Rfit_tot_LM_v9(x,ydata,const,yerr),start,opt);
x = Estimates
[errtot,err,ymodel] = Rfit_tot_LM_v9(x,ydata,const,yerr);
if(MCFLAG)
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MC loop                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nparm = length(start);
N_warmup = 1500;
N_MCMC   = 5000;
N_sub    = 200;
%%
%%
xmin(1,1)=1e-7;
xmin(1,2:4)=0.0;
xmin(1,5) = (102-3)*1e3; %% halle J. Magn. Reson. 100 267-281, 1992
xmin(1,6) = 5e12; %% tauf0^-1
xmin(1,7) = 10e3; %% J/(K⋅mol)
%%
xmax(1,1) = 1.0;
xmax(1,2:4) = 1.0;
xmax(1,5) = (102+15)*1e3; %% halle J. Magn. Reson. 100 267-281, 1992
xmax(1,6) = 1e14; %% s-1
xmax(1,7) = 40e3; %% J/(K⋅mol)
%%
for iimc=1:Nloop,
const.beta=beta(1,iimc);
fprintf('inv temp %g\n',const.beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find suitable steplength    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(MCSTART)
  X0warmup = zeros(N_warmup,Nparm);
  for ii=1:Nparm,
    X0warmup(N_warmup,ii) = xmin(1,ii) + rand(1,1)*(xmax(1,ii)-xmin(1,ii) );
  end;
  X0warmup(1,:) = X0warmup(N_warmup,:);
  Xout = X0warmup(1,:)
  %X0(1,:) = abs(Estimates(1,:));
  Xbest = Xout(1,:);
  errbest = err;
  h=zeros(Nparm,1); %% MC step-size
  h=(xmax-xmin)./1e4;
else
  Ntmp=N_MCMC;
  load(savefile);
  disp(savefile);
  %const.D0=2.3e-8;
  %const.ED=3.5e4;
  const.beta=beta;
  disp('continue from previous outputfile');
  Xout=abs(Xbest);
  Nt = length(X0);
  Xout=abs(X0(Nt,1:Nparm));
  x=abs(Xbest);
  N_MCMC=Ntmp;
  errbest=1e10;
  h=zeros(Nparm,1); %% MC step-size
  h=(xmax-xmin)./1e4;
  const
end;
for ii=1:Nparm,
   if( Xout(1,ii) < xmin(1,ii) ) 
       disp('lower bondary');
       Xout(1,ii)=xmin(1,ii); 
       disp([ii,Xout(1,ii),xmin(1,ii)]);return; 
   end;
   if( Xout(1,ii)>xmax(1,ii) ) 
       disp('upper bondary');
      Xout(1,ii)=xmax(1,ii); 
      disp([ii,Xout(1,ii),xmax(1,ii)]);return;
   end;
end;
accave = zeros(1,Nparm);
errbest=err;
MM=0;
while( ((sum(accave<ones(1,Nparm).*0.4) || sum(accave>ones(1,Nparm).*0.6))) && MM<50 )
  accave = zeros(1,Nparm);
  %Xout=x;
  for i_sub=1:N_warmup,
    %indNparm=ones(1,Nparm);
    [Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_LM_v9);
    accave = accave + acc;
    Xbest = newXbest;
  end;
  accave=accave./(N_warmup);
  for ll=1:Nparm, fprintf('%g  ',accave(1,ll)); end;  fprintf(' MCMC-step\n');
  for ll=1:Nparm,
     if(accave(1,ll)<0.4)
       %fprintf('decrease %g  %d\n',h(1,ll),ll);
       h(1,ll)=h(1,ll)-h(1,ll)*0.1;
     end;
     if(accave(1,ll)>0.6)
       %fprintf('increase %g  %d\n',h(1,ll),ll);
       h(1,ll)=h(1,ll)+h(1,ll)*0.1;
     end;
  end; %ll
  MM=MM+1;
end; % while kk
%%if(h(1,7)>xmax(1,7)/10.0) 
%%  h(1,7)=xmax(1,7)/100.0;
%%end;
fprintf('time step callibrated, acc %g\n',accave);
fprintf('numb calibration rounds: %d\n',MM-1);
fprintf('err best : %g\n',errbest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Warmup sequence            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = tic;
%% run warm up block several times
accave = zeros(1,Nparm);
n  = 1;
err_warmup = zeros(N_warmup,1);
err_warmup(1,1) = errtot;
%indNparm=ones(1,Nparm);
[Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_LM_v9);
%
while(n < N_warmup),
  for i_sub=1:N_sub,
    %indNparm=find(substeps(i_sub,:));
    [Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_LM_v9);
    accave = accave + acc;
    Xbest = newXbest;
  end;
  n = n + 1;
  err_warmup(n,1) = err;
  X0warmup(n,:) = Xout;
end; %while n<N_warmup
fprintf('accave: %5.3e\n',accave./(N_sub.*N_warmup));
fprintf(' Warmaup done\n');
accave = zeros(1,Nparm);
errbest=err;
MM=0;
while( ((sum(accave<ones(1,Nparm).*0.4) || sum(accave>ones(1,Nparm).*0.6))) && MM<50 )
  accave = zeros(1,Nparm);
  %Xout=x;
  for i_sub=1:N_warmup,
    indNparm=ones(1,Nparm);
    [Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_LM_v9);
    accave = accave + acc;
    Xbest = newXbest;
  end;
  accave=accave./(N_warmup);
  for ll=1:Nparm, fprintf('%g  ',accave(1,ll)); end;  fprintf(' MCMC-step\n');
  for ll=1:Nparm,
     if(accave(1,ll)<0.4)
       %fprintf('decrease %g  %d\n',h(1,ll),ll);
       h(1,ll)=h(1,ll)-h(1,ll)*0.1;
     end;
     if(accave(1,ll)>0.6)
       %fprintf('increase %g  %d\n',h(1,ll),ll);
       h(1,ll)=h(1,ll)+h(1,ll)*0.1;
     end;
  end; %ll
  MM=MM+1;
end; % while kk
%%if(h(1,7)>xmax(1,7)/10.0) 
%%  h(1,7)=xmax(1,7)/100.0;
%%end;
fprintf('time step callibrated, acc %g\n',accave);
fprintf('numb calibration rounds: %d\n',MM-1);
fprintf('err best : %g\n',errbest);
%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = tic;
X0 = zeros(N_MCMC,Nparm);
X0(1,:) = Xout(1,:);
err_MCMC = zeros(N_MCMC,1);
err_MCMC(1,1) = err;
n=1;
accave = zeros(1,Nparm);
while(n < N_MCMC),
  for i_sub=1:N_sub,
    %indNparm=find(substeps(i_sub,:));
    [Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_LM_v9);
    accave = accave + acc;
    Xbest = newXbest;
  end;
  n = n + 1;
  if(mod(n,1000)==0 )
    fprintf('MCMC save no: %d\n',n);
 end;
  err_MCMC(n,1) = err;
  X0(n,:) = Xout;
end; %while n<N_MCMC
%Xstat = X0;
%
accave=accave./(N_sub.*N_MCMC);  
fprintf('MCMC in %g (min)\n',toc(t0)/60.0);
fprintf('accave: %5.3e\n',accave);
%histogram2(Xstat(:,1),Xstat(:,2),'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability');
%colorbar;
%xlabel('d_1 (Ångström)','fontsize',20);
%ylabel('d_2 (Ångström)','fontsize',20);
%figure;
%histogram2(Xstat(:,4),Xstat(:,3),'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability');
%colorbar;
%xlabel('C','fontsize',20);
%ylabel('r_{eff} (Angstrom)','fontsize',20);
save(savefile,'h','X0','accave','err_MCMC','Xbest','errbest','N_MCMC','const','ydata','yerr');
end; %% MCFLAG
end;
