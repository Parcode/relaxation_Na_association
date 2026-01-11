
%%  R1_R2model_LM_v2(xx_cl,xx_mic,ydata,const,jj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ymodel,R1interA,R2interA,R1interB,R2interB] = R1_R2model_LM_v9(A2,pcl,ChiQ,tauf0,xx_mic,ydata,const,jj);
  R2 =zeros(1,const.NNtypes/3);
  R1 =zeros(1,const.NNtypes/3);
  R1f=zeros(1,const.NNtypes/3);
  R2f=zeros(1,const.NNtypes/3);
  R1s=zeros(1,const.NNtypes/3);
  R2s=zeros(1,const.NNtypes/3);
  Ds=zeros(1,const.NNtypes/3);
  Df=zeros(1,const.NNtypes/3);
  Diff=zeros(1,const.NNtypes/3);
%%  tt=[295/283,1.0].*tauf;
  tauc_s = const.tauc_s(jj);
%%  tauc_f = const.tauc_f(jj);
  kk=tauf0(1).*exp(-tauf0(2)./(const.RR.*const.T)); %% s-1 , Arrenius
  tau_fast = 1./kk;
  tauc_f = tau_fast(jj);
  %tauc_f=const.tauc_s(jj);
%%  tauc_f=tt(jj);
  tauc_cl5=const.tauc_cl5(jj);
  omega0=const.omega0;
  ResA=A2(1);
  ResA=0.0;
%%  Q = const.ChiQ^2*(3*pi^2/4.0);
  Q = ChiQ^2*(3*pi^2/4.0);
  Q = ChiQ^2*(pi^2/5.0);
  [J0s,J1s,J2s]       = JDD3(tauc_s,omega0);
  [J0f,J1f,J2f]       = JDD3(tauc_f,omega0);
  [J0cl5,J1cl5,J2cl5] = JDD3(tauc_cl5,omega0);
%%
  R1f  = 2.0*Q*(1-ResA)*(J1f + 4.0*J2f);
  R2f  = Q*(1-ResA)*(3.0*J0f + 5.0*J1f + 2.0*J2f);
  R1s  = 2.0*Q*( ResA )*(J1s + 4.0*J2s);
  R2s  = Q*( ResA )*(3.0*J0s + 5.0*J1s + 2.0*J2s);
  R1cl5= 2.0*Q*(J1cl5 + 4.0*J2cl5);
  R2cl5= Q*(3.0*J0cl5 + 5.0*J1cl5 + 2.0*J2cl5);
%%  R1f5  = 2.0*Q*(1-A2(2))*(J1f + 4.0*J2f);
%%  R2f5  = Q*(1-A2(2))*(3.0*J0f + 5.0*J1f + 2.0*J2f);
  R1b  = const.R1b;
  R2b  = const.R2b;
  Db   = const.Db;
  Dmic = const.Dmic;
  Dcl5 = const.Dcl5;
  pcl_conc=[pcl(1),pcl(2),pcl(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total contributons: %%
%%%%%%%%%%%%%%%%%%%%%%%%%
  %%ii=1;
  %%R2(1,ii) = abs(xx_cl(jj))*(const.R2mono(ii,1) + R2interA) + xx_mic(jj)*(const.R2cl(ii,1) + R2interB)
  %%const.R2cl(ii,1)
  %%R2interB
  for ii=1:const.NNtypes/3,
%%  There are free Na, known from experiment Rib
%%  and fraction following "cl" cluster
%%  The cluster consists of micelle part and small pentamer part
    R1cl     = xx_mic(ii)*(R1s+R1f) + (1.0-xx_mic(ii))*(R1cl5);
    R2cl     = xx_mic(ii)*(R2s+R2f) + (1.0-xx_mic(ii))*(R2cl5);
    Dcl      =xx_mic(ii)*Dmic(ii,jj) + (1.0-xx_mic(ii))*Dcl5(ii,jj);
    R1(1,ii) = pcl_conc(ii)*R1cl + (1.0 - pcl_conc(ii))*R1b(jj);
    R2(1,ii) = pcl_conc(ii)*R2cl + (1.0 - pcl_conc(ii))*R2b(jj);
    Diff(1,ii)=pcl_conc(ii)*Dcl  +(1.0 - pcl_conc(ii))*Db(jj);
  end;
  nn=1;
  ymodel(1,:)=[R1,R2,Diff];
end %% R1_R2model function
%%
