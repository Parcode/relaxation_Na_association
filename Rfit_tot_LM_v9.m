%% call: [err] = Rfit_tot_B(x,ydata,const,yerr);
%% Input: x: array of fitting parameters
%%        ydata: experimental relaxation rates (DeltaR = 1/T2 - 1/T1 )
%%        const: fixed constants in the model (translational diffusion ...)
%%        yerr: standard deviation of ydata
%% Output: err = sqrt( 1/(Nexp - Nparm) *  sum_i (ydata_i - model_i)/yerr_i )
function [errtot,err,ymodel] = Rfit_tot_LM_v9(x,ydata,const,yerr)
%%
%%
%% start = [R,tau_ex,Sscale,m_cmc];
  Mw=172.3; %% g/mol
  Vtot=1.0; %
  err = zeros(const.NN,1);
  errtot=0.0;
  ymodel = zeros(const.NN,const.NNtypes);
  slope=0.0;
  xx_mic  =  (const.Conc(1:3) - (1.0+slope).*const.m_cmc)./const.Conc(1:3)+slope;
  xx_mic(1)=0.0;
%%  xx_mic  =  (const.Conc - (1.0).*const.m_cmc)./const.Conc;
  ii=find(xx_mic>1);
  if(length(ii)>0)
    xx_mic(ii)=1.0
  end;
%%
  for jj=1:const.NN, %% over number of temperatures (currently two) and field-strengths (currently one)
    [ymodel(jj,1:const.NNtypes)] = R1_R2model_LM_v9(x(1),x(2:4),x(5),x(6:7),xx_mic,ydata,const,jj);
    for kk=1:const.NNtypes,
      err(jj,1) = err(jj,1) + (ydata(jj,kk)-ymodel(jj,kk) )^2/yerr(jj,kk)^2;
    end;
  end;
  errtot = sum(err(:));
  %% errtot = 0.5*errtot/(Nobs-Nfit);
  errtot = (const.beta)*0.5*errtot;
end

