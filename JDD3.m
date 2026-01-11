% Par Hakansson (2016/11/05, Oulu university)
% call: [J0,J1w0,J2w0] = JDD3(tauc,omega0)
%
function [J0,J1w0,J2w0] = JDD3(tauc,omega0)
  J0 = tauc;
  n=1;
  J1w0 = tauc/(1.0+(n*omega0*tauc)^2);
  n=2;
  J2w0 = tauc/(1.0+(n*omega0*tauc)^2);
end

