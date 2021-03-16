nu=50; 
nd = 10;;
%A = randn(m,n);
%b = randn(m,1);
Wthd=W_threshold_d';
Wthu=W_threshold_u';
Pum=200*ones(nu,1);
Pdm=10000;
Wm=200;
phi=1.5;
N0=10.^(-17.3);
u=160;
gdstar=46;
gustar=46;
alphad=((alpha_users_linear).^-1)';
alphau=((alpha_sensors_linear).^-1)';
tao=0.01*10^(-3);
taor=tao.^-1;
srtaor=sqrt(tao).^-1;
qinv=qfuncinv(2e-8);
lu=zeros(nu,1);
ld=zeros(nd,1);
cvx_begin
variable wd(nd)
variable wu(nu)
  
  for i=1:nd
      x=taor*u*log(2)*EB(i)*+qinv*srtaor;
      e(i)=exp(x);
  end
  temp1=phi*N0.*(pow_p(wd,6)).*alphad.*(exp(taor*u*log(2).*EB.*(pow_p(wd,-1))+qinv*srtaor.*(pow_p(wd,-1)))-1);
  temp2=phi*N0.*(power(wu,6)).*alphau.*(exp(taor*u*log(2).*(power(wu,-1))+qinv*srtaor.*(pow_p(wu,-1)))-1);
minimize(max(temp2)+sum(temp1))
subject to 
temp2.*gustar-Pum<=0
sum((temp1.*gdstar))-Pdm<=0
sum(wd)+sum(wu)-Wm<=0
ld <= wd <= Wthd
 lu <= wu <= Wthu
cvx_end
