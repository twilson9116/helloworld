
gamma=0.00294;
freq=3:18;
alpha=zeros(1,11);
dH_intercept=zeros(1,11);
alpha_error=zeros(1,11)

for i=1:11
[dH_prepared,frequency_prepared]=prepareCurveData(dH_Cu(i+1,:),freq);
dH_fit=fittype('m*x+b','independent','x','dependent','y');
% opts3 = fitoptions( dH_fit );
% opts3.Display = 'Off';
% opts3.Lower = [-Inf -Inf];
% opts3.StartPoint = [0 0];
% opts3.Upper = [Inf Inf];
[fitresult_dH_vs_f,gof_dh_vs_f]=fit(frequency_prepared,dH_prepared,dH_fit);%,opts3);

%Assign Values for results matrix
dH_slope=fitresult_dH_vs_f.m;
dH_intercept(i)=fitresult_dH_vs_f.b;
confidence_dH_coeffs=confint(fitresult_dH_vs_f);
dH_slope_high=confidence_dH_coeffs(2,2);
dH_intercept_high=confidence_dH_coeffs(2,1);
alpha(i)=gamma*dH_slope/(2);
alpha_high=gamma*dH_slope_high/2;
alpha_error(i)=(alpha_high-alpha(i))/1.96;

%plot
fighand=figure('Name',['Inhomogeneous Broadening']);
hold on
plot(fitresult_dH_vs_f,frequency_prepared,dH_prepared);
prediction = predint(fitresult_dH_vs_f,frequency_prepared);
plot(frequency_prepared,prediction,'m--');
grid on; legend off; xlabel( 'f, (GHZ)' ); ylabel( '\DeltaH' );
hold off
end