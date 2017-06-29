

[r,~]=size(Hr_M); param=zeros(r,4);
f1=4:18;

for i=1:r;
    Hr_v=Hr_M(i,:);
  
[Hr,f]=prepareCurveData(Hr_v,f1);
Hr_fit=fittype( 'g*sqrt(x*(x+M))', 'independent', 'x', 'dependent', 'y');
opts2 = fitoptions( Hr_fit );
opts2.Display = 'Off';
opts2.Lower = [0 -Inf];
opts2.StartPoint = [0.09755 0.278];
opts2.Upper = [Inf Inf];
[fitted_Hr,~]=fit(Hr,f,Hr_fit,opts2);


M=fitted_Hr.M
g=fitted_Hr.g
confidence_f_v_Hr_coeffs=confint(fitted_Hr);
M_high=confidence_f_v_Hr_coeffs(2,1);
g_high=confidence_f_v_Hr_coeffs(2,2);
M_sigma=(M_high-M)/(1.96);
g_sigma=(g_high-g)/(1.96);

param_M(i,:)=[M M_sigma g g_sigma]
end