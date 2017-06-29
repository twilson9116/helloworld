% % 
%     p=zeros(length(Ir_Cu),9);
%     p(:,1)=Ir_Cu;
%    for i=1:length(Ir_Cu)    
%     p(i,2:end)=process_Hres_dH(Hres_Cu(i,:),dH_Cu(i,:));
%         
%     end 

function parameters=process_Hres_dH(Hr_v,dH_v)
%%take a vectors of Hres and deltaH for 2:18 GHZ and output params
f1=3:18;

%[M, M_sigma, g, g_sigma, dH0, dH0_sigma, alpha, alpha_sigma];
%first M and g by fit Hres vs freq
[Hr,f]=prepareCurveData(Hr_v,f1);
Hr_fit=fittype( 'g*sqrt(x*(x+M))', 'independent', 'x', 'dependent', 'y');
opts2 = fitoptions( Hr_fit );
opts2.Display = 'Off';
opts2.Lower = [0 -Inf];
opts2.StartPoint = [0.09755 0.278];
opts2.Upper = [Inf Inf];
[fitted_Hr,~]=fit(Hr,f,Hr_fit,opts2);


M=fitted_Hr.M;
g=fitted_Hr.g;
confidence_f_v_Hr_coeffs=confint(fitted_Hr);
M_high=confidence_f_v_Hr_coeffs(2,1);
g_high=confidence_f_v_Hr_coeffs(2,2);
M_sigma=(M_high-M)/(1.96);
g_sigma=(g_high-g)/(1.96);


[dH, f]=prepareCurveData(dH_v,f1);

dH_fit=fittype('m*x+b','independent','x','dependent','y');
opts3 = fitoptions( dH_fit );
opts3.Display = 'Off';
opts3.Lower = [-Inf -Inf];
opts3.StartPoint = [0 0];
opts3.Upper = [Inf Inf];
[fitted_dH,~]=fit(f,dH,dH_fit,opts3);
dH_slope=fitted_dH.m;
dH0=fitted_dH.b;
conf_dH=confint(fitted_dH);
dH_slope_high=conf_dH(2,2);
dH_intercept_high=conf_dH(2,1);
alpha=g*dH_slope/(2);                                            
alpha_sigma= alpha*sqrt((g_sigma/g)^2+(((dH_slope_high-dH_slope)/(3.92/2))/dH_slope)^2);
dH0_sigma=(dH_intercept_high-dH0)/(3.92/2);


%alpha_high=g_high*dH_slope_high/2; 
parameters=[M, M_sigma, g, g_sigma, dH0, dH0_sigma, alpha, alpha_sigma];

end






