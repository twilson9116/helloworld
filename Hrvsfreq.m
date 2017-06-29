
gamma=0.00294;
freq=3:18;
M=zeros(1,11);
M_error=zeros(1,11)

%% frequency vs H_resonance
for i=1:11;

[Hr_prepared,frequency_prepared]=prepareCurveData(Hr_Cu(i+1,:),freq);
Hr_fit=fittype( '0.00294*sqrt(x*(x+M))', 'independent', 'x', 'dependent', 'y');
%order of outputs (M,g,x)
% opts2 = fitoptions( Hr_fit );
% opts2.Display = 'Off';
% opts2.Lower = [0 -Inf];
% opts2.StartPoint = [0.0975 0.278];
% opts2.Upper = [Inf Inf];
[fitresult_freq_vs_Hr,gof_f_vs_Hr]=fit(Hr_prepared,frequency_prepared,Hr_fit);
%Assign values
M(i)=fitresult_freq_vs_Hr.M;
confidence_f_v_Hr_coeffs=confint(fitresult_freq_vs_Hr);
M_high=confidence_f_v_Hr_coeffs(2,1);
M_error(i)=(M_high-M(i))/1.96;
%plot
fighand=figure('Name', ['f vs Hr' ' ']) ;
plot(fitresult_freq_vs_Hr,Hr_prepared,frequency_prepared);
grid on; legend off; xlabel( 'Hr' ); ylabel( 'frequency' );
end