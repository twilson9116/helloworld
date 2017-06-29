% Created by Tom White, 3/26/2017 
% tomw9116@gmail.com
function Results_Table=FMR_results(Table,plots_true_false)
%% v.1.0 Fit H_res and Linewidth to calculate properties of samples. 
% Calculates effective magnetization, gyromagnetic ratio, damping, and
% inhomogenous broadening by fitting resonance field and linewidth table.
% Input is the output from the function 'fit_spectrum.m' . Enter this
% either as the matlab variable directly or as the text file in the current
% directory enclosed with single quotes.


%% settings
if nargin<2;
    plots_true_false=false ;
end
 
%% FUNCTION CALCULATE PARAMETERS, PLOT, AND SAVE TO TABLE
if istable(Table);
else
    Table=readtable(Table);
end

Sample_list=unique(Table.Sample);
Results=zeros(length(Sample_list),8);
for i=1:length(Sample_list);
index=ismember(Table.Sample,Sample_list(i));

if nnz(index)>3
f=Table.frequency(index);
Hr=Table.H_res(index);
dH=Table.deltaH(index);

[kfit, kgof]=kittel_fit(Hr,f);
[dHfit,dHgof]=dH_fit(dH,f);

kfit_bounds=confint(kfit,0.68);
dHfit_bounds=confint(dHfit,0.68);

dHslope=dHfit.p1;
alpha=kfit.g*dHslope/2;

M_stdv=kfit.M-kfit_bounds(1,1); 
g_stdv=kfit.g-kfit_bounds(1,2);
dH_0_stdv=dHfit.p2-dHfit_bounds(1,1);
slope_stdv=dHslope-dHfit_bounds(1,2);
alpha_stdv=alpha*.5*sqrt((g_stdv/kfit.g)^2+(slope_stdv/dHfit.p1)^2);

Results(i,:)=[kfit.M,M_stdv,kfit.g,g_stdv,dHfit.p2,dH_0_stdv,alpha,alpha_stdv];

%% PLOTS
if plots_true_false
fighand=figure('Name', ['Kittel for ' Sample_list{i}]) ;
plot(kfit,Hr,f);
grid on; legend off; xlabel( 'Hr' ); ylabel( 'frequency' );

fighand=figure('Name',['Linewidth vs. Frequency for' Sample_list{i}]);
hold on
plot(dHfit,f,dH);
prediction = predint(dHfit,f);
plot(f,prediction,'m--');
grid on; legend off; xlabel( 'f, (GHZ)' ); ylabel( '\DeltaH' );
hold off
end %ends PLOTS
else 
    sprintf('Too few sweeps to fit Kittel for %s.',Sample_list{i})
end %ends if statement requiring enough data points to do kittel or dH fits
end %ends for loop i
Results_Table=[Sample_list array2table(Results)];
Results_Table.Properties.VariableNames = {'Sample' 'M_eff' 'M_sigma' 'gamma' 'g_sigma' 'dH_0' 'dH_sigma' 'alpha' 'a_sigma'};
writetable(Results_Table,strcat('Results_FMR_',datestr(now,'mmdd_HH_MM'))); 
end %ends the entire script
%% FIT: FREQUENCY VS. H_RES
function [kfit,kgof]=kittel_fit(Hres,freq);
[Hr2,f2]=prepareCurveData(Hres,freq);
kittel=fittype( 'g*sqrt(x*(x+M))', 'independent', 'x', 'dependent', 'y');
%order of outputs (M,g,x)
opts2 = fitoptions( kittel );
opts2.Display = 'Off';
opts2.Lower = [0 -Inf];
opts2.StartPoint = [0.0975 0.278];
[kfit,kgof]=fit(Hr2,f2,kittel,opts2);
end  %ends kittel function

%% FIT: DELTA H VS. FREQUENCY
    function [dH_fit,dHgof]=dH_fit(dH,freq);
[dH3,f3]=prepareCurveData(dH,freq);
[dH_fit,dHgof]=fit(f3,dH3,'poly1');
    end %ends lindwidth function