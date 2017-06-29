%% Calculate paramaters from FMR Data
%If dH
function FMR_tom(input_file,dH_multiplier,gof_setting)
clc;
close all;
fclose all;

gof_for_estimate=0.8; 
if nargin<3
gof_setting=0.999;     
end
gof_setting_R=gof_setting;

if nargin<2
dH_multiplier=1;
end
%warningID='curvefit:cfit:subsasgn:coeffsClearingConfBounds';
% warning('off',warningID); %this supressess the warning that is generated when we set the slope and shift to zero; line 
%% This scans the file and seperates the data colums into their respective vectors
fileID=fopen(input_file); %opens text file
%this while loop finds the start of the data
startline=0;
stop=0;
while stop~=1
    tline=fgetl(fileID);
    s1=tline;
    s2='[Data]';
    stop=strcmp(s1,s2);
    startline=startline+1;
end
fclose(fileID);
startline=startline+2;
%below: this creates a matrix with the raw data file to matrix
data_matrix=dlmread(input_file,'\t',startline,0);%(filename,type,rowoffset,columnoffset)
%assign each row to its respective value frequency , field, I data, Q..
row=sortrows(data_matrix);
frequency=row(:,1);
H=row(:,2);
%I=row(:,3);
%Q=row(:,4);
R=row(:,5); %use this to compare to results

%% SEPERATE FREQUENCIES: this places the raw data into corresponding arrays
%each loop will process a new frequency sweep
freq_array=unique(frequency);                                               %create array of the unique frequencies
parameter_matrix=zeros(length(freq_array),9);                           %pre-allocates matrix that will store the data from each frequency
%Below: this loop seperates each frequency sweep into an array, fits I & Q
% creates R data, then plots the data
for jj=1:length(freq_array);                                                %
    freq_index=find(freq_array(jj)==frequency);                             % creates an array of the indecies where a frequency begins and starts- so this code will not work if frequency data is not grouped I am proud of this one
    startindex=freq_index(1);
    lastindex=freq_index(end);
    current_frequency=freq_array(jj);
    current_H_data=H(startindex:lastindex);
  %  current_I_data=I(startindex:lastindex);
  %  current_Q_data=Q(startindex:lastindex);
    current_R_data=R(startindex:lastindex);
   
   [H_res_estimate,delta_H_estimate,k1_estimate,k2_estimate,R_fitresult_estimate,R_GOF_estimate]=fit_data(current_H_data,current_R_data,current_frequency,gof_for_estimate);
   Hr1=H_res_estimate;
   dH1=delta_H_estimate;
   H_up=Hr1+dH_multiplier*dH1;
   H_low=Hr1-dH_multiplier*dH1;
   Htmpup=abs(current_H_data-H_up);
   Htmplow=abs(current_H_data-H_low);
   [~, indexup]=min(Htmpup);
   [~ ,indexlow]=min(Htmplow);
   R_chopped=current_R_data(indexlow:indexup);
   H_chopped=current_H_data(indexlow:indexup);
 %Modified R fit
    [H_res,delta_H,k1_out,k2_out,fitresult,GOF_out]=fit_data(H_chopped,R_chopped,current_frequency,gof_setting_R);
    
    confidence95=confint(fitresult);
    sigma=([H_res delta_H]-confidence95(1,1:2));
    sigma=(sigma./(3.92/2));
    Hrsigma=sigma(1);
    dH_sigma=sigma(2);
    SNR_dH=delta_H/dH_sigma;
    %save results from each fit into a common matrix
    parameter_matrix(jj,:)=[current_frequency H_res Hrsigma delta_H dH_sigma SNR_dH GOF_out k1_out k2_out  ];
    
%     %Plot Data
%     fighand=figure('Name',[num2str(current_frequency) 'GHZ' ' ' input_file]) ;
%      hold on    
%      plot(current_H_data,current_R_data,'.')
%      plot(fitresult,H_chopped,R_chopped);
%     grid on; legend off; xlabel( 'H' ); ylabel( 'R' );
%   hold off
end
%% frequency vs H_resonance
[Hr_prepared,frequency_prepared]=prepareCurveData(parameter_matrix(:,2),parameter_matrix(:,1));
Hr_fit=fittype( 'g*sqrt((x+Hk)*(x+Hk+M))', 'independent', 'x', 'dependent', 'y');
%order of outputs (Hk,M,g,x)
opts2 = fitoptions( Hr_fit );
opts2.Display = 'Off';
opts2.Lower = [0 0 -Inf];
opts2.StartPoint = [0 0.0975404049994095 0.278498218867048];
opts2.Upper = [Inf Inf Inf];
[fitresult_freq_vs_Hr,gof_f_vs_Hr]=fit(Hr_prepared,frequency_prepared,Hr_fit,opts2);
%Assign values
M=fitresult_freq_vs_Hr.M;
gamma=fitresult_freq_vs_Hr.g;
Hk=fitresult_freq_vs_Hr.Hk;
confidence_f_v_Hr_coeffs=confint(fitresult_freq_vs_Hr);
M_high=confidence_f_v_Hr_coeffs(2,2);
M_low=confidence_f_v_Hr_coeffs(1,2);
g_high=confidence_f_v_Hr_coeffs(2,3);
g_low=confidence_f_v_Hr_coeffs(1,3);
Hk_high=confidence_f_v_Hr_coeffs(2,1);
Hk_low=confidence_f_v_Hr_coeffs(1,1);

%plot
% fighand=figure('Name', ['f vs Hr' ' ' input_file]) ;
% plot(fitresult_freq_vs_Hr,Hr_prepared,frequency_prepared);
% grid on; legend off; xlabel( 'Hr' ); ylabel( 'frequency' );
%% delta_H vs. frequency
[dH_prepared,frequency_prepared]=prepareCurveData(parameter_matrix(:,4),parameter_matrix(:,1));
dH_fit=fittype('m*x+b','independent','x','dependent','y');
opts3 = fitoptions( dH_fit );
opts3.Display = 'Off';
opts3.Lower = [-Inf -Inf];
opts3.StartPoint = [0 0];
opts3.Upper = [Inf Inf];
[fitresult_dH_vs_f,gof_dh_vs_f]=fit(frequency_prepared,dH_prepared,dH_fit,opts3);

%Assign Values for results matrix
dH_slope=fitresult_dH_vs_f.m;
dH_intercept=fitresult_dH_vs_f.b;
confidence_dH_coeffs=confint(fitresult_dH_vs_f);
dH_slope_high=confidence_dH_coeffs(2,2);
dH_slope_low=confidence_dH_coeffs(1,2);
dH_intercept_high=confidence_dH_coeffs(2,1);
dH_intercept_low=confidence_dH_coeffs(1,1);
alpha=gamma*dH_slope/(2);
alpha_low=g_low*dH_slope_low/2;                                               
alpha_high=g_high*dH_slope_high/2;
%plot
fighand=figure('Name',['Inhomogeneous Broadening']);
hold on
plot(fitresult_dH_vs_f,frequency_prepared,dH_prepared);
prediction = predint(fitresult_dH_vs_f,frequency_prepared);
plot(frequency_prepared,prediction,'m--');
grid on; legend off; xlabel( 'f, (GHZ)' ); ylabel( '\DeltaH' );
hold off

%% WRITE TO OUTPUT TEXT FILE
%make folder
time=clock;
timestamp=horzcat(num2str(time(2)),num2str(time(3)),num2str(time(4)),num2str(time(5)));
output_folder=horzcat('Out',input_file,'_','dHm=',num2str(dH_multiplier),'_',timestamp); %name of output file
mkdir(output_folder);

%Create file
output_file=horzcat('RESULT','_',input_file,'_','dHmult=',num2str(dH_multiplier),'_',timestamp,'.dat'); % <- name for output data file
fileID=fopen([output_file],'w');
results=[M gamma alpha dH_intercept Hk];
results_low=[M_low g_low alpha_low dH_intercept_low Hk_low];
results_high=[M_high g_high alpha_high dH_intercept_high Hk_high];
%fprintf(fileID','%14s\r\n','Top row= results, 2nd row=low, 3rd row=high')
fprintf(fileID,'%12s %12s %12s %12s %12s\r\n','M(Oe)', '\Gamma(GHZ/Oe)', 'alpha()', 'dH_0(Oe)','H_k(Oe)');
fprintf(fileID,'%12f %12f %12f %12f %12f\r\n',results);
fprintf(fileID,'%12f %12f %12f %12f %12f\r\n',results_low);
fprintf(fileID,'%12f %12f %12f %12f %12f\r\n',results_high);
fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s\r\n','f(GHZ)','Hr(Oe)','Hr_sigma','dH(Oe)','dH_sigma','SNR_dH','rsquare','k1','k2');
[rr,~]=size(parameter_matrix);
%current_frequency H_res_R Hrsigma delta_H dH_sigma k1_out k2_out GOF_out SNR_dH
for qq=1:rr;
    fprintf(fileID,'%12f %12f %12f %12f %12f %12f %12f %12f %12f\r\n',parameter_matrix(qq,:));
end
fclose(fileID);
%move file to folder
%movefile(output_file,output_folder);
% copyfile(output_file,output_folder);
% mkdir('results')
% movefile(output_file,'results')
% movefile(input_file,output_folder);
% %copyfile(input_file,output_folder);
% save_all_figs(output_folder);

%% %% f vs Hr: FIT FUNCTION:this creates the fits to frequencies
    function [H_res,delta_H,k1_out,k2_out,fit_result,GOF]=fit_data(current_H_data,Amplitude_in,current_frequency,gof_setting)
        [H_prepared,Amplitude_prepared]=prepareCurveData(current_H_data,Amplitude_in);
        %initial estimates
        [~,index1]=min(Amplitude_prepared);
        [~,index2]=max(Amplitude_prepared);
        min_index=min(index1,index2);
        max_index=max(index1,index2);
        delta_H_initial=(H_prepared(max_index)-H_prepared(min_index))*sqrt(3)/2;
        H_res_initial=((H_prepared(max_index)-H_prepared(min_index))/2)+H_prepared(min_index);
        ep_initial=0;
        k1=1.5.*(H_res_initial./delta_H_initial).^3;
        k2=(H_res_initial./delta_H_initial).^2;
        
        fit_equation=fittype('k1*(4*(x-Hr))/(delta_H^2+4*(x-Hr)^2)^2+k2*((delta_H^2-4*(x-Hr)^2)/(delta_H^2+4*(x-Hr)^2)^2)','independent', 'x', 'dependent', 'y' );
        %fit options, order of variables is (H,deltaH,k1,k2,slope,offset)
        opts = fitoptions(fit_equation);
        opts.Display = 'Off';
        opts.Lower = [0 0 -Inf -Inf];
        opts.StartPoint = [H_res_initial delta_H_initial k1 k2];
        opts.Upper = [Inf Inf Inf Inf];
        [fit_result,gof]=fit(current_H_data,Amplitude_prepared,fit_equation,opts);
        
        %check Goodness of Fit
            nn=0;mm=0;
            if gof.rsquare<gof_setting;
            k1_fit = k1;
            k2_fit = k2.*[10 100 1e3 -1e3 -1e4];
            var=0;
            for ii = 1: length(k1_fit)
                for j=1:length(k2_fit)
                    opts.Lower = [0, 0, -Inf, -Inf];
                    opts.Upper = [inf, inf, Inf, Inf];
                    opts.StartPoint = [H_res_initial delta_H_initial k1_fit(ii) k2_fit(j)];
                    [fitresult1, gof1] = fit(current_H_data, Amplitude_prepared, fit_equation,opts);
                    nn=nn+1;
                    if gof1.rsquare>var
                        var=gof1.rsquare;
                        fit_result = fitresult1;
                        gof=gof1;
                        mm=mm+1;
                        rsquarecheck=gof.rsquare;
                    end
                    if gof.rsquare>=gof_setting
                        break
                    end
                end
                if gof.rsquare>=gof_setting
                    break
                end
                gof.rsquare;
            end%end for loop-iteration to improve GOF
            end
            nn;
            mm;
        H_res=fit_result.Hr;
        delta_H=fit_result.delta_H;
        k1_out=fit_result.k1;
        k2_out=fit_result.k2;
        GOF=gof.rsquare;
    end%end fit function
%%
%     function save_all_figs(save_location);
%         cd(save_location);
%         h=get(0,'children');
%         for ii=1:length(h);
%             saveas(h(ii),['figure',num2str(ii)],'fig');
%         end
%         cd('../');
%     end
end%this ends the whole show