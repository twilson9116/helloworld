function Out=process_sweeps(input_file,gof_setting,dH_multiplier);
%% Processes sweeps and give matrix of H_res and delta_H in columns from 2:18GHZ each
% note that maxiter=8000 and there are external iterations for initial
% guesses
%clc;
close all;
fclose all;
fileID=fopen(input_file); %opens text file
disp(input_file)
Out=zeros(1,102);

gof_for_estimate=0.99;
if nargin<3
    gof_setting=0.99999;
end
GOF=gof_setting;

if nargin<2
    dH_multiplier=1;
end
% if nargout()=2
%Out2=zeros(1,51);
%end
startline=0;
stop=0;
while stop~=1
    tline=fgetl(fileID);
    s1=tline;
    s2='[Data]';
    stop=strcmp(s1,s2);
    startline=startline+1;
end%this while loop finds the start of the data
fclose(fileID);
startline=startline+2;

data_matrix=dlmread(input_file,'\t',startline,0);
row=sortrows(data_matrix);
frequency=row(:,1);
frequency=round(frequency,1);
H=row(:,2);
R=row(:,5);
H=-H;
freq_uniques=unique(frequency);

%disp(length(freq_uniques))

for jj=1:length(freq_uniques);                                                %
    %disp(jj)
    if any([2:18]==freq_uniques(jj))
        freq_index=find(freq_uniques(jj)==frequency);                             % creates an array of the indecies where a frequency begins and starts- so this code will not work if frequency data is not grouped I am proud of this one
        startindex=freq_index(1);
        lastindex=freq_index(end);
        current_frequency=freq_uniques(jj);
        current_H_data=H(startindex:lastindex);
        current_R_data=R(startindex:lastindex);
        
        %initial estimates
        
        [~,index_min]=min(current_R_data);
        [~,index_max]=max(current_R_data);
        low_index=min(index_min,index_max);
        high_index=max(index_min,index_max);
        Hr_1=((current_H_data(high_index)-current_H_data(low_index))/2)+current_H_data(low_index);
        dH_1=(sqrt(3)/2).*(current_H_data(high_index)-current_H_data(low_index));
        
        k1_guess1=1.5.*(Hr_1./dH_1).^3;
        k2_guess1=(Hr_1./dH_1).^2;
       
[initial_fit,initial_gof]=fit_data(current_H_data,current_R_data,gof_for_estimate,Hr_1,dH_1,k1_guess1,k2_guess1);
         
            Hr_2=initial_fit.Hr;
            dH_2=initial_fit.delta_H;
             k1_guess2=initial_fit.k1;
            k2_guess2=initial_fit.k2;
            
            %chop the data that will be fit
            H_up=Hr_2+dH_multiplier*dH_2; % magnitude of upper and lower fields that will be chopped
            H_low=Hr_2-dH_multiplier*dH_2;
            Upper_chop=abs(current_H_data-H_up);% subtract scalar field value of chop locations so that it equals minimum absoluite value
            Lower_chop=abs(current_H_data-H_low);
            [~, indexup]=min(Upper_chop);
            [~ ,indexlow]=min(Lower_chop);
            R_chopped=current_R_data(indexup:indexlow);
            H_chopped=current_H_data(indexup:indexlow);
  %actual fit is right here    
            try
            [fitresult,GOF_out]=fit_data(H_chopped,R_chopped,GOF,initial_fit.Hr,initial_fit.delta_H,initial_fit.k1,initial_fit.k2);
            catch 
             fitresult=initial_fit;
            end
            H_res=fitresult.Hr;
            delta_H=fitresult.delta_H;
          k1=fitresult.k1;
          k2=fitresult.k2;
       fighand=figure('Name',[input_file, 'freq=' num2str(freq_uniques(jj))]);     
        hold on
        plot(current_H_data,current_R_data,'o r')
       plot(H_chopped,R_chopped,'o b');
        plot(fitresult);
        hold off;
            
            conf_95=confint(fitresult);
            Hr_sigma=(conf_95(2,1)-H_res)/1.96;
            dH_sigma=(conf_95(2,2)-delta_H)/1.96;
            Out(freq_uniques(jj)-1)=H_res;
            Out(freq_uniques(jj)-1+17)=delta_H;
            Out(freq_uniques(jj)-1+34)=Hr_sigma;
            Out(freq_uniques(jj)-1+51)=dH_sigma;
            Out(freq_uniques(jj)-1+68)=k1;
            Out(freq_uniques(jj)-1+85)=k2;
            
            
            % if nargout()=2
            %         k1_out=fitresult.k1;
            %         k2_out=fitresult.k2;
            %Out2(jj)=k1_out;
            %Out2(jj+17)=k2_out;
            %Out2(jj+34)=GOF_out;
            %end
     
       end
end
fig_array=get(0,'children');
mkdir('figures with dH');
cd('figures with dH');
for i=1:length(fig_array);
savefig(fig_array(i),[fig_array(i).Name '.fig'],'compact')
end
cd('../')
close all
function [best_fit,best_gof]=fit_data(current_H_data,Amplitude_in,gof_setting,Hr_guess,dH_guess,k1_guess,k2_guess)
       
    [H_prepared,Amplitude_prepared]=prepareCurveData(current_H_data,Amplitude_in);

        fit_equation=fittype('k1*(delta_H*4*(x-Hr))/(delta_H^2+4*(x-Hr)^2)^2+k2*((delta_H^2-4*(x-Hr)^2)/(delta_H^2+4*(x-Hr)^2)^2)','independent', 'x', 'dependent', 'y' );
        %fit options, order of variables is (H,deltaH,k1,k2,slope,offset)
        opts = fitoptions(fit_equation);
        opts.Display = 'Off';
        opts.Lower = [0 0 -Inf -Inf];
        opts.StartPoint = [Hr_guess dH_guess k1_guess k2_guess];
        opts.Upper = [Inf Inf Inf Inf];
       opts.MaxIter=8000;
        [fit_try,gof_current]=fit(H_prepared,Amplitude_prepared,fit_equation,opts);
        
        %check Goodness of Fit
         best_gof=gof_current.rsquare;
         best_fit=fit_try;
         if best_gof<gof_setting;
           
            %pause 
             k1_array = k1_guess.*[1 -100 100];
             k2_array = k2_guess.*[100 1e3 -1e3];
             
             for ii = 1: length(k1_array)
         
                 for j=1:length(k2_array)
                   sprintf('Best GOF= %.3f \nCurrent GOF= %.3f',best_gof,gof_current.rsquare)
                     opts.StartPoint = [Hr_guess dH_guess k1_array(ii) k2_array(j)];
                     [fit_try, gof_current] = fit(H_prepared, Amplitude_prepared,fit_equation,opts);
                     if gof_current.rsquare>best_gof
                         best_gof=gof_current.rsquare;
                         best_fit = fit_try;
                     end                   
                     if best_gof>=gof_setting
                         break
                     end
                 end
                 if best_gof>=gof_setting
                     break
                 end
             end
         end

    end%end fit function
    
end




