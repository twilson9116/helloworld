%% I am re-processing all Hres and dH data
%I want to do 2 things
% 1) remove all saturated points
% 2)save the plots so i can see the quality of the fits

%1) get raw data 

function Data_Matrix=fun924(inputfolder)

clc;
close all;
fclose all;
listing=dir(pwd)
listing=listing(3:end);

mkdir('output');

names_list=cell(length(listing),1);
filenumber=zeros(length(listing),1);
gof_for_estimate=0.5;
gof_setting=0.9;

Data_Matrix=zeros(length(listing),(17*6+1));
warning_matrix=zeros(length(listing),(17+1));
dH_multiplier=1.3;
for file_count=1:length(listing)
    filenumber(file_count)=file_count;
    names_list{file_count}=listing(file_count).name;
    input_file=listing(file_count).name;
disp(listing(file_count).name)
fileID=fopen(input_file);
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
data_matrix=dlmread(input_file,'\t',startline,0);
row=sortrows(data_matrix);

f=row(:,1);
H=row(:,2);
if mean(H)<0
    H=-H;
end
R=row(:,5);

% 2- Seperate into frequencies
%use if statement to delete rows that are x distance away from f integer
%values
f=f(abs(f-round(f))<.2);

f=round(f);
f=f(f<18.1);
freq_array=unique(f);  

for jj=1:length(freq_array);  
    freq_index=find(freq_array(jj)==f);                             % creates an array of the indecies where a frequency begins and starts- so this code will not work if frequency data is not grouped I am proud of this one
    startindex=freq_index(1);
    lastindex=freq_index(end);
    current_frequency=freq_array(jj);
    current_H_data=H(startindex:lastindex);
    current_R_data=R(startindex:lastindex);
if length(current_R_data)>20
   [H_res_estimate,delta_H_estimate,k1_estimate,k2_estimate,R_fitresult_estimate,GOF_estimate]=fit_data(current_H_data,current_R_data,current_frequency,gof_for_estimate);
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

        if length(R_chopped)>10
            [H_res,delta_H,k1_out,k2_out,fitresult,GOF_out]=fit_data(H_chopped,R_chopped,current_frequency,gof_setting);
        else 
            H_res=H_res_estimate; delta_H=delta_H_estimate; k1_out=k1_estimate; k2_out=k2_estimate; fitresult=R_fitresult_estimate; GOF_out=GOF_estimate;
            warning_matrix(file_count,current_frequency)=1;
        end 
    
confidence95=confint(fitresult);
    sigma=([H_res delta_H]-confidence95(1,1:2));
    sigma=(sigma./(3.92/2));
    Hrsigma=sigma(1);
    dH_sigma=sigma(2);
    SNR_dH=delta_H/dH_sigma;
    %save results from each fit into a common matrix
    Data_Matrix(file_count,current_frequency)=H_res;
   
    Data_Matrix(file_count,current_frequency+17)=delta_H;
    
    Data_Matrix(file_count,current_frequency+34)=Hrsigma;
    
    Data_Matrix(file_count,current_frequency+51)=dH_sigma;
     
    Data_Matrix(file_count,current_frequency+68)=k1_out; 
    
    Data_Matrix(file_count,current_frequency+85)=k2_out;
    
    %Plot Data
    fighand=figure('Name',['sample=' num2str(filenumber(file_count)) '; f=' num2str(current_frequency) '_GHZ' ]) ;
hold on    
     plot(current_H_data,current_R_data,'.')
     plot(fitresult,H_chopped,R_chopped);
    grid on; legend off; xlabel( 'H' ); ylabel( 'R' );
hold off

cd('output')
saveas(fighand,fighand.Name)
cd('../');
close all
end
end
end

Data_Matrix(:,1)=filenumber;
cd('output')
save('Results.mat',Data_Matrix,filenumber,names_list);

end
 

%%FIT FUNCTION--------------------------
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