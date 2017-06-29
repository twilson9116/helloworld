%% FMR_spectrum: reads '-.log' files created by NanOsc FMR software. 
%The input can be a single text file or a folder with multiple files. The
%input is file name, folder name, or folder directory as string ('encolse
%with single quotes'). The output is a table that is automaticaly written
%to the current folder.

function H_Table=FMR_spectrum(input,plot_sweeps)
clc; close all; fclose all;
%% SETTINGS:
if nargin<2
    plot_sweeps=false;%'plot_sweeps' is logical that controls if the data and fits are plotted
end
rsquare_min=0.9;
Use_rigorous_fit=false;
rsquare_for_estimate=0.5;
fit_width=1;

%SETTINGS DESCRIPTIONS
%'rsquare_min' controls the minimum allowable quality of the fit before
%proceeding. Smaller values will make program run faster. Larger values
%improve reliability of the fits. Max value=1.
%'Use_rigorous_fit' is a logical if set to true, preliminary fit is done to
%find an approximate linewidth. The final fit is then done on
%data only near the spectrum. The number of data used in final fit is
%controled by 'fit_width'. See below. Use this if there is a large
%amount of data not in the spectrum
%'rsquare_for_estimate' is used for the initial fit that is used to find
%linewidth and initial guesses.
%fit_width controls how much of the data will be fit from the spectrum.
%The number is a multiplier of the linewidth.
%For example if this =1 than and the linewidth is 20 Oe, Only data 20 Oe
%above and below the resonance field will be used for the final fitting.
%If all data must be included, you may use an arbitrarily large number(ex. fit_width=10000)

if plot_sweeps
    mkdir('plots_FMR');
end

%% DISTINGUISH FILE OR FOLDER; IF FOLDER PROCESS ALL FILES
if isdir(input);
    file_structure=dir(input); %Outputs all files and folders in the input directory as a structure. the first two entries are holders for parent directories and use '.' and '..'
    file_list={file_structure(3:end).name}';%this saves just the file names as strings 
    addpath(input); %This places the input directory in Matlabs search path. it cna now find the files
    %H_Matrix=zeros(4000,6); %!!!!!I need to pre-allocate H_table to save
    %space. 
    for j=1:length(file_list);
        display('Processing:')
        display(file_list(j))
        T=FMR_single_sample(file_structure(j+2).name,plot_sweeps, rsquare_min, Use_rigorous_fit, rsquare_for_estimate, fit_width); %Sample,frequency,H_res,deltaH,k_1,k_2,R_square_fit
        if j==1
            [r,~]=size(T);
            H_Table=T;
        else
            H_Table=vertcat(H_Table,T);
        end
    end
else
    H_Table=FMR_single_sample(input,plot_sweeps, rsquare_min, Use_rigorous_fit, rsquare_for_estimate, fit_width);
end
writetable(H_Table,strcat('H_table_',datestr(now,'mmdd_HH_MM')));
end
function SingleSampleTable=FMR_single_sample(input_file,plot_sweeps, rsquare_min, Use_rigorous_fit, rsquare_for_estimate, fit_width);
%% FIND DATA AND SAVE TO MATRIX
%the below for-loop finds the start of the data & bypasses heading
fileID=fopen(input_file); %opens text file written by NanOsc software
startline=0;
stop=0;
while stop~=1
    tline=fgetl(fileID);
    s1=tline;
    s2='[Data]';
    stop=strcmp(s1,s2);
    startline=startline+1;
end
fclose(fileID); %Closes file. 'fopen' is used for 'fgetl' to find start of data. 'dlmread' is used to read data into Matlab.
startline=startline+2;

data_matrix=dlmread(input_file,'\t',startline,0); %reads in data and saves as matrix.
rows=sortrows(data_matrix);% This sorts all rows in ascending order. This does not destinguish seperate sweeps using the same frequency.
f=rows(:,1); %frequency, 1st column
H=rows(:,2); %field data, 2nd column
R=rows(:,5); %Absorption data corrected phase and drift, 5th column
%I=row(:,3); %In-phase data, 3rd column
%Q=row(:,4); %Quadrature data, 4th column
%This code currently uses the phase and drift corrected data. I
% and Q are used in a version that corrects phase shift of quadrature data
% and also for the linear drift that is observed with some data.
%Contact Tom White for the additional code.

unique_f=unique(f);
%Below, all output variables are pre-allocated to matrices or cell;
frequency=unique_f;
H_res=zeros(length(unique_f),1);
H_r_stdv=zeros(length(unique_f),1);
deltaH=zeros(length(unique_f),1);
dH_stdv=zeros(length(unique_f),1);
k_1=zeros(length(unique_f),1);
k_2=zeros(length(unique_f),1);
R_square_fit=zeros(length(unique_f),1);
Sample=cell(length(unique_f),1);
Sample(:)={input_file};
% stop FIND DATA AND SAVE TO MATRIX
%% ISOLATE DATA BY FREQUENCY AND FIT
for i=1:length(unique_f);
    try
        current_rows=find(unique_f(i)==f);%creates an array of the indecies where a frequency begins and starts- so this code will not work if frequency data is not grouped.
        current_frequency=unique_f(i);
        display(strcat(num2str(unique_f(i)),'GHz'))
        current_H_data=H(current_rows(1):current_rows(end));
        current_R_data=R(current_rows(1):current_rows(end));
        
        [~,index_min]=min(current_R_data);%finds the index for the minimum absorption value
        [~,index_max]=max(current_R_data);%find index for the max value
        low_index=min(index_min,index_max);%finds which index is higher value, i.e. higher field
        high_index=max(index_min,index_max);
        
        %Intitial Guesses for fits:
        if i==1
            Hr=((current_H_data(high_index)-current_H_data(low_index))/2)+current_H_data(low_index);
            dH=(sqrt(3)/2).*(current_H_data(high_index)-current_H_data(low_index));
            k1=-1.5.*(Hr./dH).^3;
            k2=(Hr./dH).^2;
        else
            Hr=fit2.Hr;
            dH=fit2.delta_H;
            k1=fit2.k1;
            k=fit2.k2;
        end
        if Use_rigorous_fit;
            %If 'Use_rigorous_fit' is set to 'true', does an initial fit to
            %improve guesses and to focus fit only on data in the spectrum.
            %If set to false, skips to an initieal fit and all data is fit instead of
            %only data near the resonance field and the above initial guesses are used.
            
            [fit1,~]=fit_FMR_spectrum(current_H_data,current_R_data,Hr,dH,k1,k2,rsquare_for_estimate);
            %Initial fit to obtain initial guesses and to narrow the data so that only
            %data near the resoance is fit to the function.
            
            Hr=fit1.Hr; %Overwrites initial guesses for parameters based on initial fit.
            dH=fit1.delta_H;
            k1=fit1.k1;
            k2=fit1.k2;
            
            H_up=Hr+fit_width*dH; % magnitude of upper and lower fields that will be chopped
            H_low=Hr-fit_width*dH;
            Upper_chop=abs(current_H_data-H_up);% subtract scalar field value of chop locations so that it equals minimum absoluite value
            Lower_chop=abs(current_H_data-H_low);
            [~, indexup]=min(Upper_chop);
            [~ ,indexlow]=min(Lower_chop);
            R_to_fit=current_R_data(indexup:indexlow);
            H_to_fit=current_H_data(indexup:indexlow);
        else
            R_to_fit=current_R_data;
            H_to_fit=current_H_data;
        end
        
        [fit2,gof2]=fit_FMR_spectrum(H_to_fit,R_to_fit,Hr,dH,k1,k2,rsquare_min)  ;
       
        fit2_bounds=confint(fit2,0.68);
       
        H_r_stdv(i)=fit2.Hr-fit2_bounds(1,1);    
        H_res(i)=fit2.Hr;
        deltaH(i)=fit2.delta_H;
        dH_stdv(i)=fit2.delta_H-fit2_bounds(1,2);
        k_1(i)=fit2.k1;
        k_2(i)=fit2.k2;
        R_square_fit(i)=gof2;
    
    
   
    %% PLOT DATA
    if plot_sweeps
        figure('Name',[input_file(1:end-22), '_f_' num2str(current_frequency)]);
        %input_file(1:end-22) removes the last 22 charechters so only
        %sample name remains. delete (1:end-22) for entire filename to be
        %in figure title
        hold on;
        plot(current_H_data,current_R_data,'.k');
        plot(H_to_fit,R_to_fit,'ob')
        plot(fit2,'-r');
        ylabel 'Absorption [Arb.]';
        xlabel 'Field [Oe]';
        legend(['All data'],['Data for fit'] , ['Fit Line'])
        grid on;
        hold off ;
    end %ends plot_sweeps if statment
    % stop PLOT DATA
  
    
    catch
        display('^^Fit failed for this frequency. Values stored as 0.')
    end % stop ISOLATE DATA BY FREQUENCY AND FIT
end %Ends loop when all sweeps have been processed


SingleSampleTable=table(Sample,frequency,H_res,H_r_stdv,deltaH,dH_stdv,k_1,k_2,R_square_fit); %the output is a Table.

if plot_sweeps
    figs=get(0,'children');
for i=1:length(figs)
    saveas(figs(i),[pwd '/plots_FMR/' figs(i).Name])
end
close all %If you have a large number of sweeps and MAtlab crashes due to memory issues- save and close each figure at each iteration. Just move this saveas function into the primary for loop
end
end %End of Code
function [best_fit,best_gof]=fit_FMR_spectrum(H_data,Absorption,Hr_guess,dH_guess,k1_guess,k2_guess,gof_setting)
%% FIT FUNCTION AND LOOP
%This function contains the fit equation and iteration loop

[H_prepared,Amplitude_prepared]=prepareCurveData(H_data,Absorption);

fit_equation=fittype('k1*(4*(x-Hr))/(delta_H^2+4*(x-Hr)^2)^2+k2*((delta_H^2-4*(x-Hr)^2)/(delta_H^2+4*(x-Hr)^2)^2)','independent', 'x', 'dependent', 'y' );
%Lorenzian derivative fit equation. If a drift is seen in data add a slope
%and intercept term to the equation.

%Fit options: order of variables is (Hr,deltaH,k1,k2)
opts = fitoptions(fit_equation);
opts.Display = 'Off';
opts.Lower = [0 0 -Inf -Inf];
opts.StartPoint = [Hr_guess dH_guess k1_guess k2_guess];
opts.Upper = [20000 1500 Inf Inf];
opts.MaxIter=8000; %Controls maximum number of internal loops for Matlab's fit function.
[fit_try,gof_current]=fit(H_prepared,Amplitude_prepared,fit_equation,opts);

%check Goodness of Fit
best_gof=gof_current.rsquare;
best_fit=fit_try;
if best_gof<gof_setting;
    k1_array = k1_guess.*[-1 10 -10];
    k2_array = k2_guess.*[1 -1 100 -100];
    
    for ii = 1: length(k1_array)
        
        for j=1:length(k2_array)
            % sprintf('Best GOF= %.3f \nCurrent GOF= %.3f',best_gof,gof_current.rsquare)
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

