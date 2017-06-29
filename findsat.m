
function warning_matrix=findsat(inputfolder)

clc;
close all;
fclose all;
listing=dir(pwd)
listing=listing(3:end);

names_list=cell(length(listing),1);
filenumber=zeros(length(listing),1);

warning_matrix=zeros(length(listing)+1,(17+1));
warning_matrix(1,2:end)=2:18;

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
f=f(abs(f-round(f))<.2);
f=round(f);
f=f(f<18.1);
freq_array=unique(f);  
for ii=1:length(freq_array);  
    freq_index=find(freq_array(ii)==f);                             % creates an array of the indecies where a frequency begins and starts- so this code will not work if frequency data is not grouped I am proud of this one
    startindex=freq_index(1);
    lastindex=freq_index(end);
    current_frequency=freq_array(ii);
    current_H_data=H(startindex:lastindex);
    current_R_data=R(startindex:lastindex);
R1=current_R_data;
    Rmax=max(R1);
    index_max=find(R1==Rmax);
    Rmin=min(R1);
    index_min=find(R1==Rmin);
    if index_min>1 & index_max>1 & index_min<length(R1) &index_max>length(R1)
    if ~(abs(R1(index_min)-R1(index_min-1))<.1|abs(R1(index_min)-R1(index_min+1))<0.1| abs(R1(index_max)-R1(index_max-1))<.1 | abs(R1(index_max)-R1(index_max+1))<0.1)
   warning_matrix(filenumber+1,current_frequency);
    end
    end
end
end
end
