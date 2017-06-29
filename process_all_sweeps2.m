function [Out,names]=process_all_sweeps2(Folder)
listing=dir(Folder);%list the contesnts of the folder 
listing=listing(3:end); %this makes the list only the text files. listing(1:2) are "." and ".." 
cd(Folder); %changes the directory into the folder
errorlist=cell(length(listing)+1,1);
completed_list=cell(length(listing)+1,1);
Matrix=zeros(length(listing)+1,68);
Matrix(1,:)=[2:18 2:18 2:18 2:18];
listing
for k=1:length(listing);
disp(k)
try 
Matrix(k+1,:)=process_sweeps_trouble(listing(k).name)
completed_list{k+1}=listing(k).name;
save('all_sweeps_proccessed.mat','Matrix','completed_list','errorlist')
catch 
 errorlist{k+1}=listing(k).name;   
end 
end
