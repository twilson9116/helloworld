function [Cout,fit_all]=doit(Cin)
Cout=Cin;
for i=1:8
    
fit=fit_single_sweep(Cin{i});
fit_all{i}=fit;
Cout{i}(:,3)=fit(Cin{i}(:,1));

end 

end 