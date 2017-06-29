

M_nodH=zeros(8,102);
M_wdH=zeros(8,102);

for i=1:7;
    
    M_nodH=withoutdH(names(i).name);
    M_wdH=withdH(names(i).name);
    
end 