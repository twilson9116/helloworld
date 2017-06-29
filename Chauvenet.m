%%uses Chauvenet's criterion to remove outliers
function new_v=Chauvenet(v)
new_v=v;
new_v(new_v==0)=NaN;


stddv=std(new_v,'omitnan');
avg=mean(new_v,'omitnan');
n=nnz(~isnan(new_v));
Dmax=norminv(1-1/4/n);
v_test=(abs(new_v-avg)/stddv);%ratio of deviation to std. dev. 

iterations=0;
while any(v_test>Dmax)
  
stddv=std(new_v,'omitnan');
avg=mean(new_v,'omitnan');
n=nnz(~isnan(new_v));
Dmax=norminv(1-1/4/n);
v_test=(abs(new_v-avg)/stddv);
    
   
new_v(v_test>Dmax)=NaN;


iterations=1+iterations
end

new_v

end