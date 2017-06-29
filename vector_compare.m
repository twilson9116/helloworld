
[rr cc]=size(Out_new)
adress=zeros(size(Out_new));
for c=1:cc;
for r=1:rr;
 if    any(any(Out_all==Out_new(r,c)));
 else
     adress(r,c)=1;
 end
end
end

test=zeros(164,1);
for i=1:164
   test(i)=any(adress(i,:))
end