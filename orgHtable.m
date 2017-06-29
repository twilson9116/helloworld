function [M S]=organizeHtable(Htable);
T=Htable;
[S,Sii,Si]=unique(T.Sample);
[f,fii,fi]=unique(T.frequency);

M=zeros(length(S),length(f));

for i=1:length(S);
  for   j=1:length(f);
    if isempty(T.H_res(Si==i&T.frequency==f(j)));
     M(i,j)=0;
    else
        M(i,j)= T.H_res(Si==i&T.frequency==f(j));
end
  end
end


M=[f';M]

end 
% 
% 
% [Samples,i1st,iall]=unique(T.Sample);
% [a1,b1,c1]=unique(T.frequency);
% 
% results = cell2mat(accumarray(names, T.value, [], @(x) {x([1 2]).'}));