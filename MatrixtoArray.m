


%r by 17 matrix
BigM
[r c]=size(M);
r=r-1;

L=12;
sets=6;

for jj=1:sets;
    M=BigM(L*(jj-1)+1:jj*L,:)
    M_Hr=M(:,2:18);
    M_dH=M(:,19:35);
    
    for ii=2:18
    Array_Hr(L
    Array_Hr(L*(ii-1),2)=M_Hr(:,ii);
    Array_dH((L*(ii-1),2)=M_dH(:,ii);

    
    
    
    
    Array=zeros(r*18,2);
for i=1:17;  
   % Array((r*(i-1)+1):(r*i),1)=M(2:end,1);
    Array((r*(i-1)+1):(r*i),1)=M(1,i);
    Array((r*(i-1)+1):(r*i),2)=M(2:end,i);
end
