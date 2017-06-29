
function [M conf_M]=exp_fit(t,dH_M)
[r,c]=size(dH_M);
M=zeros(3,c);
conf_M=zeros(4,c);

func=fittype('A-B*exp(-t/lambda)','independent','t','dependent','y');

%   func(A,B,lambda,t)
opt=fitoptions(func);
opt.Lower=[0 0 0];
opt.Upper=[300 300 10];
opt.Startpoint=[40 .5 2];
opt.Weights=[1 1 1 1 1 .1 .1 .1 .1 .1 .1];

for i=1:c;
    i
h=dH_M(:,i);

[fitout,gof]=fit(t,h,func,opt);

lambda=fitout.lambda;A=fitout.A; B=fitout.B;
M(1:3,i)=[A; B; lambda];
conf_M(1,i)=gof.rsquare;
conf_temp=confint(fitout,0.682);
conf_M(2:end,i)=M(1:3,i)-conf_temp(1,:)';
fighand=figure('Name', ['f=' num2str(i+1) '10 Weight']);
plot(fitout,t,h);
legend 'off';
end
end