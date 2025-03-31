data=xlsread('WindDataInExcel_Open.xls', 'Dati','D3:D62134');
data=data(data ~= 0);
rng('shuffle');

% Number of sampling
ns=100000;

% Number of Monitored Outputs
nmo=1;

%Number of Inputs
ni=4;

% Distributions of inputs
% A
mean1=100;
sigma1=0.2;
% V
Burr=fitdist(data,"Burr");

OUT_M=zeros(ns,nmo);
In_M=zeros(ns,ni);

tic

for h=1:ns

    Cp=0.39;
    ro=1.225;
    A=mean1+sigma1*randn;
    v=random(Burr);
    if v<0
        v=0;
    end
    P=Cp*0.5*ro*A*v^3;

    OUT_M(h,1)=P;
    In_M(h,1)=v;
    In_M(h,2)=Cp;
    In_M(h,3)=A;
    In_M(h,4)=ro;

end %MC
m=mean(OUT_M(:,1));
sigma=std(OUT_M(:,1));
figure(1);
hold on;
histogram(OUT_M(:,1),'normalization','pdf')
xlabel('power')
ylabel('probability')
title(['4th run for ' num2str(ns) ' samples'])
hold off;