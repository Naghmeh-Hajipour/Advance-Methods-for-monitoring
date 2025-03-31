data=xlsread('WindDataInExcel_Open.xls', 'Dati','E3:E62134');
data=data(data ~= 0);
m=mean(data);
sigma=std(data);
figure(1);
hold on;
histogram(data,'normalization','pdf','DisplayName','real')
Normal=fitdist(data,"Normal");
Normal_pdf=pdf(Normal,0:250);
plot(0:250,Normal_pdf,'DisplayName',"Normal")
Poisson=fitdist(data,"Poisson");
Poisson_pdf=pdf(Poisson,0:250);
plot(0:250,Poisson_pdf,'DisplayName',"Poisson")
Lognormal=fitdist(data,"Lognormal");
Lognormal_pdf=pdf(Lognormal,0:250);
plot(0:250,Lognormal_pdf,'DisplayName',"Lognormal")
Rician=fitdist(data,"Rician");
Rician_pdf=pdf(Rician,0:250);
plot(0:250,Rician_pdf,'DisplayName',"Rician")
Burr=fitdist(data,"Burr");
Burr_pdf=pdf(Burr,0:250);
plot(0:250,Burr_pdf,'DisplayName',"Burr")
Logistic=fitdist(data,"Logistic");
Logistic_pdf=pdf(Logistic,0:250);
plot(0:250,Logistic_pdf,'DisplayName',"Logistic")
Lognormal=fitdist(data,"Lognormal");
Lognormal_pdf=pdf(Lognormal,0:250);
plot(0:250,Lognormal_pdf,'DisplayName',"Lognormal")
Weibull=fitdist(data,"Weibull");
Weibull_pdf=pdf(Weibull,0:250);
plot(0:250,Weibull_pdf,'DisplayName',"Weibull")
xlabel('velocity');
ylabel('probability');
title('velocity pdf');
legend("show")
hold off