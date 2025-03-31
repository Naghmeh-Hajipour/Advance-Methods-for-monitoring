% Variable Declaration - Means and Uncertainties

rng('shuffle');

% Number of sampling
f=10^5;

% Number of Monitored Outputs
nmo=1;

%Number of Inputs
ni=3;

% Distributions of inputs
% m hot
mean1=580;
sigma1=5;

% m cold
mean2=130;
sigma2=2;

% T hot
mean3=523;
sigma3=2;

OUT_M=zeros(f,nmo);
In_M=zeros(f,ni);

tic

for h=1:f
    m_dot_hot = mean1+sigma1*randn; % kg/s
    m_dot_cold = mean2+sigma2*randn; % kg/s
    T_hot_in = mean3+sigma3*randn; % K
    c_p_hot = 1.005; % KJ/kg-K
    c_p_cold = 4.18; % KJ/kg-K
    T_cold_in = 300;
    U=20; %kW/m2.K
    A=10;  %m2
    type ="Parallel Flow";
    [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
    Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
    Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
    LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    Q_exchange = U*A*LMTD;

    OUT_M(h,1) = Q_exchange;
    In_M(h,1) = m_dot_cold;
    In_M(h,2) = T_hot_in;
    In_M(h,3) = T_cold_in;

end

m=mean(OUT_M(:,1));
sigma=std(OUT_M(:,1));

figure();
hold on;
histogram(OUT_M(:,1), 'normalization', 'pdf', 'DisplayName', 'Q exchange distribution')
Normal = fitdist(OUT_M(:,1),'Normal');
Normal_pdf = pdf(Normal,32000:5:34500);
plot (32000:5:34500, Normal_pdf,'LineWidth', 2, 'DisplayName', 'Normal pdf')
xlabel('Q exchange')
ylabel('probability')
title(['distribution and pdf'])
legend("show")
hold off;
