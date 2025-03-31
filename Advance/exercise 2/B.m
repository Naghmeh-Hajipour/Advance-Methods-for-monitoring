
% Stochastic input
nv=3;                    % number of stochastic input variables

% m hot
mean1=580;
sigma1=5;

% m cold
mean2=130;
sigma2=2;

% T hot
mean3=523;
sigma3=2;


% Means vector creator
m_values=zeros(1,nv);
m_values(1,:)=[mean1,mean2,mean3];

% Standard deviation vector creator
sigmavalues=zeros(1,nv);
sigmavalues(1,:)=[sigma1,sigma2,sigma3];

% Variances vector creator
v_values=zeros(1,nv);
v_values(1,:)=sigmavalues.^2;

step=10;
T_cold_in = 300;
type ="Parallel Flow";
c_p_hot = 1.005; % KJ/kg-K
c_p_cold = 4.18; % KJ/kg-K
U=20; %kW/m2.K
A=10;  %m2

delta=zeros(1,nv);
delta(1,:)=[step*mean1,step*mean2,step*mean3]; % Step size of variables

% Calculator of the first order mean for outputs values
% Mean values
m_dot_hot = m_values(1);
m_dot_cold = m_values(2);
T_hot_in = m_values(3);

% running Main
[T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
Q_exchange = U*A*LMTD;

Q = Q_exchange;

CD=zeros(1,nv);
CD2=zeros(1,nv);

%Matrix initialization
Q_p=zeros(1,nv);
Q_pp=zeros(1,nv);
Q_m=zeros(1,nv);
Q_mm=zeros(1,nv);

for t=1:nv

    disp(t)

    %Output plus
    m_values(1,:)=[mean1,mean2,mean3];
    m_values(t)=m_values(t)+delta(t);
    m_dot_hot = m_values(1);
    m_dot_cold = m_values(2);
    T_hot_in = m_values(3);

    % running Main
    [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
    Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
    Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
    LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    Q_exchange = U*A*LMTD;

    Q_p(t) = Q_exchange;

    %Output minus
    m_values(1,:)=[mean1,mean2,mean3];
    m_values(t)=m_values(t)-delta(t);
    m_dot_hot = m_values(1);
    m_dot_cold = m_values(2);
    T_hot_in = m_values(3);

    % running Main
    [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
    Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
    Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
    LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    Q_exchange = U*A*LMTD;

    Q_m(t)= Q_exchange;

    %Output plus_plus
    m_values(1,:)=[mean1,mean2,mean3];
    m_values(t)=m_values(t)+2*delta(t);
    m_dot_hot = m_values(1);
    m_dot_cold = m_values(2);
    T_hot_in = m_values(3);

    % running Main
    [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
    Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
    Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
    LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    Q_exchange = U*A*LMTD;

    Q_pp(t) = Q_exchange;

    %Output minus_minus
    m_values(1,:)=[mean1,mean2,mean3];
    m_values(t)=m_values(t)-2*delta(t);
    m_dot_hot = m_values(1);
    m_dot_cold = m_values(2);
    T_hot_in = m_values(3);

    % running Main
    [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,type);
    Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out);
    Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out);
    LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    Q_exchange = U*A*LMTD;

    Q_mm(t) = Q_exchange;

    % Central-Difference Formulas for first order derivative
    %     Dg_kp_min=(kp_min_plus(t)-kp_min_minus(t))/(2*delta(t)); %O(h2)
    Dg_kp_min=(Q_mm(t)-8*Q_m(t)+8*Q_p(t)-Q_pp(t))/(12*delta(t)); %O(h4)

    % Central-Difference Formulas for second order derivative
    %     Dg2_kp_min=(kp_min_plus(t)-2*kp_min+kp_min_minus(t))/(delta(t)^2); %O(h2)
    Dg2_kp_min(1,1)=((16*Q_m(t))-(Q_mm(t))-(30*Q)+(16*Q_p(t))-(Q_pp(t)))/(12*delta(t)*delta(t)); %O(h4)

    CD(1,t)=Dg_kp_min;
    CD2(1,t)=Dg2_kp_min;

    if t==nv
        save(['workspace_total_s', num2str(step),'.mat'])
    end

end

v_Q = 0;

for r=1:nv
    v_values(1,:)=sigmavalues.^2;
    v_Q=v_Q+(((CD(1,r))^2)*v_values(1,r));
end

sigmaQ=sqrt(v_Q);

% Second Order mean
s_Q=0;

for r=1:nv
    v_values(1,:)=sigmavalues.^2;
    s_Q=s_Q+((CD2(1,r))*v_values(1,r));
end

m_Q=Q+(0.5*s_Q);     %second order
% m_kp_min=kp_min;                     %first order

save(['Prob_values_total_s', num2str(step),'.mat'],'m_Q','sigmaQ')
% pdf plotter of kp
X1=((m_Q-5*sigmaQ):(m_Q+5*sigmaQ));
figure(1)
plot(X1,normpdf(X1,m_Q,sigmaQ),'b');
title(['normal pdf of Q, step = ',num2str(step)])
h1=xlabel('Q exchange');
h2=ylabel('probability');


% sensitivity analysis
nominal_in=[mean1,mean2,mean3];
nominal_out=Q;

sensitivity_Q=CD(1,:).*(nominal_in/nominal_out);

figure(2)
bar(sensitivity_Q);
title(['Q Sensitivity, step = ',num2str(step)])
set(gca,'XTickLabel',{'m hot','m cold','T hot'})