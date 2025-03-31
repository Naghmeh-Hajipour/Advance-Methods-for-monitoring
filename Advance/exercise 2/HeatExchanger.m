function [T_hot_out,T_cold_out,epsilon]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,HE_Type);
% [T_hot_out,T_cold_out]=HeatExchanger(c_p_hot,m_dot_hot,T_hot_in,c_p_cold,m_dot_cold,T_cold_in,U,A,HE_Type);
% This function calculates the outlet temperatures of a heat exchanger
% using Epsilon-NTU method. This function uses effectiveness.m as a
% function and should have access to that function.
% 
% The inputs are as follows:
% Hot Flow: c_p_hot, m_dot_hot, T_hot_in.
% Cold Flow: c_p_cold, m_dot_cold, T_cold_in.
% Heat exchanger design parameters: U,A, HE_Type.
% 
% HE_Type defines the type of heat exchanger: (see reference)
%   'Parallel Flow'
%   'Counter Flow'
%   'One Shell Pass'
%   'N Shell Pass'
%   'Cross Both Unmixed'
%   'Cross Cmax Mixed'
%   'Cross Cmin Mixed'
%
% Reference:
% Frank P. Incropera, Introduction to heat transfer. New York:Wiley, 1985, Section 11.4.
% Programmer: Seyyed Ali Hedayat Mofidi (seyyed4li@yahoo.com)
C_hot = m_dot_hot*c_p_hot;
C_cold = m_dot_cold*c_p_cold;
C_min = min(C_hot,C_cold); % finds the flow with lower heat capacity and higher temperature change.
C_max = max(C_hot,C_cold); % finds the flow with higher heat capacity and lower temperature change. 
C_r=C_min/C_max;
NTU = U*A/C_min;
epsilon = effectiveness (NTU,C_r,HE_Type);
Q_max = C_min*(T_hot_in-T_cold_in);
Q = epsilon * Q_max ;
T_hot_out = T_hot_in - Q/C_hot ;
T_cold_out = T_cold_in + Q/C_cold ;
end
