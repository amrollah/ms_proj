function [x_new, y, K] = one_KF_step(z_data)
%% function: one_KF_step.m
%
% A. Cortinovis, andrea.cortinovis@ch.abb.com
% ABB CHCRC, 12.08.2013, Project: VISMO
%
% Summary: EKF implementation 
%
% Inputs: z_data = [GHI Temp PV_output] e R^{3x1}
%
% Outputs: x_new e R^{2x1}: KF states
%          y e R: estimated PV output
%          K e R^{1x2}: Kalman Filter matrix 
%
% Dependencies: none
%
% KF Details:
% static parameter evolution using nonlinear output equation and 
% linearization for C around past input. A, R, Q are constant.
%
% x_k+1 = [1 0; 0 1] * x_k with x_k = [alpha1; alpha2]
% y = alpha1 * u + alpha2 * u^2 
% u = GHI
%
% -------------------------------------------------------------------------

persistent x_post P_post R Q A u_old

% initialize function at first function call
if isempty(x_post)
    % Define matrices
    A = eye(2);
    R = 10; % 10;
    Q = eye(2)*1e-10; %eye(2)*1e-10;
    P_post = eye(2)*1e-5; % eye(2)*1e-5;
    P_post(2,2) = 1e-9; %1*1e-9;
    x_post = [3; -0.00053];
    u_old = z_data(1);
end


if ~isnan(z_data(3)) 
    % measured history
    x_prior = A*x_post;
    P_prior = A*P_post*A' + Q;
    C = [u_old u_old^2];
    K = (P_prior * C')/(C*P_prior*C' + R);
    z = z_data(3);
    u_new = z_data(1);
    x_post = x_prior + K * (z - (x_prior(1)*u_new + x_prior(2)*u_new^2));
    P_post = (eye(2) - K*C)*P_prior;
    y = x_post(1)*u_new + x_post(2)*u_new^2;
    x_new = x_post';
    u_old = u_new;
else
    % Predicted horizon
    u_new = z_data(1);
    x_new = A*x_post;
    y = x_new(1)*u_new + x_new(2)*u_new^2;
    K = [NaN NaN];
end