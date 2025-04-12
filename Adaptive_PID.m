clear all;
close all;

% Crane parameters
len1 = 0.5;   
len2 = 0.5;   
mass = 0.5;   
mass1 = 0.25; 
mass2 = 0.5;  
g = 9.81;     

% Initial values
s = 0;      
v = 0;      
th1 = 0;    
th2 = 0;    
w1 = 0;     
w2 = 0;     
f = 0;      
t = 0;     

dt = 0.01;  
xend = 1;   

% PID initial parameters
Kp = 0.5;     
Ki = 0.005;   
Kd = 1;     

% Adaptation rates 
beta_p = 0.0001;
beta_i = 0.00005;
beta_d = 0.0001;

% Storage of states
sdata = s;
vdata = v;
th1data = th1;
th2data = th2;
w1data = w1;
w2data = w2;
fdata = 0;
tdata = t;

% Storage for PID parameters
Kp_data = Kp;
Ki_data = Ki;
Kd_data = Kd;

% Simulation loop
for i = 1:1:2000
    q = [s; th1; th2];
    dq = [v; w1; w2];
    
    M = [mass+mass1+mass2, (mass1+mass2)*len1*cos(th1), mass2*len2*cos(th2);
         (mass1+mass2)*len1*cos(th1), (mass1+mass2)*len1^2, mass2*len1*len2*cos(th1-th2);
         mass2*len2*cos(th2), mass2*len1*len2*cos(th1-th2), mass2*len2^2];
    
    C = [0, -1*(mass1+mass2)*len1*w1*sin(th1), -1*mass2*len2*w2*sin(th2);
         0, 0, mass2*len1*len2*w2*sin(th1-th2);
         0, -1*mass2*len1*len2*w1*sin(th1-th2), 0];
    
    G = [0; (mass1+mass2)*g*len1*sin(th1); mass2*g*len2*sin(th2)];
    
    [f, Kp, Ki, Kd] = adaptivePIDcontroller(Kp, Kd, Ki, xend, s, dt, beta_p, beta_d, beta_i);
    
    T = [f; 0; 0];
    
    ddq = pinv(M) * (T - G - C*dq);
    
    dq_next = dq + dt * ddq;
    q_next = q + dt * dq_next;
    
    s = q_next(1);
    v = dq_next(1);
    th1 = q_next(2);
    w1 = dq_next(2);
    th2 = q_next(3);
    w2 = dq_next(3);
    t = t + dt;
    
    tdata = [tdata, t];
    fdata = [fdata, f];
    vdata = [vdata, v];
    th1data = [th1data, th1];
    th2data = [th2data, th2];
    w1data = [w1data, w1];
    w2data = [w2data, w2];
    sdata = [sdata, s];
    
    Kp_data = [Kp_data, Kp];
    Ki_data = [Ki_data, Ki];
    Kd_data = [Kd_data, Kd];
end

figure;
subplot(3,2,1);
plot(tdata, sdata, 'b');
xlabel('Time (s)');
ylabel('Cart Position (m)');
title('Cart Position s(t)');
grid on;

subplot(3,2,2);
plot(tdata, vdata, 'r');
xlabel('Time (s)');
ylabel('Cart Velocity (m/s)');
title('Cart Velocity v(t)');
grid on;

subplot(3,2,3);
plot(tdata, th1data, 'g'); 
xlabel('Time (s)');
ylabel('\theta_1 (rad)');
title('Pendulum 1 Angle \theta_1(t)');
grid on;

subplot(3,2,4);
plot(tdata, w1data, 'm'); 
xlabel('Time (s)');
ylabel('\omega_1 (rad/s)');
title('Pendulum 1 Angular Velocity \omega_1(t)');
grid on;

subplot(3,2,5);
plot(tdata, th2data, 'c'); 
xlabel('Time (s)');
ylabel('\theta_2 (rad)');
title('Pendulum 2 Angle \theta_2(t)');
grid on;

subplot(3,2,6);
plot(tdata, w2data, 'k'); 
xlabel('Time (s)');
ylabel('\omega_2 (rad/s)');
title('Pendulum 2 Angular Velocity \omega_2(t)');
grid on;

figure;
subplot(2,2,1);
plot(tdata, fdata, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Input f (N)');
title('Control Force over Time');
grid on;

subplot(2,2,2);
plot(tdata, Kp_data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Kp');
title('Proportional Gain Adaptation');
grid on;

subplot(2,2,3);
plot(tdata, Ki_data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Ki');
title('Integral Gain Adaptation');
grid on;

subplot(2,2,4);
plot(tdata, Kd_data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Kd');
title('Derivative Gain Adaptation');
grid on;

function [f, Kp, Ki, Kd] = adaptivePIDcontroller(Kp, Kd, Ki, xend, s, dt, beta_p, beta_d, beta_i)
    persistent e_prev e_int 
    
    if isempty(e_prev)
        e_prev = 0;
        e_int = 0;
    end
    
    e = xend - s;
    e_int = e_int + e * dt;
    e_der = (e-e_prev)/dt;
    
    e_prev = e;
    
    f = Kp * e + Ki * e_int + Kd * e_der;
    
    e0 = 0.01;
    
    e_squared = e^2;
    if (e_squared - e0) == 0
        sgn_term = 0;
    elseif (e_squared - e0) > 0
        sgn_term = 1;
    elseif (e_squared - e0) < 0
        sgn_term = -1;
    end
    
    Kp = Kp + beta_p * (e_squared * sgn_term) * dt;
    Ki = Ki + beta_i * (e^(1/3) * e_int * sgn_term)* dt;
    Kd = Kd + beta_d * (e^(5/3) * e_der * sgn_term)*dt;
    

end