close all; clear all;clc;
%% desired formation + speed given
q1d = [2*sqrt(3);6];
q2d = [0;0];
q3d = [4*sqrt(3);0];

speed_leader_1 = 20;
speed_leader_1_alpha = 20 + 6*rand(); % deg
omega_leader_1 = 1;
speed_leader_1_alpha = speed_leader_1_alpha/180*pi;
%% calculate v2 v3 omega
fun_v2 = @(v2_y) [norm([-2*omega_leader_1;v2_y] - omega_leader_1*[6;-2*sqrt(3)])-speed_leader_1;
                  norm([-2*omega_leader_1;v2_y + 4*sqrt(3)*omega_leader_1] - omega_leader_1*[6;2*sqrt(3)])-speed_leader_1];
v2_y0 = speed_leader_1;
options = optimoptions('fsolve','Display','off');
v2_y = fsolve(fun_v2,v2_y0,options);

v2 = [-2*omega_leader_1;v2_y]; v3 = [-2*omega_leader_1;v2_y + 4*sqrt(3)*omega_leader_1];px_ = [2*sqrt(3);0;4*sqrt(3)]; py_ = [6;0;0];
fun = @(x_) [(x_+x_')*[px_,py_]+[-(v2+v3)';v2';v3'],(x_+x_')*[1;1;1]]; 
x0 = ones(3,3);
x_ = fsolve(fun,x0,options);
omega_rotate = x_ + x_';

%%
omega_leader_1 = 0;
fun_v2 = @(v2_y) [norm([-2*omega_leader_1;v2_y] - omega_leader_1*[6;-2*sqrt(3)])-speed_leader_1;
                  norm([-2*omega_leader_1;v2_y + 4*sqrt(3)*omega_leader_1] - omega_leader_1*[6;2*sqrt(3)])-speed_leader_1];
              v2_y0 = speed_leader_1;
options = optimoptions('fsolve','Display','off');
v2_y = fsolve(fun_v2,v2_y0,options);

v2 = [-2*omega_leader_1;v2_y]; v3 = [-2*omega_leader_1;v2_y + 4*sqrt(3)*omega_leader_1];px_ = [2*sqrt(3);0;4*sqrt(3)]; py_ = [6;0;0];
fun = @(x_) [(x_+x_')*[px_,py_]+[-(v2+v3)';v2';v3'],(x_+x_')*[1;1;1]]; 
x0 = ones(3,3);
x_ = fsolve(fun,x0,options);
omega_straight = x_ + x_';
%%
p_exp = [q1d,q2d,q3d]';
delta_omega = omega_rotate - omega_straight;
v_straight = -omega_straight*p_exp
v_rotate = -delta_omega*p_exp
v_final = -omega_rotate*p_exp
