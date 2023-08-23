%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic: Tensegrity-based leader speed follow
%  
%  Time: 2022.9.11
%  
%  Random desired shape
%  Given leader speed (scalar)
%
% 6 agents
% 4 agents undirected graph + 2 follower directed graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;

draw = 0;
%% desired formation + speed given
q1d = [0;0];
q2d = [2;0];
q3d = [2;2];
q4d = [0;2];

q5d = [-4;-2];
q6d = [1;-4];

speed_leader_3 = 1;
speed_leader_3_alpha = 30 + 6*rand(); % deg
speed_leader_3_alpha = speed_leader_3_alpha/180*pi;
%% omega calculation
qxd = [q1d(1);q2d(1);q3d(1);q4d(1)];
qgd =  (q1d+q2d+q3d+q4d)/4;% gravity center
speed_leader_cos = (q3d(1) - qgd(1))/norm(q3d - qgd);
speed_leader_sin = (q3d(2) - qgd(2))/norm(q3d - qgd);

k = sqrt(-speed_leader_cos*speed_leader_3/([1 1 -3 1]*qxd));
D = k*[1;1;-3;1];
omega_4x4 = D*D';
omega = zeros(6,6);
omega(1:4,1:4) = omega_4x4;

dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d) norm(q1d-q5d) norm(q1d-q6d);
        0               0     norm(q2d - q3d) norm(q2d - q4d) norm(q2d-q5d) norm(q2d-q6d);
        0               0             0       norm(q3d - q4d) norm(q3d-q5d) norm(q3d-q6d);
        0               0             0                 0     norm(q4d-q5d) norm(q4d-q6d);
        0               0             0                 0                0  norm(q5d-q6d);
        0               0             0                 0                0              0];
dijd = dijd+dijd';

syms z1 z2 z3; % agent 5 follows agent 1 2 4
wijs_5 =solve(z2 == 0.5, [z1 z2 z3]*[q5d(1) - q1d(1); q5d(1) - q2d(1); q5d(1) - q4d(1)] == speed_leader_3*speed_leader_cos, [z1 z2 z3]*[q5d(2) - q1d(2); q5d(2) - q2d(2); q5d(2) - q4d(2)] == speed_leader_3*speed_leader_sin,z1*z2*z3<0,z1,z2,z3);
omega(5,1) = subs(wijs_5.z1);omega(5,2) = subs(wijs_5.z2);omega(5,4) = subs(wijs_5.z3);

syms zz1 zz2 zz3; % agent 6 followes agent 1 2 5
wijs_6 =solve(zz2 == 0.5, [zz1 zz2 zz3]*[q6d(1) - q1d(1); q6d(1) - q2d(1); q6d(1) - q5d(1)] == speed_leader_3*speed_leader_cos, [zz1 zz2 zz3]*[q6d(2) - q1d(2); q6d(2) - q2d(2); q6d(2) - q5d(2)] == speed_leader_3*speed_leader_sin,zz1*zz2*zz3<0,zz1,zz2,zz3);
omega(6,1) = subs(wijs_6.zz1);omega(6,2) = subs(wijs_6.zz2);omega(6,5) = subs(wijs_6.zz3);


%%
x1_ini = q1d+ 2*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 2*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 2*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 2*[rand()-0.5;rand()-0.5];

x5_ini = q5d+ 2*[rand()-0.5;rand()-0.5];
x6_ini = q6d+ 2*[rand()-0.5;rand()-0.5];

xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini, x5_ini, x6_ini];

%%
dt = 0.001;
steps = 1000;


% draw
if draw ==1
    x_plot = 0;y_plot = 0;figure(1);set(figure(1),'Position',[500,500,550,500]);p = plot(x_plot, y_plot,'o');
end
% for i=1:100
%     pause(0.05);
% end

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4); x5 = xs_ini(:,5); x6 = xs_ini(:,6);

global expd_c expd_s
expd_c = 1;
expd_s = -0.5;

x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];
x13s = [0;norm(x1-x3)-dijd(1,3)];x24s = [0;norm(x2-x4)-dijd(2,4)];

x51s = [0;norm(x1-x5)-dijd(1,5)];x52s = [0;norm(x2-x5)-dijd(2,5)];x54s = [0;norm(x4-x5)-dijd(4,5)];

x61s = [0;norm(x1-x6)-dijd(1,6)];x62s = [0;norm(x2-x6)-dijd(2,6)];x65s = [0;norm(x5-x6)-dijd(5,6)];


x1s = x1;x2s = x2;x3s = x3;x4s = x4;x5s = x5;x6s = x6;


% main loop
for t=1:steps
    x = [x1,x2,x3,x4,x5,x6]; 
    v = zeros(2,6);
    for i=1:6
        for j=1:6
            if i==j
                continue
            end
            if omega(i,j)>0 % strut
                v(:,i) = v(:,i)+force_strut(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            elseif omega(i,j)<0 % cable
                v(:,i) = v(:,i)+force_cable(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            end
        end
    end 
    
%     % speed direction change
%     if mod(t,30)==0
%         speed_leader_3_alpha = speed_leader_3_alpha + 10*(rand()-0.5)/180*pi;
%     end
    
    v(:,3) = speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)];
    v = v*40;
    
%     d_theta = 0.08;
%         speed_leader_3_alpha = speed_leader_3_alpha+(d_theta)/180*pi;
    
    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt; x5 = x5+v(:,5)*dt; x6 = x6+v(:,6)*dt;
    
    % draw
    if draw ==1
        if mod(t,100)==0
        axis equal;
        axis([-2 15 -2 15]);
        x_plot =  [x1(1), x2(1), x3(1), x4(1), x5(1), x6(1)];
        y_plot =  [x1(2), x2(2), x3(2), x4(2), x5(2), x6(2)];
        set (p, 'XData', x_plot, 'YData', y_plot);    
    
        lines = [];
    
        for i=1:6
            for j=1:6
                if omega(i,j)>0 % strut
                    lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)], 'linestyle','-','color','r','LineWidth',1.5)];
                elseif omega(i,j)<0 % cable
                    lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)])];
                end
            end
        end 
        drawnow;
        delete(lines)
        hold off;
        end
    end
   
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4]; x5s = [x5s,x5]; x6s = [x6s,x6];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x24s = [x24s, [t*dt;norm(x2-x4)-dijd(2,4)]]; 
    
    x51s = [x51s, [t*dt;norm(x1-x5)-dijd(1,5)]];x52s = [x52s, [t*dt;norm(x2-x5)-dijd(2,5)]];x54s = [x54s, [t*dt;norm(x4-x5)-dijd(4,5)]];
    x61s = [x61s, [t*dt;norm(x1-x6)-dijd(1,6)]];x62s = [x62s, [t*dt;norm(x2-x6)-dijd(2,6)]];x65s = [x65s, [t*dt;norm(x5-x6)-dijd(5,6)]];
end

figure(1);
set(figure(1),'Position',[200,100,1250,300]);
% 
subplot(1,2,1);
plot(x12s(1,:),x12s(2,:),'LineWidth',1.5);
hold on;
plot(x23s(1,:),x23s(2,:),'LineWidth',1.5);
plot(x34s(1,:),x34s(2,:),'LineWidth',1.5);
plot(x41s(1,:),x41s(2,:),'LineWidth',1.5);
plot(x13s(1,:),x13s(2,:),'LineWidth',1.5);
plot(x24s(1,:),x24s(2,:),'LineWidth',1.5);

plot(x51s(1,:),x51s(2,:),'LineWidth',1.5);
plot(x52s(1,:),x52s(2,:),'LineWidth',1.5);
plot(x54s(1,:),x54s(2,:),'LineWidth',1.5);

plot(x61s(1,:),x61s(2,:),'LineWidth',1.5);
plot(x62s(1,:),x62s(2,:),'LineWidth',1.5);
plot(x65s(1,:),x65s(2,:),'LineWidth',1.5);

legend('edge (1,2)','edge (2,3)','edge (3,4)','edge (4,1)','edge (1,3)','edge (2,4)','d edge (5,1)','d edge (5,2)','d edge (5,4)','d edge (6,1)','d edge (6,2)','d edge (6,5)');
xlabel('t/s'); ylabel('║rij║-║rij*║');


subplot(1,2,2);
xs=zeros(6,2,steps+1);
xs(1,:,:) = x1s; xs(2,:,:) = x2s; xs(3,:,:) = x3s; xs(4,:,:) = x4s; xs(5,:,:) = x5s; xs(6,:,:) = x6s;

for i=1:6
    for j = 1:steps+1
        px(j) = xs(i,1,j);
        py(j) = xs(i,2,j);
        
        if j==1
            G = plot(px(j),py(j),'.','color','g','Markersize',25);
            hold on;
        end
        
        if j==steps+1
            K = plot(px(j),py(j),'.','color','k','Markersize',25);
        end
    end
        H(i) = plot(px,py,'--','LineWidth',2.5);
end
    lines = [];  
    for i=1:6
        for j=1:i
            if omega(i,j)>0 % strut
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                if i == 5 || i ==6 || i==3 || j==3
                    lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], 'linestyle','--','color','r','LineWidth',1.5)];
                    continue
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], 'linestyle','-','color','r','LineWidth',1.5)];
            elseif omega(i,j)<0 % cable
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                if i == 5 || i ==6 || i==3 || j==3
                    lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'linestyle','--','LineWidth',1.5)];
                    continue;
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'LineWidth',1.5)];
            end
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');
    legend([H([1 2 3 4 5 6]), G, K] ,'agent 1','agent 2','agent 3','agent 4','agent5', 'agent6','initial configuration','final configuration');
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(dij^(-2*expd_c))*(norm(x1-x2)^(2*expd_c));
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(dij^(-2*expd_s))*(norm(x1-x2)^(2*expd_s));
end
