%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic: Tensegrity-based leader speed follow
%  
%  Time: 2022.11.16
%  
% 3 agents (1 leader + 2 follower) 
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;

%% desired formation + speed given
q1d = [0;0];
q2d = [6;-2*sqrt(3)];
q3d = [6;2*sqrt(3)];

speed_leader_3 = 20;
speed_leader_3_alpha = 30; % deg
speed_leader_3_alpha = speed_leader_3_alpha/180*pi;
%% omega calculation
qxd = [q1d(1);q2d(1);q3d(1)];
qgd =  (q1d+q2d+q3d)/3;% gravity center
speed_leader_cos = (q3d(1) - qgd(1))/norm(q3d - qgd);

k = sqrt(-speed_leader_cos*speed_leader_3/([1 1 -2]*qxd));
D = k*[1;1;-2];
omega = D*D';

dijd = [0     norm(q1d - q2d) norm(q1d - q3d);
        0               0     norm(q2d - q3d);
        0               0             0       ];
dijd = dijd+dijd';
%%
% x1_ini = q1d+ 10*[rand()-0.5;rand()-0.5];
% x2_ini = q2d+ 10*[rand()-0.5;rand()-0.5];
% x3_ini = q3d+ 10*[rand()-0.5;rand()-0.5];

x1_ini = q1d+ 1*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 1*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 1*[rand()-0.5;rand()-0.5];

xs_ini = [x1_ini, x2_ini, x3_ini];

%%
dt = 0.001;
steps = 4000;

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3);

global expd_c expd_s 
expd_c = 1;
expd_s = -0.5;



x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x13s = [0;norm(x1-x3)-dijd(1,3)];


x1s = x1;x2s = x2;x3s = x3;

dv1s =[];
dv2s =[];

% main loop
for t=1:steps
    x = [x1,x2,x3]; 
    v = zeros(2,3);
    for i=1:3
        for j=1:3
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
%     % velocity direction change
%     if mod(t,steps/4)==0
%         speed_leader_3_alpha = speed_leader_3_alpha + 80*(rand()-0.5)/180*pi;
%     end
    
    v(:,3) = speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)];
    
%     speed_leader_3_alpha = speed_leader_3_alpha + 0.1/180*pi;
    

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt;
  
   
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; 
    
    
    t0=0;
    if t>=t0
        dv1s = [dv1s,[(t-t0)*dt;norm(v(:,1) - speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)])]];
        dv2s = [dv2s,[(t-t0)*dt;norm(v(:,2) - speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)])]];
    end
end

figure(1);
set(figure(1),'Position',[200,100,1550,300]);
% 
subplot(1,3,1);
plot(x12s(1,:),x12s(2,:),'LineWidth',1.5);
hold on;
plot(x23s(1,:),x23s(2,:),'LineWidth',1.5);
plot(x13s(1,:),x13s(2,:),'LineWidth',1.5);

grid on;
legend('edge (1,2)','edge (2,3)','edge (1,3)');
xlabel('t/s'); ylabel('║rij║-║rij*║');


subplot(1,3,2);
plot(dv1s(1,:),dv1s(2,:),'LineWidth',1.5);
hold on;
plot(dv2s(1,:),dv2s(2,:),'LineWidth',1.5);
grid on;
legend('agent 1', 'agent 2');
xlabel('t/s'); ylabel('║vi - v*║');



subplot(1,3,3);
xs=zeros(4,2,steps+1);
xs(1,:,:) = x1s; xs(2,:,:) = x2s; xs(3,:,:) = x3s;
for i=1:3
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
    for i=1:3
        for j=1:3
            if omega(i,j)>0 % strut
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], 'linestyle','-','color','r','LineWidth',1.5)];
            elseif omega(i,j)<0 % cable
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'LineWidth',1.5)];
            end
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');zlabel('z');
    legend([H([1 2 3]), G, K] ,'agent 1','agent 2','agent 3','initial configuration','final configuration');
    grid on;
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(dij^(-2*expd_c))*(norm(x1-x2)^(2*expd_c));
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(dij^(-2*expd_s))*(norm(x1-x2)^(2*expd_s));
end
