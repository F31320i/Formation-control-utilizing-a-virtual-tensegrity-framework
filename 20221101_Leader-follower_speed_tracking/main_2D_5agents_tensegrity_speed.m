%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic: Tensegrity-based leader speed follow
%  
%  Time: 2022.8.22
%  
%  Random desired shape
%  Given leader speed (scalar)
% 5 agents
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;
%%

q1d = [0;0];
q2d = [2;0];
q3d = [1;1];
q4d = [0.8;1];
q5d = [1.3;1.5];

speed_leader_3 = 1;
speed_leader_3_alpha = 20 + 6*rand(); % deg
speed_leader_3_alpha = speed_leader_3_alpha/180*pi;

%% omega calculation
qxd = [q1d(1);q2d(1);q3d(1);q4d(1); q5d(1)];
qgd =  (q1d+q2d+q3d+q4d+q5d)/5;% gravity center
speed_leader_cos = (q3d(1) - qgd(1))/norm(q3d - qgd);

k = sqrt(-speed_leader_cos*speed_leader_3/([1 1 -4 1 1]*qxd));
D = k*[1;1;-4;1;1];
omega = D*D';

dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d) norm(q1d-q5d);
        0               0     norm(q2d - q3d) norm(q2d - q4d) norm(q2d-q5d);
        0               0             0       norm(q3d - q4d) norm(q3d-q5d);
        0               0             0                 0     norm(q4d-q5d);
        0               0             0                 0                0];
dijd = dijd+dijd';
%%

x1_ini = q1d+ 1*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 1*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 1*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 1*[rand()-0.5;rand()-0.5];
x5_ini = q5d+ 1*[rand()-0.5;rand()-0.5];


xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini, x5_ini];

%%
dt = 0.001;
steps = 10000;

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4); x5 = xs_ini(:,5);

global expd_c expd_s
expd_c = 1;
expd_s = -1;

x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];
x13s = [0;norm(x1-x3)-dijd(1,3)];x24s = [0;norm(x2-x4)-dijd(2,4)];

x51s = [0;norm(x5-x1)-dijd(5,1)]; x52s = [0;norm(x5-x2)-dijd(5,2)];
x53s = [0;norm(x5-x3)-dijd(5,3)]; x54s = [0;norm(x5-x4)-dijd(5,4)];


x1s = x1;x2s = x2;x3s = x3;x4s = x4; x5s = x5;

% main loop
for t=1:steps
    x = [x1,x2,x3,x4,x5]; 
    v = zeros(2,5);
    for i=1:5
        for j=1:5
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
%     v(:,3) = speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)];
    v = v*5;

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt; x5 = x5+v(:,5)*dt;

    
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4]; x5s = [x5s,x5];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x24s = [x24s, [t*dt;norm(x2-x4)-dijd(2,4)]]; 
    
    x51s = [x51s, [t*dt;norm(x5-x1)-dijd(5,1)]]; x52s = [x52s, [t*dt;norm(x5-x2)-dijd(5,2)]];
    x53s = [x53s, [t*dt;norm(x5-x3)-dijd(5,3)]]; x54s = [x54s, [t*dt;norm(x5-x4)-dijd(5,4)]];
end



figure(1);
set(figure(1),'Position',[200,100,1250,300]);
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
plot(x53s(1,:),x53s(2,:),'LineWidth',1.5);
plot(x54s(1,:),x54s(2,:),'LineWidth',1.5);
legend('edge (1,2)','edge (2,3)','edge (3,4)','edge (4,1)','edge (1,3)','edge (2,4)','edge (5,1)','edge (5,2)','edge (5,3)','edge (5,4)');
xlabel('t/s'); ylabel('║rij║-║rij*║');


subplot(1,2,2);
xs=zeros(4,2,steps+1);
xs(1,:,:) = x1s; xs(2,:,:) = x2s; xs(3,:,:) = x3s; xs(4,:,:) = x4s; xs(5,:,:) = x5s;
for i=1:5
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
    for i=1:5
        for j=1:5
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
    legend([H([1 2 3 4 5]), G, K] ,'agent 1','agent 2','agent 3','agent 4','agent 5','initial configuration','final configuration');
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(dij^(-2*expd_c))*(norm(x1-x2)^(2*expd_c));
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(dij^(-2*expd_s))*(norm(x1-x2)^(2*expd_s));
end
