%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic:  constant repulsion on struts, cables have expl as index
%  Time: 2022.4.14
%
%
%
%
%  Conclusion: when expl>0, converge, when expl<0, diverge (proved)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;
%%
q1d = [-sqrt(3)/2;-1/2;0];
q2d = [sqrt(3)/2;-1/2;0];
q3d = [0;1;0];
q4d = [-1/2;-1/2;3];
q5d = [1;0;3];
q6d = [-2/3;8/3-sqrt(3);3*sqrt(3)-2];

scale = 3;

% q1d = [0;0;0]*scale                + scale*[rand()-0.5;rand()-0.5;rand()-0.5];
% q2d = [1;0;0]*scale                + scale*[rand()-0.5;rand()-0.5;rand()-0.5];
% q3d = [1/2;sqrt(3)/2;0]*scale      + scale*[rand()-0.5;rand()-0.5;rand()-0.5];
% q4d = [1/2;-1/(2*sqrt(3));1]*scale + scale*[rand()-0.5;rand()-0.5;rand()-0.5];
% q5d = [1;1/sqrt(3);1]*scale        + scale*[rand()-0.5;rand()-0.5;rand()-0.5];
% q6d = [0;1/sqrt(3);1]*scale        + scale*[rand()-0.5;rand()-0.5;rand()-0.5];

q1d = [0;0;0]*scale;
q2d = [1;0;0]*scale;
q3d = [1/2;sqrt(3)/2;0]*scale;
q4d = [1/2;-1/(2*sqrt(3));1]*scale;
q5d = [1;1/sqrt(3);1]*scale;
q6d = [0;1/sqrt(3);1]*scale;


qd = [q1d,q2d,q3d,q4d,q5d,q6d];
Q = [qd; 1 1 1 1 1 1];
D = null(Q);
% d1 = D(:,1); d2 = D(:,2);
% d1new = d1*d2(1) - d2*d1(1); d1new = d1new/norm(d1new);
% d2new = d1*d2(6) - d2*d1(6); d2new = d2new/norm(d2new);
% D = [d1new,d2new];
omega = D*D';
dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d) norm(q1d - q5d) norm(q1d - q6d);
        0               0     norm(q2d - q3d) norm(q2d - q4d) norm(q2d - q5d) norm(q2d - q6d);
        0               0             0       norm(q3d - q4d) norm(q3d - q5d) norm(q3d - q6d);
        0               0             0                 0     norm(q4d - q5d) norm(q4d - q6d);
        0               0             0                 0              0      norm(q5d - q6d);
        0               0             0                 0              0                   0];
dijd = dijd+dijd';
%%
x1_ini = q1d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];
x2_ini = q2d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];
x3_ini = q3d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];
x4_ini = q4d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];
x5_ini = q5d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];
x6_ini = q6d+ 5*[rand()-0.5;rand()-0.5;rand()-0.5];

% x1_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];
% x2_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];
% x3_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];
% x4_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];
% x5_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];
% x6_ini = [0;0;0]+ 1*[rand()-0.5;rand()-0.5;rand()-0.5];

% x1_ini = q1d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
% x2_ini = q2d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
% x3_ini = q3d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
% x4_ini = q4d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
% x5_ini = q5d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
% x6_ini = q6d+ 10*[rand()-0.5;rand()-0.5;rand()-0.5];
%%
dt = 0.0001;
steps = 10000;

x1 = x1_ini; x2 = x2_ini; x3 = x3_ini; x4 = x4_ini; x5 = x5_ini; x6 = x6_ini;

global expd_c expd_s
expd_c = 1;
expd_s = -0.5;

x12s = [0;norm(x1-x2)-dijd(1,2)];x13s = [0;norm(x1-x3)-dijd(1,3)];x14s = [0;norm(x1-x4)-dijd(1,4)];x23s = [0;norm(x2-x3)-dijd(2,3)];x25s = [0;norm(x2-x5)-dijd(2,5)];
x35s = [0;norm(x3-x5)-dijd(3,5)];x36s = [0;norm(x3-x6)-dijd(3,6)];x45s = [0;norm(x4-x5)-dijd(4,5)];x46s = [0;norm(x4-x6)-dijd(4,6)];x56s = [0;norm(x5-x6)-dijd(5,6)];

x15s = [0;norm(x1-x5)-dijd(1,5)];x26s = [0;norm(x2-x6)-dijd(2,6)];x34s = [0;norm(x3-x4)-dijd(3,4)];

x1s = x1;x2s = x2;x3s = x3;x4s = x4; x5s = x5;x6s = x6;

% main loop
for t=1:steps
    x = [x1,x2,x3,x4,x5,x6]; 
    v = zeros(3,6);
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
    v = v*50;
    
    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt; x5 = x5+v(:,5)*dt; x6 = x6+v(:,6)*dt;
   
    
    
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];x5s = [x5s,x5];x6s = [x6s,x6];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x14s = [x14s, [t*dt;norm(x1-x4)-dijd(1,4)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; x25s = [x25s, [t*dt;norm(x2-x5)-dijd(2,5)]];
    x35s = [x35s, [t*dt;norm(x3-x5)-dijd(3,5)]]; x36s = [x36s, [t*dt;norm(x3-x6)-dijd(3,6)]]; x45s = [x45s, [t*dt;norm(x4-x5)-dijd(4,5)]]; x46s = [x46s, [t*dt;norm(x4-x6)-dijd(4,6)]]; x56s = [x56s, [t*dt;norm(x5-x6)-dijd(5,6)]];

    x15s = [x15s, [t*dt;norm(x1-x5)-dijd(1,5)]]; x26s = [x26s, [t*dt;norm(x2-x6)-dijd(2,6)]]; x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]];
end

figure(1);
set(figure(1),'Position',[200,100,1250,400]);

% 
subplot(1,2,1);
plot(x12s(1,:),x12s(2,:),'LineWidth',1.5);
hold on;
plot(x13s(1,:),x13s(2,:),'LineWidth',1.5);plot(x14s(1,:),x14s(2,:),'LineWidth',1.5);plot(x23s(1,:),x23s(2,:),'LineWidth',1.5);plot(x25s(1,:),x25s(2,:),'LineWidth',1.5);plot(x35s(1,:),x35s(2,:),'LineWidth',1.5);plot(x36s(1,:),x36s(2,:),'LineWidth',1.5);
plot(x45s(1,:),x45s(2,:),'LineWidth',1.5);plot(x46s(1,:),x46s(2,:),'LineWidth',1.5);plot(x56s(1,:),x56s(2,:),'LineWidth',1.5);
legend('cable (1,3)','cable (1,4)','cable (2,3)','cable (2,5)','cable (3,5)','cable (3,6)','cable (4,5)','cable (4,6)','cable (5,6)');
xlabel('t/s'); ylabel('║rij║-║rij*║');

subplot(1,2,2);
xs=zeros(6,3,steps+1);
xs(1,:,:) = x1s; xs(2,:,:) = x2s; xs(3,:,:) = x3s; xs(4,:,:) = x4s; xs(5,:,:) = x5s; xs(6,:,:) = x6s;
for i=1:6
    for j = 1:steps+1
        px(j) = xs(i,1,j);
        py(j) = xs(i,2,j);
        pz(j) = xs(i,3,j);
        
        if j==1
            G = plot3(px(j),py(j),pz(j),'.','color','g','Markersize',25);
            hold on;
        end
        
        if j==steps+1
            K = plot3(px(j),py(j),pz(j),'.','color','k','Markersize',25);
        end
    end
        H(i) = plot3(px,py,pz,'--','LineWidth',2.5);
end
    lines = [];  
    for i=1:6
        for j=1:6
            if omega(i,j)>0 % strut
                for k=1:3
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], [pi(3) pj(3)], 'linestyle','-','color','r','LineWidth',1.5)];
            elseif omega(i,j)<0 % cable
                for k=1:3
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], [pi(3) pj(3)],'LineWidth',1.5)];
            end
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');zlabel('z');
    legend([H([1 2 3 4 5 6]), G, K] ,'agent 1','agent 2','agent 3','agent 4','agent 5', 'agent 6','initial configuration','final configuration');
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(norm(x1-x2)/dij)^(2*expd_c);
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(norm(x1-x2)/dij)^(2*expd_s);
end