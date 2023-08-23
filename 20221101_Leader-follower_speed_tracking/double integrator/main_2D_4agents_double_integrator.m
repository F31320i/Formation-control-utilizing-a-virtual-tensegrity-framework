%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic:  constant repulsion on struts, cables have expl as index
%  On 2D plane four agents, for compariasion with Naomi
%  Time: 2022.10.31 
%  Double integrator model
%
%
%
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;
%%

q1d = [0;0];
q2d = [2;0];
q3d = [2;2];
q4d = [0;2];

qd = [q1d,q2d,q3d,q4d];
Q = [qd; 1 1 1 1];
D = null(Q);
% d1 = D(:,1); d2 = D(:,2);
% d1new = d1*d2(1) - d2*d1(1); d1new = d1new/norm(d1new);
% d2new = d1*d2(6) - d2*d1(6); d2new = d2new/norm(d2new);
% D = [d1new,d2new];
omega = D*D';
dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d);
        0               0     norm(q2d - q3d) norm(q2d - q4d);
        0               0             0       norm(q3d - q4d);
        0               0             0                 0   ];
dijd = dijd+dijd';
%%

x1_ini = q1d+ 5*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 5*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 5*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 5*[rand()-0.5;rand()-0.5];


xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];


%%
dt = 0.001;
steps = 1000;

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4);
v1=[0;0];v2=[0;0];v3=[0;0];v4=[0;0];

global expd_c expd_s
expd_c = 1;
expd_s = -0.5;

% damping
damping = 5;


x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];
x13s = [0;norm(x1-x3)-dijd(1,3)];x24s = [0;norm(x2-x4)-dijd(2,4)];


x1s = x1;x2s = x2;x3s = x3;x4s = x4;


%% main loop

for t=1:steps
    x = [x1,x2,x3,x4];
    v = [v1,v2,v3,v4];
    
    
        a = zeros(2,4);
        for i=1:4
            for j=1:4
                if i==j
                    continue
                end
                if omega(i,j)>0 % strut
                    a(:,i) = a(:,i)+force_strut(x(:,i),x(:,j),omega(i,j),dijd(i,j));
                elseif omega(i,j)<0 % cable
                    a(:,i) = a(:,i)+force_cable(x(:,i),x(:,j),omega(i,j),dijd(i,j));
                end
            end
        end
        a = a*200;
        
        for i=1:4
            a(:,i) = a(:,i) - damping*v(:,i);
        end
        
    %v
    for i=1:4
        v(:,i) = v(:,i) + a(:,i)*dt;
    end 
    v = v*200;
        
    
    % x
    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt;
    
    % l_{ij}
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];
    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x24s = [x24s, [t*dt;norm(x2-x4)-dijd(2,4)]]; 
end





















%% post processing
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
legend('edge (1,2)','edge (2,3)','edge (3,4)','edge (4,1)','edge (1,3)','edge (2,4)');
xlabel('t/s'); ylabel('║rij║-║rij*║');


subplot(1,2,2);
xs=zeros(4,2,steps+1);
xs(1,:,:) = x1s; xs(2,:,:) = x2s; xs(3,:,:) = x3s; xs(4,:,:) = x4s;
for i=1:4
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
    for i=1:4
        for j=1:4
            if omega(i,j)>0.0001 % strut
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], 'linestyle','-','color','r','LineWidth',1.5)];
            elseif omega(i,j)<-0.0001 % cable
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'LineWidth',1.5)];
            end
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');zlabel('z');
    legend([H([1 2 3 4]), G, K] ,'agent 1','agent 2','agent 3','agent 4','initial configuration','final configuration');
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(norm(x1-x2)/dij)^(2*expd_c);
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(norm(x1-x2)/dij)^(2*expd_s);
end
