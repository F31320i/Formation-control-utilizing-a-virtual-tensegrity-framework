%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic:  constant repulsion on struts, cables have expl as index
%  On 2D plane four agents, for compariasion with Naomi
%  Time: 2022.3.22
%
%
% Time 2022.8.1
%  Adjusted into distance-based
% find out if distance-based can track leader speed
% 
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;
%%

q1d = [0;0];
q2d = [1;0];
q3d = [1;1];
q4d = [0;1];

dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d);
        0               0     norm(q2d - q3d) norm(q2d - q4d);
        0               0             0       norm(q3d - q4d);
        0               0             0                 0   ];
dijd = dijd+dijd';
%%

x1_ini = q1d+ 0.2*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 0.2*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 0.2*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 0.2*[rand()-0.5;rand()-0.5];

% x1_ini = q1d;
% x2_ini = q2d;
% x3_i3i = q3d;
% x4_ini = q4d;

% rot = [cos(0.1) -sin(0.1);sin(0.1) cos(0.1)];
% x1_ini = rot*q1d+[0.1;0.2];
% x2_ini = rot*q2d+[0.1;0.2];
% x3_ini = rot*q3d+[0.1;0.2];
% x4_ini = rot*q4d+[0.1;0.2];

% x1_ini = [1;0];
% x2_ini = [2;0];
% x3_ini = [3;0];
% x4_ini = [4;0.0];



xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];

% xs_ini = [ 0    1.0000    1.0000   -0.0000;
%             0         0    1.0000    0.9997];
% xs_ini = [3.1472   -2.7301    2.3236   -2.2150;
%           4.0579    4.1338   -3.0246    1.4688];
% 
% xs_ini = [-2.6272    5.6309    1.2114   -0.1110;  % for aij=+-5
%           -0.4115    0.4681   -1.6841    2.2406];

%%
dt = 0.001;
steps = 2000;

% draw
% x_plot = 0;y_plot = 0;z_plot = 0;figure(1);set(figure(1),'Position',[500,500,550,500]);p = plot3(x_plot, y_plot,z_plot,'o');

% for i=1:100
%     pause(0.05);
% end

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4);

x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];
x13s = [0;norm(x1-x3)-dijd(1,3)];x24s = [0;norm(x2-x4)-dijd(2,4)];


x1s = x1;x2s = x2;x3s = x3;x4s = x4;

% main loop
for t=1:steps
    x = [x1,x2,x3,x4]; 
    v = zeros(2,4);
    for i=1:4
        for j=1:4
            if i==j
                continue
            end
            if i==1 && j==3
                continue
            end
            
            if i==3 && j==1
                continue
            end
            
            v(:,i) = v(:,i)+force(x(:,i),x(:,j),dijd(i,j));
        end
    end 
    v(:,3) = [0.1;0];
    v = v*20;

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt;
   
    
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];


    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x24s = [x24s, [t*dt;norm(x2-x4)-dijd(2,4)]]; 
end

figure(1);
set(figure(1),'Position',[200,100,1250,400]);
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
                for k=1:2
                    pi(k) = xs(i,k,steps+1);pj(k) = xs(j,k,steps+1);
                end
                
            if i==1 && j==3
                continue
            end
            
            if i==3 && j==1
                continue
            end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'LineWidth',1.5)];
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');zlabel('z');
    legend([H([1 2 3 4]), G, K] ,'agent 1','agent 2','agent 3','agent 4','initial configuration','final configuration');
    


%% functions
function f = force(x1,x2,dij)
    f = (x2-x1)*(norm(x2-x1)^2 - dij^2);
end
