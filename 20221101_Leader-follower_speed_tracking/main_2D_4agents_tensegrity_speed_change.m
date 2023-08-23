%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic: Tensegrity-based leader speed follow
%  
%  Time: 2022.8.19
%  
%  Random desired shape
%  Given leader speed (scalar)
%
%  Speed direction change
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
clc;

%% desired formation + speed given
q1d = [0;0];
q2d = [2;0];
q3d = [2;2];
q4d = [0;2];

speed_leader_3 = 1;
speed_leader_3_alpha = 20 + 6*rand(); % deg
speed_leader_3_alpha = speed_leader_3_alpha/180*pi;
%% omega calculation
qxd = [q1d(1);q2d(1);q3d(1);q4d(1)];
qgd =  (q1d+q2d+q3d+q4d)/4;% gravity center
speed_leader_cos = (q3d(1) - qgd(1))/norm(q3d - qgd);

k = sqrt(-speed_leader_cos*speed_leader_3/([1 1 -3 1]*qxd));
D = k*[1;1;-3;1];
omega = D*D';

dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d);
        0               0     norm(q2d - q3d) norm(q2d - q4d);
        0               0             0       norm(q3d - q4d);
        0               0             0                 0   ];
dijd = dijd+dijd';
%%
x1_ini = q1d+ 5*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 3*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 2*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 3*[rand()-0.5;rand()-0.5];

xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];

%%
dt = 0.001;
steps = 4000;


% draw
x_plot = 0;y_plot = 0;figure(1);set(figure(1),'Position',[500,500,550,500]);p = plot(x_plot, y_plot,'o');

% for i=1:100
%     pause(0.05);
% end

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4);

global expd_c expd_s
expd_c = 1;
expd_s = -0.5;

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
            if omega(i,j)>0 % strut
                v(:,i) = v(:,i)+force_strut(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            elseif omega(i,j)<0 % cable
                v(:,i) = v(:,i)+force_cable(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            end
        end
    end 
    
    if mod(t,300)==0
        speed_leader_3_alpha = speed_leader_3_alpha + 40*(rand()-0.5)/180*pi;
    end
    v(:,3) = speed_leader_3*[cos(speed_leader_3_alpha);sin(speed_leader_3_alpha)];
    v = v*20;

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt;
    
%     % draw
%     if mod(t,100)==0
%     axis equal;
%     axis([-2 35 -2 35]);
%     x_plot =  [x1(1), x2(1), x3(1), x4(1)];
%     y_plot =  [x1(2), x2(2), x3(2), x4(2)];
%     set (p, 'XData', x_plot, 'YData', y_plot);    
%     
%     lines = [];
%     
%     for i=1:4
%         for j=1:4
%             if omega(i,j)>0 % strut
%                 lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)], 'linestyle','-','color','r','LineWidth',1.5)];
%             elseif omega(i,j)<0 % cable
%                 lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)])];
%             end
%         end
%     end
%     drawnow;
%     delete(lines)
%     hold off;
%     end
   
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
    x13s = [x13s, [t*dt;norm(x1-x3)-dijd(1,3)]]; x24s = [x24s, [t*dt;norm(x2-x4)-dijd(2,4)]]; 
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
    legend([H([1 2 3 4]), G, K] ,'agent 1','agent 2','agent 3','agent 4','initial configuration','final configuration');
    
    


%% functions
function pull = force_cable(x1,x2,wij,dij)
    global expd_c;
    pull = (x2-x1)*-wij*(dij^(-2*expd_c))*(norm(x1-x2)^(2*expd_c));
end

function push = force_strut(x1,x2,wij,dij)
    global expd_s;
    push = (x2-x1)*-wij*(dij^(-2*expd_s))*(norm(x1-x2)^(2*expd_s));
end
