%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Topic:  constant repulsion on struts, cables have expl as index
%  On 2D plane four agents, for compariasion with Naomi
%  Time: 2022.3.22
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
q2d = [1;0];
q3d = [1;1];
q4d = [0;1];

% q1d = [0;0] + 1*[rand()-0.5;rand()-0.5];
% q2d = [1;0] + 1*[rand()-0.5;rand()-0.5];
% q3d = [1;1] + 1*[rand()-0.5;rand()-0.5];
% q4d = [0;1] + 1*[rand()-0.5;rand()-0.5];

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

x1_ini = q1d+ 10*[rand()-0.5;rand()-0.5];
x2_ini = q2d+ 10*[rand()-0.5;rand()-0.5];
x3_ini = q3d+ 10*[rand()-0.5;rand()-0.5];
x4_ini = q4d+ 10*[rand()-0.5;rand()-0.5];

% x1_ini = q1d;
% x2_ini = q2d;
% x3_ini = q3d;
% x4_ini = q4d;

% rot = [cos(0.1) -sin(0.1);sin(0.1) cos(0.1)];
% x1_ini = rot*q1d+[0.1;0.2];
% x2_ini = rot*q2d+[0.1;0.2];
% x3_ini = rot*q3d+[0.1;0.2];
% x4_ini = rot*q4d+[0.1;0.2];

xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];

% xs_ini = [3.1472   -2.7301    2.3236   -2.2150;
%           4.0579    4.1338   -3.0246    1.4688];
% 
% xs_ini = [-2.6272    5.6309    1.2114   -0.1110;  % for aij=+-5
%           -0.4115    0.4681   -1.6841    2.2406];
%%
dt = 0.0005;
steps = 2000;

% draw
% x_plot = 0;y_plot = 0;z_plot = 0;figure(1);set(figure(1),'Position',[500,500,550,500]);p = plot3(x_plot, y_plot,z_plot,'o');

% for i=1:100
%     pause(0.05);
% end

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4);

global expd
expd = 1;

x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];
x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];


x1s = x1;x2s = x2;x3s = x3;x4s = x4;

x_now = [x1,x2,x3,x4];
v_ini = 0;
    for i=1:4
        for j=1:4
            if omega(i,j)>0 % strut
                v_ini = v_ini+ (-omega(i,j))*dijd(i,j)*norm(x_now(:,i) - x_now(:,j));
            elseif omega(i,j)<0 % cable
                v_ini = v_ini+ 0.5*1/(expd+1)*(-omega(i,j)/(dijd(i,j)^(2*expd)))*(norm(x_now(:,i) - x_now(:,j))^2)^(expd+1);
            end
        end
    end 
vs = [0;v_ini];

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
    v = v*20;

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt;
    x_now = [x1,x2,x3,x4];

    xs = [x1, x2, x3, x4];
    Hessian = cal_hessian(xs);
    
    v_ini = 0;
    for i=1:4
        for j=1:4
            if omega(i,j)>0 % strut
                v_ini = v_ini+ (-omega(i,j))*dijd(i,j)*norm(x_now(:,i) - x_now(:,j));
            elseif omega(i,j)<0 % cable
                v_ini = v_ini+ 0.5*1/(expd+1)*(-omega(i,j)/(dijd(i,j)^(2*expd)))*(norm(x_now(:,i) - x_now(:,j))^2)^(expd+1);
            end
        end
    end 
    vs = [vs, [t*dt;v_ini]];
    
    
    
    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];

%     % draw
%     if mod(t,10)==0
%     axis equal;
%     x_plot =  [x1(1), x2(1), x3(1), x4(1), x5(1), x6(1)];
%     y_plot =  [x1(2), x2(2), x3(2), x4(2), x5(2), x6(2)];
%     z_plot =  [x1(3), x2(3), x3(3), x4(3), x5(3), x6(3)];
%     set (p, 'XData', x_plot, 'YData', y_plot, 'ZData', z_plot);    
%     
%     lines = [];
%     
%     for i=1:6
%         for j=1:6
%             if omega(i,j)>0 % strut
%                 lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)], [x(3,i) x(3,j)], 'linestyle','-','color','r','LineWidth',1.5)];
%             elseif omega(i,j)<0 % cable
%                 lines = [lines; line([x(1,i) x(1,j)], [x(2,i) x(2,j)], [x(3,i) x(3,j)])];
%             end
%         end
%     end 
%     drawnow;
%     delete(lines)
%     hold off;
%     end

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
end

figure(1);
set(figure(1),'Position',[200,100,1550,300]);
subplot(1,3,1);
% raise to zero
v_min = min(vs(2,:));
for i=1:length(vs)
    vs(2,i) = vs(2,i)-v_min;
end
plot(vs(1,:),vs(2,:),'LineWidth',1.5);
hold on;
legend('total potential');
xlabel('t/s'); ylabel('V');

% 
subplot(1,3,2);
plot(x12s(1,:),x12s(2,:),'LineWidth',1.5);
hold on;
plot(x23s(1,:),x23s(2,:),'LineWidth',1.5);
plot(x34s(1,:),x34s(2,:),'LineWidth',1.5);
plot(x41s(1,:),x41s(2,:),'LineWidth',1.5);
legend('cable (1,2)','cable (2,3)','cable (3,4)','cable (4,1)');
xlabel('t/s'); ylabel('║rij║-║rij*║');


subplot(1,3,3);
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
    global expd;
    pull = (x2-x1)*-wij*(norm(x1-x2)/dij)^(2*expd);
end

function push = force_strut(x1,x2,wij,dij)
    push = -1*(x2-x1)/norm(x2-x1)*wij*dij; %constant push
end

function Hessian = cal_hessian(xs)
    x1 = xs(:,1); x2 = xs(:,2); x3 = xs(:,3); x4 = xs(:,4);
%     lijc = 3/4; lijs = 5/4*sqrt(2);
    Hessian = zeros(8,8);

    c12 = eye(2)*0.25*(norm(x1-x2)^2) + 0.5*(x1-x2)*(x1-x2)';
    c23 = eye(2)*0.25*(norm(x2-x3)^2) + 0.5*(x2-x3)*(x2-x3)';
    c34 = eye(2)*0.25*(norm(x3-x4)^2) + 0.5*(x3-x4)*(x3-x4)';
    c41 = eye(2)*0.25*(norm(x1-x4)^2) + 0.5*(x1-x4)*(x1-x4)';

    s13 = eye(2)*-sqrt(2)/4*power(norm(x1-x3)^2, -0.5) + sqrt(2)/4*power(norm(x1-x3)^2, -1.5)*(x1-x3)*(x1-x3)';
    s24 = eye(2)*-sqrt(2)/4*power(norm(x2-x4)^2, -0.5) + sqrt(2)/4*power(norm(x2-x4)^2, -1.5)*(x2-x4)*(x2-x4)';

    Hessian(1:2, 1:2) = c12+s13+c41; Hessian(1:2, 3:4) = -c12; Hessian(1:2, 5:6) = -s13; Hessian(1:2, 7:8) = -c41; 
    Hessian(3:4, 1:2) = -c12; Hessian(3:4, 3:4) = c12+c23+s24; Hessian(3:4, 5:6) = -c23; Hessian(3:4, 7:8) = -s24; 
    Hessian(5:6, 1:2) = -s13; Hessian(5:6, 3:4) = -c23; Hessian(5:6, 5:6) = s13+c23+c34; Hessian(5:6, 7:8) = -c34; 
    Hessian(7:8, 1:2) = -c41; Hessian(7:8, 3:4) = -s24; Hessian(7:8, 5:6) = -c34; Hessian(7:8, 7:8) = c41+s24+c34; 
    
    eigens = eig(Hessian);
    eigens = real(eigens);
    if min(eigens)<-0.1
        xiaoyule = 1
        eigens
    end

end