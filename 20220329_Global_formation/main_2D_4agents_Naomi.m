%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Time: 2022.3.23
%  Topic:  Naomi 2009 reproduce
%  On 2D plane four agents, for compariasion with Mine
% 
%
%  A PROBLEM: there exists some other equivalent configuration
%
%  Time: 2022.3.24
%  Add Hessian calculation, at specific points
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
omega = D*D';
dijd = [0     norm(q1d - q2d) norm(q1d - q3d) norm(q1d - q4d);
        0               0     norm(q2d - q3d) norm(q2d - q4d);
        0               0             0       norm(q3d - q4d);
        0               0             0                 0   ];
dijd = dijd+dijd';
%%

% x1_ini = q1d+ 10*[rand()-0.5;rand()-0.5];
% x2_ini = q2d+ 10*[rand()-0.5;rand()-0.5];
% x3_ini = q3d+ 10*[rand()-0.5;rand()-0.5];
% x4_ini = q4d+ 10*[rand()-0.5;rand()-0.5];
% 
% xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];

% rot = [cos(0.1) -sin(0.1);sin(0.1) cos(0.1)];
% x1_ini = rot*q1d+[0.1;0.2];
% x2_ini = rot*q2d+[0.1;0.2];
% x3_ini = rot*q3d+[0.1;0.2];
% x4_ini = rot*q4d+[0.1;0.2];
% 
% xs_ini = [x1_ini, x2_ini, x3_ini, x4_ini];

% xs_ini = [3.1472   -2.7301    2.3236   -2.2150;  % for origin pi/atan
%           4.0579    4.1338   -3.0246    1.4688];
%%
xs_ini = [-2.6272    5.6309    1.2114   -0.1110;  % for aij=+-4
          -0.4115    0.4681   -1.6841    2.2406];

% xs_ini = [0    1    1  1.6149;  % 
%           0   0   1    -0.6148];

dt = 0.001;
steps = 1000;

x1 = xs_ini(:,1); x2 =xs_ini(:,2); x3 = xs_ini(:,3); x4 = xs_ini(:,4);

x12s = [0;norm(x1-x2)-dijd(1,2)];x23s = [0;norm(x2-x3)-dijd(2,3)];x34s = [0;norm(x3-x4)-dijd(3,4)];x41s = [0;norm(x4-x1)-dijd(4,1)];

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
                v(:,i) = v(:,i)+force(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            elseif omega(i,j)<0 % cable
                v(:,i) = v(:,i)+force(x(:,i),x(:,j),omega(i,j),dijd(i,j));
            end
        end
    end 
    v = v*20;

    x1 = x1+v(:,1)*dt; x2 = x2+v(:,2)*dt; x3 = x3+v(:,3)*dt; x4 = x4+v(:,4)*dt;

    xs = [x1, x2, x3, x4];
%     Hessian = cal_hessian(xs);

    x1s = [x1s,x1];x2s = [x2s,x2];x3s = [x3s,x3];x4s = [x4s,x4];

    x12s = [x12s, [t*dt;norm(x1-x2)-dijd(1,2)]]; x23s = [x23s, [t*dt;norm(x2-x3)-dijd(2,3)]]; 
    x34s = [x34s, [t*dt;norm(x3-x4)-dijd(3,4)]]; x41s = [x41s, [t*dt;norm(x4-x1)-dijd(4,1)]]; 
end

figure(1);
set(figure(1),'Position',[200,100,1350,450]);

subplot(1,2,1);
plot(x12s(1,:),x12s(2,:),'LineWidth',1.5);
hold on;
plot(x23s(1,:),x23s(2,:),'LineWidth',1.5);
plot(x34s(1,:),x34s(2,:),'LineWidth',1.5);
plot(x41s(1,:),x41s(2,:),'LineWidth',1.5);
legend('cable (1,2)','cable (2,3)','cable (3,4)','cable (4,1)');
xlabel('t/s'); ylabel('║rij║-║rij*║');

subplot(1,2,2);
xs_p=zeros(4,2,steps+1);
xs_p(1,:,:) = x1s; xs_p(2,:,:) = x2s; xs_p(3,:,:) = x3s; xs_p(4,:,:) = x4s;
for i=1:4
    for j = 1:steps+1
        px(j) = xs_p(i,1,j);
        py(j) = xs_p(i,2,j);
        
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
                    pi(k) = xs_p(i,k,steps+1);pj(k) = xs_p(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)], 'linestyle','-','color','r','LineWidth',1.5)];
            elseif omega(i,j)<0 % cable
                for k=1:2
                    pi(k) = xs_p(i,k,steps+1);pj(k) = xs_p(j,k,steps+1);
                end
                lines = [lines; line([pi(1) pj(1)], [pi(2) pj(2)],'LineWidth',1.5)];
            end
        end
    end 
    axis equal;
    xlabel('x');ylabel('y');zlabel('z');
    legend([H([1 2 3 4]), G, K] ,'agent 1','agent 2','agent 3','agent 4','initial configuration','final configuration');
    
%     H =[(8*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(1/2) + (2*(2*xx - 2)^2)/((xx - 1)^2 + (yy - 1)^2) - (8*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(1/2) + (8*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(1/2) + (8*xx^2)/(xx^2 + yy^2) + (2*(2*xx - 2)^2)/((xx - 1)^2 + yy^2) + (2*(2*xx - 2)^2*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(3/2) - (8*xx^2*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(3/2) - (2*(2*xx - 2)^2*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(3/2),                                                                                                                                                                           (8*xx*yy)/(xx^2 + yy^2) + (4*yy*(2*xx - 2))/((xx - 1)^2 + yy^2) + (2*(2*xx - 2)*(2*yy - 2))/((xx - 1)^2 + (yy - 1)^2) + (4*yy*(2*xx - 2)*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(3/2) - (2*(2*xx - 2)*(2*yy - 2)*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(3/2) - (8*xx*yy*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(3/2);(8*xx*yy)/(xx^2 + yy^2) + (4*yy*(2*xx - 2))/((xx - 1)^2 + yy^2) + (2*(2*xx - 2)*(2*yy - 2))/((xx - 1)^2 + (yy - 1)^2) + (4*yy*(2*xx - 2)*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(3/2) - (2*(2*xx - 2)*(2*yy - 2)*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(3/2) - (8*xx*yy*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(3/2), (8*yy^2)/((xx - 1)^2 + yy^2) + (8*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(1/2) + (2*(2*yy - 2)^2)/((xx - 1)^2 + (yy - 1)^2) - (8*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(1/2) + (8*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(1/2) + (8*yy^2)/(xx^2 + yy^2) + (8*yy^2*((5*2^(1/2))/4 - ((xx - 1)^2 + yy^2)^(1/2)))/((xx - 1)^2 + yy^2)^(3/2) - (8*yy^2*((xx^2 + yy^2)^(1/2) - 3/4))/(xx^2 + yy^2)^(3/2) - (2*(2*yy - 2)^2*(((xx - 1)^2 + (yy - 1)^2)^(1/2) - 3/4))/((xx - 1)^2 + (yy - 1)^2)^(3/2)];


%% functions
function f = force(x1,x2,wij,dij)
%     aij = pi/atan(-wij);
%     lij = dij*(1-1/aij);
    if wij<0 %cable
        aij = 4;
        lij = dij*(1-1/aij);
    else % strut
        aij = -4;
        lij = dij*(1-1/aij);
    end
    f = aij*(-wij)*((x2-x1)/1)*(1-lij/norm(x2-x1));
end

function Hessian = cal_hessian(xs)
    x1 = xs(:,1); x2 = xs(:,2); x3 = xs(:,3); x4 = xs(:,4);
    lijc = 3/4*0.95710678118654790580421831691638; lijs = 1.25*1.3535533905932739529021091584582;
%     lijc = 3/4; lijs = 5/4*sqrt(2);
    Hessian = zeros(8,8);

    c12 = eye(2)*(1-lijc*power(norm(x1-x2)^2, -0.5)) + lijc*power(norm(x1-x2)^2, -1.5)*(x1-x2)*(x1-x2)';
    c23 = eye(2)*(1-lijc*power(norm(x2-x3)^2, -0.5)) + lijc*power(norm(x2-x3)^2, -1.5)*(x2-x3)*(x2-x3)';
    c34 = eye(2)*(1-lijc*power(norm(x3-x4)^2, -0.5)) + lijc*power(norm(x3-x4)^2, -1.5)*(x3-x4)*(x3-x4)';
    c41 = eye(2)*(1-lijc*power(norm(x4-x1)^2, -0.5)) + lijc*power(norm(x4-x1)^2, -1.5)*(x4-x1)*(x4-x1)';

    s13 = eye(2)*(1-lijs*power(norm(x1-x3)^2, -0.5)) + lijs*power(norm(x1-x3)^2, -1.5)*(x1-x3)*(x1-x3)';
    s24 = eye(2)*(1-lijs*power(norm(x2-x4)^2, -0.5)) + lijs*power(norm(x2-x4)^2, -1.5)*(x2-x4)*(x2-x4)';

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

% ddp = [-xs(2,1);xs(1,1);-xs(2,2);xs(1,2);-xs(2,3);xs(1,3);-xs(2,4);xs(1,4);]

 