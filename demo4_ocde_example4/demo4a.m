function demo4
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example:
%dx=Ax+Bv
%min 1/2 vT*M*v
%s.t. L*z >= Kx

%Show the situation of CASE D in Fig. 2 of OCDE paper

clear;
clc;
close all;
currentFolder=pwd;
cd ..\
addpath(strcat(pwd,'\OCDEtool'));
cd(currentFolder);

%%Start: Defined by user
syms x1 x2  v1 v2
x=[x1;x2]; %state of ODE
v=[v1;v2]; % var. for inner opt.

A=[-1 1 
   0 -2];
B=[1 0
   0 1];
M=[1 0
   0 -1]; 
L=[1 1 
   1 0];
K=[1 0
   0 2];
f=A*x+B*v; %excellular equation, dx=f(x,v)
obj=v'*M*v/2; %objective function, sym 
heq=[]; %equality constraint, heq=0
hieq=L*v-K*x; % inequality constraints, hieq>=0
x0=[0.5;1]; %initial condition of x
%x0=[3.9521;11.1138];
v0=[]; %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=1; %maximum number of updating active set
tstart=0;
tfinal=20;
%End: Defined by user

%%Main part
[tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);

%backward integartion starting from a point where SCC fail
x0=yout(end,1:2)';
v0=yout(end,3:4)';
tstart=tout(end)';
tfinal=0.01;
[tout_b,yout_b,teout_b,yeout_b,ieq_actout_b,miu_ieq_b, miu_eq_b]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);


figure
plot(yout(:,1),yout(:,2),'blue');
hold on
plot(yout_b(:,1),yout_b(:,2),'red');
xlabel('x_1(t)');
ylabel('x_2(t)');
set(gca,'FontSize', 14)

tout=[tout;tout_b];
yout=[yout;yout_b];
%%Plot simulation results
% figure
% h1=plot(tout,yout(:,1),'red');
% hold on
% h2=plot(tout,yout(:,2),'blue');
% legend('x_1(t)','x_2(t)')
% xlabel('t');
% set(gca,'FontSize', 14)

figure
plot(yout(:,1),yout(:,2),'blue');
xlabel('x_1(t)');
ylabel('x_2(t)');
set(gca,'FontSize', 14)





figure
plot(tout,yout(:,4),'blue');
xlabel('z_1(t)');
ylabel('z_2(t)');
set(gca,'FontSize', 14)

% 
% figure
% h1=plot(tout,yout(:,3),'red');
% hold on
% h2=plot(tout,yout(:,4),'blue');
% legend('v_1(t)','v_2(t)')
% xlabel('t');
% plot(teout(1:end-1),yeout(1:end-1,3),'o');
% plot(teout(1:end-1),yeout(1:end-1,4),'o');

% figure
% h1=plot(tout,yout(:,5),'red');
% hold on
% %h2=plot(tout,yout(:,6),'blue');
% %legend('\mu_1 (t)','\mu_2(t)')
% 
% matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
% hieq_fun=@(allvar)hieq_gen(allvar);
% hieq_v=[];
% for i=1:length(tout)
% hieq_v(:,i)=hieq_fun(yout(i,:)');
% end
% h1=plot(tout,hieq_v(1,:),'blue');
% legend('\mu_1 (t)','hieq(1)(t)')
% xlabel('t');
% 
% figure
% h2=plot(tout,hieq_v(2,:),'blue');
% legend('hieq(2)(t)')
% xlabel('t');


end
