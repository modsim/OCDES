function demo4b
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example:
%dx=Ax+Bv
%min 1/2 vT*M*v
%s.t. L*z >= Kx

%This example fulfil the sufficient condition of reviewer. It has a
%solution

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
   1 -2];
B=[1 0
   0 1];
M=[5 -2
   -2 5]; 
L=[5 16 
   2 -10];
K=[1 0
   0 2];
f=A*x+B*v; %excellular equation, dx=f(x,v)
obj=v'*M*v/2; %objective function, sym 
heq=[]; %equality constraint, heq=0
hieq=L*v-K*x; % inequality constraints, hieq>=0
x0=[-10;10]; %initial condition of x
%x0=[3.9521;11.1138];
v0=[]; %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=100; %maximum number of updating active set
tstart=0;
tfinal=2.5;
%End: Defined by user

%%Main part
[tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);



% figure
% plot(tout,yout(:,1),'red');
% hold on
% plot(tout,yout(:,2),'blue');
% plot(tout,yout(:,3),'red');
% plot(tout,yout(:,4),'blue');
% legend('x_1(t)','x_2(t)','z_1(t)','z_2(t)')
% xlabel('t');
% set(gca,'FontSize', 14)

figure
plot(tout,yout(:,5),'-','LineWidth',1.5);
hold on
plot(tout,yout(:,6),'--','LineWidth',1.5);

matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
hieq_fun=@(allvar)hieq_gen(allvar);
hieq_v=[];
for i=1:length(tout)
hieq_v(:,i)=hieq_fun(yout(i,:)');
end
h1=plot(tout,hieq_v(1,:),'-.','LineWidth',1.5);
hold on
h2=plot(tout,hieq_v(2,:),':','LineWidth',1.5);
legend('\mu_1 (t)','\mu_2(t)','l_1(t)','l_2(t)')
xlabel('t');
set(gca,'FontSize', 14)



%%Plot simulation results
figure
subplot(2,2,1)
h1=plot(tout,yout(:,1),'red');
hold on
h2=plot(tout,yout(:,2),'blue');
legend('x_1(t)','x_2(t)')
xlabel('t');
set(gca,'FontSize', 14)

subplot(2,2,2)
h1=plot(tout,yout(:,3),'red');
hold on
h2=plot(tout,yout(:,4),'blue');
legend('z_1(t)','z_2(t)')
xlabel('t');
set(gca,'FontSize', 14)

subplot(2,2,3)
h1=plot(tout,yout(:,5),'red');
hold on
h2=plot(tout,yout(:,6),'blue');
legend('\mu_1 (t)','\mu_2(t)')
xlabel('t');
set(gca,'FontSize', 14)


subplot(2,2,4)
matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
hieq_fun=@(allvar)hieq_gen(allvar);
hieq_v=[];
for i=1:length(tout)
hieq_v(:,i)=hieq_fun(yout(i,:)');
end
h1=plot(tout,hieq_v(1,:),'red');
hold on
h2=plot(tout,hieq_v(2,:),'blue');
legend('l_1(t)','l_2(t)')
xlabel('t');
set(gca,'FontSize', 14)


end
