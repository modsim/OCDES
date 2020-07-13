function demo4c
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example: Extendedn Spirallus example


clear;
clc;
close all;
currentFolder=pwd;
cd ..\
addpath(strcat(pwd,'\OCDEtool'));
cd(currentFolder);

%%Start: Defined by user
syms cx cs1 cs2 cp 
syms vupt1 vupt2 v1 v2 v3 v4 v5 v6 
x=[cx;cs1;cs2;cp]; %state of ODE
v=[vupt1;vupt2;v1;v2;v3;v4;v5;v6]; % var. for inner opt.

Y_DX=1;
f=[
v5/Y_DX*cx;
-vupt1*cx;
-vupt2*cx;
v6*cx;
]; %excellular equation, dx=f(x,v)
obj=sqrt(v'*v); %objective function, sym 
%obj=-v5/sqrt(v'*v); %objective function, sym 
N=[1 0 -1 0 0 -1 0 0;
   0 1 1 -1 0 0 0 0;
   0 0 0 1 -1 0 0 0;
   0 0 -1 0 1 1 -1 0;
   0 0 0 1 1 0 0 -1];
heq=[N*v;
    vupt1-3.8*cs1/(cs1+1);
    vupt2-1.1*cs2/(cs2+1);
    ]; %equality constraint, heq=0
hieq=[v1;v4;v5]; % inequality constraints, hieq>=0

x0=[1;20;10;0]; %initial condition of x
v0=[1 1 1 1 1 1 1 1]'; %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=100; %maximum number of updating active set
tstart=0;
tfinal=1.0;
%End: Defined by user

%%Main part
[tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);


%%Plot simulation results
figure
subplot(2,1,1)
h1=plot(tout,yout(:,1),'red','LineWidth',1.5);
hold on
h2=plot(tout,yout(:,2),'blue','LineWidth',1.5);
hold on
h3=plot(tout,yout(:,3),'black','LineWidth',1.5);
hold on
h4=plot(tout,yout(:,4),'green','LineWidth',1.5);
legend('c_x(t)','c_{s1}(t)','c_{s2}(t)','c_p(t)')
xlabel('t');
set(gca,'FontSize', 12)


subplot(2,1,2)
h1=plot(tout,yout(:,5),'red','LineWidth',1.5);
hold on
h2=plot(tout,yout(:,6),'blue','LineWidth',1.5);
legend('v_{upt1}(t)','v_{upt2}(t)')
xlabel('t');
set(gca,'FontSize', 12)

figure
subplot(2,1,1)
h1=plot(tout,yout(:,13),'red','LineWidth',1.5);
hold on
h2=plot(tout,yout(:,14),'blue','LineWidth',1.5);
hold on
h2=plot(tout,yout(:,15),'green','LineWidth',1.5);
legend('\mu_1 (t)','\mu_2(t)','\mu_3(t)')
xlabel('t');
set(gca,'FontSize', 12)
ylim([-0.1 2])

subplot(2,1,2)
matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
hieq_fun=@(allvar)hieq_gen(allvar);
hieq_v=[];
for i=1:length(tout)
hieq_v(:,i)=hieq_fun(yout(i,:)');
end
h1=plot(tout,hieq_v(1,:),'red','LineWidth',1.5);
hold on
h2=plot(tout,hieq_v(2,:),'blue','LineWidth',1.5);
hold on
h2=plot(tout,hieq_v(3,:),'green','LineWidth',1.5);
legend('l_1(t)','l_2(t)','l_3(t)')
xlabel('t');
set(gca,'FontSize', 12)
ylim([-0.1 4.5])
end
