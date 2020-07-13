clear;
clc;
currentFolder=pwd;
cd ..\
cd ..\
addpath(strcat(pwd,'\OCDEtool'));
addpath(strcat(pwd,'\DFBA_addon'));
cd(currentFolder);


%addpath('C:\Users\Xiao\sciebo\MATLAB\mo_gen_TL_Xiaosave');

S = CG_WT_OCDEpaper();
O = CG_WT_DFBApaper_ode();
In =CG_WT_DFBApaper_exMdegM(S,O);
A_tot=S.Stoich; %dx(all meta)=Av
x_ch=S.PoolNames;
v_ch=S.FluxNames;
v_r=S.Reversibilities;


[exM_p,inM_p]=chpos(x_ch, [In.exM_ch,In.degM_ch]);
A_in=A_tot(inM_p,:);
x_ch_in=x_ch(inM_p);

init_x0={O.x,O.x0};
nv=length(v_ch);
nx=length(x_ch_in);

%remove equality of Sv=0 
Aeq=A_in; 
beq=zeros(length(x_ch_in),1);    
[T,B,n_idx,r_idx]=substitutex1byx2(Aeq,beq); %so that v=T*v2+B

% if any(abs(B))>=1e-4
%     fprintf('B \neq 0, some later results may be wrong. \n');
% end

%treat with v_upt
%v_upt_0=v_max_0*C_glu0/(C_glu0+k_upt_0);
v_upt_pos=chpos(v_ch, In.v_upt_ch);
n_sub=length(v_upt_pos); %sustrate number
n_sub_eq=sum(In.v_upt_iseq); %substrate with equality uptake rate
n_sub_less=n_sub-n_sub_eq; %substrate with inequality uptake rate

A_upt=zeros(n_sub,nv); %A_upt*v_upt (= or <=) v_upt0
for i=1:n_sub
A_upt(i,v_upt_pos(i))=1;
end
Aupt_eq=[]; %rows in A_upt, w.r.t =
bupt_eq=[];
Aupt_ieq=[];%rows in A_upt, w.r.t <=
bupt_ieq=[];
for i=1:length(In.v_upt_iseq)
    if In.v_upt_iseq(i)==1
       Aupt_eq=[Aupt_eq;A_upt(i,:)];
       bupt_eq=[bupt_eq;In.v_upt_rho(i)];      
    else
       Aupt_ieq=[Aupt_ieq;A_upt(i,:)];
       bupt_ieq=[bupt_ieq;In.v_upt_rho(i)];
    end
end


%Airr, birr
nieq=0;
for i=1:length(v_ch)
    if v_r(i)==0 && all(i~=v_upt_pos) %not irreversable && not fixed  
       nieq=nieq+1;       
    end    
end
Airr=zeros(nieq,size(A_in,2));
irow=1;
for i=1:length(v_ch)
    if v_r(i)==0 && all(i~=v_upt_pos) %not irreversable && not fixed  
        Airr(irow,i)=-1;
        irow=irow+1;
    end
    
end
birr=zeros(nieq,1);

% %DFBA test
% init_x0={O.x,O.x0};
% heq=[Aeq*In.v_in;
%      Aupt_eq*In.v_in-bupt_eq];
% hieq=[-Airr*In.v_in+birr];
% v0=rand(length(In.v_in),1);
% opt_opt=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5000, 'MaxFunctionEvaluations',50000);
% y=getLagrMulti(subs(In.obj,init_x0{:}), subs(heq,init_x0{:}), subs(hieq,init_x0{:}), In.v_in, v0,opt_opt);
% obj_m=subs(In.obj,In.v_in,y.y_opt);
% %DFBA test

% %solve by EF approach
% opt_opt=optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',5000, 'MaxFunctionEvaluations',50000);
% %opt_sol=odeset('RelTol', 1e-6, 'AbsTol', 1e-6,'InitialStep',1e-6); %option for integrator
% opt_sol=[];
% x=O.x;
% x0=O.x0;
% v=In.v_in;
% v0=rand(length(v),1);
% f=O.f_fun(x,v);
% obj=In.obj;
% heq=[
% A_in*v
% Aupt_eq*v-In.v_upt_rho
% ];
% hieq=[-Airr*v+birr];
% tstart=0;
% tfinal=3;
% [tout_EF,yout_EF]=sOCDE_EF_main(f, In.obj, heq, hieq,x,x0,v,v0,tstart, tfinal,opt_opt,opt_sol);
% %solve by EF approach


%%%%%%%%========using extreme rays=================
%%%%%%LP after replacement%%%%%
v2=In.v_in(n_idx);
%calculate extreme rays by 0=Sv, virr>=0
% v=ext'*y with y>=0
Aieq_n=Airr*T;
ine=-Aieq_n;
equ=[];
% options=[];
% zerotol=1e-10;
%[ext, bas, edges, totalnumrays, totalnumedges] = skeleton(ine, equ, options, zerotol);
% save('ext','ext');
load ext
global scale_mu_
scale_mu_=1;

v2_fun=@(y)ext'*y;
v_in_fun=@(y)(T*v2_fun(y)+B);

%assume that B=0
sizey=size(ext,1);
y_tom=sym('y',[sizey,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=O.x;
z=[In.v_in;y_tom];
f=O.f_fun(x,T*(ext'*y_tom)+B); %excellular equation, dx=f(x,v)
obj=In.obj+1e-3*y_tom'*y_tom; %objective function, sym 
heq=[In.v_in-T*ext'*y_tom;     
     Aupt_eq*In.v_in-bupt_eq]; %equality constraint, heq=0
hieq=y_tom; % inequality constraints, hieq>=0
x0=O.x0; %initial condition of x
z0=ones(length(z),1); %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5000, 'MaxFunctionEvaluations',50000);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14,'MaxStep',1e-3); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=100; %maximum number of updating active set
tstart=0;
tfinal=3;

tic
[tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(...
    f,obj,heq, hieq,x,x0,z,z0,tstart, tfinal, opt_init,opt_sol);
teout=[teout;tout(end)];
yeout=[yeout;yout(end,:)];
ieq_actout=[ieq_actout,ieq_actout(:,end)];
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
plot(tout,yout(:,1),'--','LineWidth',1); 
hold on
plot(tout,yout(:,2)/5,'-','LineWidth',1); 
set(gca,'fontsize',14)
grid on; 
xlabel('Time [h]', 'FontSize', 14);
ylabel('Concentration', 'FontSize', 14);
axis([0 tfinal -0.1 5])    
legend('c_{Biomass}(t)','c_{GLC}(t)/5','Location','northeast','Orientation','vertical')
fig=gcf;
saveas(fig,'CG_solu_KKT','epsc')

figure
N=length(teout);
Nieq=sum(ieq_actout);
plot(teout(2:end-1),Nieq(2:end-1)-1,'o','LineWidth',1); 
hold on
for i=1:N-1
    Nieq_i=Nieq(i);
    ts=[teout(i):1e-2:teout(i+1)];    
    plot(ts,Nieq_i*ones(length(ts),1),'-','LineWidth',2);
    hold on
end
axis([0 tfinal 90 100])    
xlabel('time');
ylabel('No. of active inequality const.')
set(gca,'fontsize',14);
fig=gcf;
saveas(fig,'CG_activeset','epsc')

ActIeq=[11,21,27,81,100];
ActIeq_scale=40;
ActIeq_all=[11,21,27,40,81,100];
figure
plot(tout,yout(:,length([O.x0;In.v_in])+ActIeq),'LineWidth',1); 
hold on
plot(tout,yout(:,length([O.x0;In.v_in])+ActIeq_scale)/10,'-','LineWidth',1); 
hold on
plot(teout(2:end-1),zeros(length(teout)-2),'o');
set(gca,'fontsize',14)
grid off; 
xlabel('Time', 'FontSize', 14);
%axis([2.5 3 -0.1 2])    
legend('a_{11}(t)','a_{21}(t)','a_{27}(t)','a_{81}(t)','a_{100}(t)','a_{40}(t)/10','Location','northeast','Orientation','vertical')
fig=gcf;
saveas(fig,'CG_ieq','epsc')

figure
plot(tout,yout(:,length([x;z])+ActIeq_all),'LineWidth',1);
hold on
plot(teout(2:end-1),zeros(length(teout)-2),'o');
set(gca,'fontsize',14)
grid off; 
xlabel('Time', 'FontSize', 14);
axis([0 3 -0.0001 0.0002])    
legend('\mu_{11}(t)','\mu_{21}(t)','\mu_{27}(t)','\mu_{40}(t)','\mu_{81}(t)','\mu_{100}(t)','Location','northwest','Orientation','vertical')
fig=gcf;
saveas(fig,'CG_mu','epsc');


%Matrix N
data=A_tot;      %Sample 2-dimensional data
col_header=v_ch;     %Row cell array (for column labels)
row_header=x_ch';     %Column cell array (for row labels)
xlswrite('Matrix_N.xlsx',data,'Tabelle1','B2');     %Write data
xlswrite('Matrix_N.xlsx',col_header,'Tabelle1','B1');     %Write column header
xlswrite('Matrix_N.xlsx',row_header,'Tabelle1','A2');      %Write row header

%Matrix B
data=-Airr;      
col_header=v_ch;    
xlswrite('Matrix_B.xlsx',data,'Tabelle1','A2');    
xlswrite('Matrix_B.xlsx',col_header,'Tabelle1','A1');    

%Matrix Ea
for i=1:length(y_tom)
    a_ch{i}=strcat('a',num2str(i));
    
end
data=T*ext';      %Sample 2-dimensional data
col_header=a_ch;     %Row cell array (for column labels)
row_header=v_ch';     %Column cell array (for row labels)
xlswrite('Matrix_Ea.xlsx',data,'Tabelle1','B2');     %Write data
xlswrite('Matrix_Ea.xlsx',col_header,'Tabelle1','B1');     %Write column header
xlswrite('Matrix_Ea.xlsx',row_header,'Tabelle1','A2');      %Write row header





