function [KKT,miu_ieq, miu_eq, ieq_act,KKT_fun,event_fun,init]=OCDE_sol_init(obj,heq, hieq,x,x0,v,v0,opt)
%This function: 
%(1) formulates the KKT optimality condition of the inner NLP 
%in Matlab symbolic format 
%(2) generates m-files of KKT optimality condition and associated function handles 
%(3) provides initial conditions for state variables, inner optimization 
%variables, Lagrange multipliers, initial active set of inequality constraints, 
%(4) define event functions for detecting the switching of active set.
 
%Inputs:
%obj: objective function (sym)
%heq: left hand side of equality constraints (sym)
%hieq: left hand side of inequality (>=) constraints (sym)
%x:  state varaibles(sym)
%x0: initial point of x (real column vector). x0 must be given and fixed.
%v:  inner optimization variables (sym)
%v0: initial guess of inner optimization variables (real column vector). v0 can be empty.
%opt: options of initilizations
%opt.MaxTry: the number of multi-start times when solving NLP for initialization (integer), default 10.
%opt.tol_act: tolerance of active inequality (real scalar), default 1e-6.
%opt.optimoptions=: options for used optimizer. default is empty.
 
 
%Outputs:
%KKT: KKT optimality condition of inner optimization (sym)
%miu_ieq: Lagrange multipliers of inequality constraints (sym)
%miu_eq: Lagrange multipliers of equality constraints (sym)
%ieq_act: index for active inequality constraint (sym). If ieq_act(j)=1: the j-th inequality constraint is active, otherwise non-active.
%KKT_fun([x;v;miu_ieq;miu_eq],ieq_act): function handle, which
%returns the evaluated value of KKT
%event_fun: event function for detecting the switching of active inequalities
%init: initial values of state variables, inner optimization variables, 
%Lagrange multipliers and the active set


if isfield(opt, 'tol_act') && ~isempty(opt.tol_act)
tol_act=opt.tol_act;
else
tol_act=1e-6;
end  

%KKT conditions
[KKT, miu_ieq, miu_eq, ieq_act,He]= KKTcon(obj,v,heq,hieq);
%matlabFunction(KKT,'File','KKT_gen','Vars',{[x;v;miu_ieq;miu_eq],ieq_act});
matlabFunction(KKT,'File','KKT_gen','Optimize',false,'Sparse',true,'Vars',{[x;v;miu_ieq;miu_eq],ieq_act});
KKT_fun=@(allvar,ieq_act)KKT_gen(allvar,ieq_act);

%event function
matlabFunction(miu_ieq-hieq,'File','event_gen','Vars',{[x;v;miu_ieq;miu_eq]});
event_fun=@(allvar)event_gen(allvar);

% %Hessian function
% matlabFunction(He,'File','Hessian_gen','Vars',{[x;v;miu_ieq;miu_eq]});
% He_fun=@(allvar)Hessian_gen(allvar);

init_x0={x,x0};

if isfield(opt, 'optimoptions') && ~isempty(opt.optimoptions)
opt_opt=opt.optimoptions;
else
opt_opt=[];
end  

y=getLagrMulti(subs(obj,init_x0{:}), subs(heq,init_x0{:}), subs(hieq,init_x0{:}), v, v0,opt_opt);
v0_new=y.y_opt;
init_v0={v,v0_new};


vhieq=subs(subs(hieq,init_x0{:}),init_v0{:});
ieq_act0=zeros(length(vhieq),1);
for i=1:length(vhieq)
    if vhieq(i)<=tol_act
    ieq_act0(i)=1;
    end
end
miu_ieq0=y.lam.hieq;
miu_eq0=y.lam.heq;

init={[x;v;miu_ieq;miu_eq;ieq_act], [x0;v0_new;miu_ieq0;miu_eq0;ieq_act0]};
      
end

%=============================================================================
%  OCDES – Solving Optimization-Constrained Dynamic Systems by optimality
%         conditions
%  
%  Copyright (c) 2020: Xiao Zhao, Forschungszentrum Juelich GmbH, Juelich,
%                      Germany. 
%
%  This code can only be used for academic purposes.                               
%  All rights reserved.
%=============================================================================
