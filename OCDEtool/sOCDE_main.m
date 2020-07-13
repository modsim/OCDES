function [tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f, obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal,opt_init,opt_sol)
%This function solves OCDE by calling OCDE_sol_init.m for initialization and
%OCDE_sol.m for integration

%Input:
%f: state equations of the dynamic part (sym)
%obj: objective function (sym)
%heq: left hand side of equality constraints (sym)
%hieq: left hand side of inequality (>=) constraints (sym)
%x:  state varaibles(sym)
%x0: initial point of x (real column vector). x0 must be given and fixed.
%v:  inner optimization variables (sym)
%v0: initial guess of inner optimization variables (real column vector). v0 can be empty.
%opt_init: options for initilization
%opt_sol:  options for integration

%Outputs:
%tout: integration steps of t (real vector)
%yout: integration results of [x;v;miu_ieq;miu_eq] (real matrix)
%teout: the time, when the switching of active set happens (real vector)
%yeout: states of [x;v;miu_ieq;miu_eq], when the switching of active set happens (real matrix)
%ieq_actout: record the different values of active sets along the solution trajectory (real matrix)
%miu_ieq: L. multipliers for inequality constraints
%miu_eq:  L. multipliers for   equality constraints

[KKT,miu_ieq, miu_eq, ieq_act,KKT_fun,event_fun, init]=OCDE_sol_init...
    (obj,heq, hieq,x,x0,v,v0,opt_init);

v0=eval(subs(v,init{:}));
miu_ieq0=eval(subs(miu_ieq,init{:}));
miu_eq0=eval(subs(miu_eq,init{:}));
ieq_act0=eval(subs(ieq_act,init{:}));
Jac_fun=[];

matlabFunction(f,'File','f_gen','Vars',{x,v});
f_fun=@(x,v)f_gen(x,v);

[tout,yout,teout,yeout,ieq_actout]=OCDE_sol(f_fun, KKT_fun, x0, v0, miu_ieq0, miu_eq0,...
    ieq_act0, tstart,tfinal,event_fun, Jac_fun,opt_sol);

end
