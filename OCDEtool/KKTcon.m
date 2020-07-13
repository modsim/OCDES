function [KKT_con,miu_ieq,miu_eq,ieq_act,He]= KKTcon(obj_sym,x,heq,hieq)
%  min_x obj_sym(x,p) s.t. heq(x,p)=0, hieq(x,p)>=0
%  return KKT condition
%  DL=0
%  heq=0
%  ieq_act_j*heq_j+(1-ieq_act_j)*miu_ieq=0, 
%  where ieq_act \in{0,1} denotes whether%  heq is active


scale_mu_=1;

n_ieq=length(hieq);
n_eq=length(heq);
ieq_act=sym('ieq_act_',[n_ieq,1]);

if n_eq>0
miu_eq =sym('miu_eq_', [n_eq, 1]); 
L=obj_sym+scale_mu_*miu_eq'*heq;
else
miu_eq=[];   
L=obj_sym;
end

if n_ieq>0
miu_ieq=sym('miu_ieq_',[n_ieq, 1]);
L=L-scale_mu_*miu_ieq'*hieq;
else
miu_ieq=[];
end

DxL=jacobian(L,x);
He=hessian(L,x); 
%save('He','He');%test

diag_ieq=diag(ieq_act);
cc=diag_ieq*hieq-(eye(n_ieq,n_ieq)-diag_ieq)*miu_ieq;

KKT_con=[
DxL';
heq;
cc;
];

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