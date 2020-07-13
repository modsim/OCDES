% function hessen=hessen_FD(f_fun, x, eps)
% %dfdy by using stepwise finite difference
% y=f_fun(x);
% 
% m=length(y);
% n=length(x);
% 
% if m>=2
%    dfdx=[];
% end
% 
% dfdx_fun=@(x)gradient_FD(f_fun, x, eps);
% hessen=gradient_FD(dfdx_fun, x, eps);
% 
%
% end