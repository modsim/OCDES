function dfdx=gradient_FD(f_fun, x, eps)
%dfdy by using stepwise finite difference
y=f_fun(x);

m=length(y);
n=length(x);

dfdx=zeros(m,n);
for j=1:n
    xr=x;
    xr(j)=xr(j)+eps;
    yr=f_fun(xr);
    
    xl=x;
    xl(j)=xl(j)-eps;
    yl=f_fun(xl);           
    dfdx(:,j)=(yr-yl)/(2*eps);          
end

end