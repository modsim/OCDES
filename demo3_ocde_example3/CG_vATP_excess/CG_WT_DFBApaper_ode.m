function O=CG_WT_DFBApaper_ode()
f_fun=@(x,v_in)f(x,v_in);
syms Biomass GLCex 
x=[Biomass;GLCex];
x0=[1;20];

O.f_fun=f_fun;
O.x=x;
O.x0=x0;

end

function dx=f(x,v_in)
v_upt=v_in(38);
v_biomass=v_in(33);

X=x(1);

dx=[v_biomass*X;
    -v_upt*X;   
    ];
end