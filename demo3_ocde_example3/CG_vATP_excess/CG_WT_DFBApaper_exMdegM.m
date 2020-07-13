function In=CG_WT_DFBApaper_exMdegM(S,O)       
In.exM_ch={'CO2ex','NH3ex','SO3ex','O2ex','H_ex','GLCex','Biomass'};
In.degM_ch={'ADP','NAD','NADP','COA','MQ'}; %so that A_in has full rank


%init conditions
x_s=O.x;
c_GLC=x_s(2);

In.v_upt_ch={'GLC_t_PEP'};
In.v_upt_tom=[sym(In.v_upt_ch{1})];

init_x0={O.x,O.x0};
% In.v_upt_rho=[vmax_g*S1/(S1+K_g);vmax_x*S2/(S2+K_x)];
rou_A=4.5*c_GLC/(c_GLC+1);

In.v_upt_rho=[rou_A];
In.v_upt_rho0=subs(In.v_upt_rho,init_x0{:});


%is v_upt eq?
In.v_upt_iseq=[1];

%obj fun
v_ch=S.FluxNames;
for i=1:length(v_ch)
v_all{i} =sym(v_ch{i});  
end
v_all=vertcat(v_all{:});
In.v_in=v_all;

%OBG1
bio=sym('bmsynth');
v_excess=sym('ATP_excess');

In.obj=-v_excess/(v_all'*v_all);

%OBG2
%In.obj=-v_all'*v_all;

%OBG3
%In.obj=In.v_upt_tom/bio; %No switch

% %OBG4 
% v_sucD=sym('sucD');
% v_pyk=sym('pyk');
% v_pgk=sym('pgk');
% v_ATPase=sym('ATPase');
% v_atp=v_sucD+v_pyk+v_pgk+v_ATPase;
% In.obj=-v_atp/(v_all'*v_all); %No sitch
end