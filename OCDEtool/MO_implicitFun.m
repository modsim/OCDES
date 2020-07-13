function varargout=MO_implicitFun(fun_cell,varargin)
%fun_cell:  a cell of function handles
%varargout=[fun_cell(1)(varargin),...,fun_cell(N)(varargin)]

N=length(fun_cell);
varargout=cell(N,1);
for i=1:N
    if ~isempty(fun_cell{i})
        varargout{i}=fun_cell{i}(varargin{:});
    else
        varargout{i}=[];
    end
end

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