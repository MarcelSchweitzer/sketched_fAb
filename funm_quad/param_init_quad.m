function [ param,modified ] = param_init_quad( param )
%PARAM_INIT_QUAD generates/checks the input parameter struct
%  param = PARAM        returns a default setting
%  param = PARAM(param) returns a valid parameter setting
%
%  This file is part of the FUNM_QUAD code described in 
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%


modified = 0;

if ~nargin,
    param = struct;
    modified = 1;
end


if ~isfield(param,'verbose'),
    param.verbose = 1;
    disp('Warning: .verbose not specified, set to 1.');
    modified = 1;
end

if ~isfield(param,'function'),
    param.function = 'invSqrt';
    if param.verbose, 
        disp('Warning: .function not specified, set to invSqrt.');
    end
    modified = 1;
end

if ~isfield(param,'tol'),
    param.tol = 1e-13;
    if param.verbose, 
        disp('Warning: quadrature tolerance .tol not specified, set to 1e-13.');
    end
    modified = 1;
end
    
if ~isfield(param,'thick'),
    param.thick = [];
    if param.verbose, 
        disp('Warning: .thick not specified, set to [].');
    end
    modified = 1;
end

if ~isfield(param,'inner_product'),
    param.inner_product = @(x,y) y'*x;
    if param.verbose, 
        disp('Warning: .inner_product not specified, set to @(x,y) y''*x.');
    end
    modified = 1;
end

if ~isfield(param,'V_full'),
    param.V_full = false;
    if param.verbose, 
        disp('Warning: .V_full not specified, set to false.');
    end
    modified = 1;
end

if ~isfield(param,'H_full'),
    param.H_full = true;
    if param.verbose, 
        disp('Warning: .H_full not specified, set to true.');
    end
    modified = 1;
end

if ~isfield(param,'restart_length'),
    param.restart_length = 50;
    if param.verbose, 
        disp('Warning: .restart_length not specified, set to 50.');
    end
    modified = 1;
end

if ~isfield(param,'max_restarts'),
    param.max_restarts = 10;
    if param.verbose, 
        disp('Warning: .max_restarts not specified, set to 10.');
    end
    modified = 1;
end

if ~isfield(param,'hermitian'),
    param.hermitian = false;
    if param.verbose, 
        disp('Warning: .hermitian not specified, set to false.');
    end
    modified = 1;
end

if isa(param.function,'function_handle') && param.hermitian == false,
    if param.verbose, 
        disp('Warning: passing a function handle as param.function is only recommended for Hermitian matrices due to possible instabilities in the diagonalization.');
    end
    modified = 1;
end


if ~isfield(param,'reorth_number'),
    param.reorth_number = 0;
    if param.verbose, 
        disp('Warning: .reorth_number not specified, set to 0.');
    end
    modified = 1;
end


if ~isfield(param,'exact'),
    param.exact = [];
    if param.verbose, 
        disp('Warning: .exact not specified, no error available.');
    end
    modified = 1;
end

if ~isfield(param,'stopping_accuracy'),
    param.stopping_accuracy = 1e-12;
    if param.verbose, 
        disp('Warning: .stopping_accuracy not specified, set to 1e-12.');
    end
    modified = 1;
end

if ~isfield(param,'min_decay'),
    param.min_decay = 0.95;
    if param.verbose, 
        disp('Warning: .min_decay not specified, set to 0.95.');
    end
    modified = 1;
end

if ~isfield(param,'waitbar'),
    param.waitbar = true;
    if param.verbose, 
        disp('Warning: .waitbar not specified, set to true.');
    end
    modified = 1;
end

if ~isempty(param.thick) && ~isfield(param,'number_thick'),
    param.number_thick = 5;
    if param.verbose, 
        disp('Warning: .number_thick not specified, set to 5.');
    end
    modified = 1;
end

if ~isfield(param,'transformation_parameter') && strcmp(param.function,'invSqrt'),
    param.transformation_parameter = 1;
    if param.verbose, 
        disp('Warning: .transformation_parameter not specified, set to 1.');
    end
    modified = 1;
end

if ~isfield(param,'truncation_length'),
    param.truncation_length = inf;
    if param.verbose, 
        disp('Warning: .truncation_length not specified, set to inf.');
    end
    modified = 1;
end
