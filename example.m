clc;
clear;
% The small example shows the minimization of the function with multiple minima

% the description of parameters is given in params.m file as well as in the
% original paper

% Guénolé, J., Nöhring, W. G., Vaid, A., Houllé, F., Xie, Z., Prakash, A., & Bitzek, E. (2020). 
% "Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization 
% in atomistic simulations and its implementation in lammps" // Computational Materials Science, 175, 109584.

% initial time step for the integrator
p.dt=1e-4;
% tolerance on the gradient norm to stop optimization
p.tol=1e-9;
% type of the integrator
p.Integrator='verlet';
% save every "p.iter_save" iteration of the algorithm during the minimization
p.iter_save=1;
% display results on every "p.iter_disp" iteration
p.iter_disp=1;
%%% initial guess
x0=[0.1;-0.4];

%%% perform the optimization
out = FIRE(@(x) fun_en(x),x0,p);

%%% calculate function values on a plane
N_grid = 100; %%% make a N_grid x N_grid of points from (-1,1)
[X,Y] = meshgrid(linspace(-1,1,N_grid),linspace(-1,1,N_grid));
X=X(:);
Y=Y(:);
Z=arrayfun(@(N) fun_en([X(N);Y(N)]),(1:length(X))');
X=reshape(X,N_grid,N_grid);
Y=reshape(Y,N_grid,N_grid);
Z=reshape(Z,N_grid,N_grid);

%%% plot the results (minima shown by a red marker)
figure
contour(X,Y,Z)
colorbar
hold on
plot(out.x(1),out.x(2),'ro','markersize',10,'markerfacecolor','r')
hold off

%%% define the function and its gradient
function [f,g] = fun_en(x_val)
    %%% values of two variables
    x = x_val(1);
    y = x_val(2);
    %%% define function values
    f    = exp(x-2*x^2-y^2)*sin(6*(x + y + x*y^2));
    %%% its gradient (column vector)
    g = zeros(length(x_val),1);
    g(1) = cos(6*x*y^2 + 6*y + 6*x)*exp(- 2*x^2 + x - y^2)*(6*y^2 + 6) - sin(6*x*y^2 + 6*y + 6*x)*exp(- 2*x^2 + x - y^2)*(4*x - 1);
    g(2) = cos(6*x*y^2 + 6*y + 6*x)*exp(- 2*x^2 + x - y^2)*(12*x*y + 6) - 2*y*sin(6*x*y^2 + 6*y + 6*x)*exp(- 2*x^2 + x - y^2);
end