function out = FIRE(FUN,x0,varargin)
%%% optimization with the fast inertial relaxation engine (FIRE) algorithm following the publication
% Guénolé, J., Nöhring, W. G., Vaid, A., Houllé, F., Xie, Z., Prakash, A., & Bitzek, E. (2020). 
% "Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization 
% in atomistic simulations and its implementation in lammps" // Computational Materials Science, 175, 109584.


%%% FUN is the energy function and its gradient
%%% example: function [e,g] = my_energy(x)
%%%          e=x^2; % energy
%%%          g=2*x % respective gradient (force f=-g)
%%%          g should be a column vector

%%% Initialise positions,speeds,energy and forces
xk = x0; %%% should be a column vector
vk = zeros(size(x0));
[~,gk] = feval(FUN,xk);
Fk=-gk;

% Create parser
par = inputParser;

% Set parameters
par = params(par);

% Parse input
par.parse(varargin{:});


%%% initialise parameters (see pararms.m)
alpha0=par.Results.alpha; %%% initial value of alpha
%%% number of steps
N_pos=0; %%% since P=v*F>0
N_neg=0; %%% since P=v*F<=0
%%% counter
i=1;

%%% set the current parameter structure
par_cur=par.Results;
%%% solution inside iteration
x_iter=[];

%%% check if the initial guess satisfies the minimization
if norm(Fk)<par.Results.tol
    %%% save the solution at current iteration
    if par.Results.iter_save>0
        if mod(i,par.Results.iter_save)==0
            x_iter=[x_iter xk];
        end
    end
    if par.Results.iter_disp>0
        if mod(i,par.Results.iter_disp)==0
            display(['Internal iterations ',num2str(i),', norm grad ', num2str(norm(norm(Fk)))]);
        end
    end
else
    while i<=par.Results.IterMax
        %%% calculate the power
        P=vk'*Fk;
        if P>0
            %%% increment of N_pos
            N_pos=N_pos+1;
            %%% drop N_neg
            N_neg=0;
            if N_pos>par.Results.Nmin
                %%% update current parameters
                par_cur.dt=min(par.Results.dt*par.Results.dt_inc,par.Results.dt_max);
                par_cur.alpha=par.Results.alpha*par.Results.f_alpha;
            end
        else
            %%% drop N_pos
            N_pos=0;
            %%% increment of N_neg
            N_neg=N_neg+1;
            %%% if P<=0 for NP_max, quit
            if N_neg>par.Results.NP_max
                break
            end
            %%% this is to adjust dt and alpha
            if ~(par.Results.initial_delay && i<par.Results.Nmin)
                if par.Results.dt*par.Results.dt_dec>par.Results.dt_min
                    %%% update current parameters
                    par_cur.dt=par.Results.dt*par.Results.dt_dec;
                end
                %%% update current parameters
                par_cur.alpha=alpha0;
            end
            %%% avoid uphill motion
            xk=xk-0.5*par_cur.dt*vk;
            vk=zeros(size(vk));
        end
        % Set parameters
        par.parse(par_cur);
        %%% integrate the equations
        [xk,vk,Fk] = MDintegrator(FUN,xk,vk,Fk,par.Results);
        %%% save the solution at current iteration
        if par.Results.iter_save>0
            if mod(i,par.Results.iter_save)==0
                x_iter=[x_iter xk];
            end
        end
        if par.Results.iter_disp>0
            if mod(i,par.Results.iter_disp)==0
                display(['Internal iterations ',num2str(i),', norm grad ', num2str(norm(norm(Fk)))]);
            end
        end
        if norm(Fk)<par.Results.tol
            break
        end
        i=i+1;
    end
end
%%% calculate the final function value
[ek,gk] = feval(FUN,xk);
%%% return the results
out.x=xk;
out.x_iter=x_iter;
out.f=ek;
out.g=gk;
out.iter=i;
end

