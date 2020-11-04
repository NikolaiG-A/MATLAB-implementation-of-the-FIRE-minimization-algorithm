function [x,v,F] = MDintegrator(FUN,xold,vold,Fold,varargin)
%%% integration sted, i.e. x=xold+dx, v=vold+dv
%%% Fold - force vector at the current step

% Create parser
par = inputParser;

% Set parameters
par = params(par);

% Parse input
par.parse(varargin{:});

%%% initialise parameters
dt=par.Results.dt; %%% time step
alpha=par.Results.alpha; %%% mixing parameter (FIRE algorithm)
m=par.Results.mass; %%% mass of the particles



if strcmpi(par.Results.Integrator,'euler')
    %%% Mixing (FIRE algorithm)
    vold=(1-alpha)*vold+alpha*norm(vold)/norm(Fold);
    %%% new locations
    x=xold+dt*vold;
    %%% new forces
    [~,F] = feval(FUN,x);
    F=-F;
    %%% new speeds
    v=vold+dt*F/m;
elseif strcmpi(par.Results.Integrator,'si_euler')
    %%% modify the speed 
    v=vold+dt*Fold/m;
    %%% Mixing (FIRE algorithm)
    v=(1-alpha)*v+alpha*norm(v)/norm(Fold);
    %%% new locations
    x=xold+dt*vold;
    %%% new forces
    [~,F] = feval(FUN,x);
    F=-F;
elseif strcmpi(par.Results.Integrator,'leapfrog')
    %%% initialise previous half-step
    v_m=vold-0.5*dt*Fold/m;
    %%% initialise the forward half-step
    v_p=v_m+dt*Fold/m;
    %%% Mixing (FIRE algorithm)
    v_p=(1-alpha)*v_p+alpha*Fold*norm(v_p)/norm(Fold);
    %%% new locations
    x=xold+dt*v_p;
    %%% new forces
    [~,F] = feval(FUN,x);
    F=-F;
    %%% new speeds
    v=vold+dt*F/m;
elseif strcmpi(par.Results.Integrator,'verlet')
    %%% initialise the forward half-step
    v_p=vold+0.5*dt*Fold/m;
    %%% Mixing (FIRE algorithm)
    v_p=(1-alpha)*v_p+alpha*Fold*norm(v_p)/norm(Fold);
    %%% new locations
    x=xold+dt*v_p;
    %%% new forces
    [~,F] = feval(FUN,x);
    F=-F;
    %%% new speeds
    v=v_p+0.5*dt*F/m;
end
end

