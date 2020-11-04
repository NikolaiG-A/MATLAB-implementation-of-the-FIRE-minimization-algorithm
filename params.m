function p = params(p)
% initialising parameters

% dt - initial time step (default dt=1e-8)
% alpha - initial value of parameter alpha from 0 to 1 (default alpha=0.1)
% f_alpha - change of alpha,i.e. alpha=f_alpha*alpha, around 1 (default f_alpha=0.99)
% dt_inc - adjusted time step incriment, i.e. dt_new=dt*dt_inc where dt_inc>1 (default dt_inc=2)
% dt_dec - decrease time step incriment, i.e. dt_new=dt*dt_dec where 0<dt_dec<1 (default dt_inc=0.5)
% dt_max - maximum time step (default dt_max=1e-1)
% dt_min - maximum time step (default dt_min=1e-10)
% IterMax - maximum number of iterations (default IterMax=1e6)
% tol - tolerance of the gradient norm, i.e. stopping criterion ||g||<tol (default IterMax=1e-7)
% Nmin - number of steps since the power P=v*F was negative to update dt and alpha
% NP_max - maximum number of iterations with P<=0
% initial_delay - 0 or 1, activates the delay in modifying dt and alpha
% iter_save - save every 'iter_save' solution during minimisation (default 0)
% iter_disp - display every 'iter_disp' result of norm grad during minimisation (default 0)
% mass - mass of the particles


% Integrator type:
%%% verlet - verlet algorithm (default)
%%% leapfrog - leapfrog scheme
%%% euler - explicit Euler algorithm
%%% si_euler - semi-implicit Euler algorithm

%%% assign parameters

addParameter(p,'dt',1e-8,@(x) x > 0);
addParameter(p,'alpha',0.1,@(x) x >=0);
addParameter(p,'f_alpha',0.99,@(x) x >=0);
addParameter(p,'dt_inc',2,@(x) x >= 1);
addParameter(p,'dt_dec',0.5,@(x) x <= 1 & x>0);
addParameter(p,'dt_max',1e-1,@(x) x >0);
addParameter(p,'dt_min',1e-10,@(x) x >0);
addParameter(p,'IterMax',1e6,@(x) x >= 1);
addParameter(p,'tol',1e-7,@(x) x > 0);
addParameter(p,'Nmin',5,@(x) x > 1);
addParameter(p,'NP_max',10,@(x) x > 1);
addParameter(p,'initial_delay',true,@islogical);
addParameter(p,'iter_save',0,@(x) x >=0);
addParameter(p,'iter_disp',0,@(x) x >=0);
addParameter(p,'mass',1,@(x) x > 0);
addParameter(p,'Integrator','verlet', @(x) ismember(x,{'verlet','leapfrog','euler','si_euler'}));

end