function RunCoulouvratSolver(xmin,xmax,IC,fluxFunc,tEnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Quasi-analytical solutions for scalar nonlinear conservation models
%                    assumed to be of the form:
%
%                      u_t + f(u)_x = 0,    {x>R, t>0}
%
%                Coded by Manuel A. Diaz, ISNA, 2018.07.20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref.: 
% [1] Coulouvrat, François. "A quasi-analytical shock solution for general
%     nonlinear progressive waves." Wave Motion 46.2 (2009): 97-107. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: The present implementation follows the sign convention of the flux
%       function as in Ref.[1]. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; close all; clc;

% Basic Parameters
    %tEnd = 20; % output time
      %IC = 11;  % 
%fluxFunc = {'buckley','burgers','cubic'};
   debug = true;

% Build base mesh
dx=1/100; x=(xmin:dx:xmax); 

% Load and initial condition
u0=CommonIC(x,IC); % <-- See more details in CommonIC.m

% Compute exact solution
[xe,ue]=quasiAnaliticalSolver(x,u0,tEnd,fluxFunc,debug);

% Plot IC and exact Solution
figure(2); plot(x,u0,'-.k',xe,ue,'-r'); 
title(sprintf('t=%g [-]',tEnd)); ylabel('u(x,t)'); xlabel('x');
legend('Initial Condition','Quasi-analytical'); legend boxoff; 

% Save Exact Solution 
save('quasiAnalitical.mat','xe','ue','tEnd');