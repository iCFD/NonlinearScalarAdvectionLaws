%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 1-D wave equation Hybrid numerical schemes, namely
%
%             Tangent Hyperbola for INterface Capturing and
%             Monotonic Upwind Scheme for Conservation Laws 
%                          (MUSCL-THINC-BVD) scheme
%
%                                  &
%
%               Tangent Hyperbola for INterface Capturing and
%                    Weighted Essentially Non-Oscilaroty
%                          (WENO5-THINC-BVD) scheme
%
%
%                 du/dt + df/dx = S, for x \in [a,b]
%                  where f = f(u): nonlinear flux only!
%                     and S = s(u): source term
%
%
%             coded by Manuel Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NHRI, 2018.06.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: 
% [1] Deng, Xi, Bin Xie, and Feng Xiao. "Some practical versions of
%     boundary variation diminishing (BVD) algorithm." arXiv preprint
%     arXiv:1708.01148 (2017).  
% [2] Deng, Xi, et al. "Limiter-free discontinuity-capturing scheme for 
%     compressible gas dynamics with reactive fronts." C & F (2018).
% [3] Deng, Xi, et al. "High fidelity discontinuity-resolving
%     reconstruction for compressible multiphase flows with moving 
%     interfaces." Journal of Computational Physics (2018). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes: 
% -------
% The present implementation serves for the purpose of comparions and a
% summary of the main hybrid algorithms reported in [1-3]. In this
% snipets, I prioritize readability rather than code speed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %close all; clc;

%% Parameters
   nx = 200;	% number of cells
  CFL = 0.70;	% Courant Number
 tEnd = 100.;   % End time
limit = 'MM';   % MC, MM, VA. (only for MUSCL Methods)
scheme = 5;     % 1:WENO5, 
                % 2:WENO5-THINC-BVD,  
                % 3:MUSCL, 
                % 4:MUSCL-THINC-BVD,
                % 5:THINC-BVD,

fluxfun='buckley+'; % select flux function
% Define our Flux function
switch fluxfun 
    case 'burgers' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^2/2; 
        dflux = @(w) w;
        ICcase=2; CFL=0.4; IC=11; tEnd=20;
    case 'cubic' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^3/3; 
        dflux = @(w) w.^2;
        ICcase=2; CFL=0.4; IC=11; tEnd=100;
    case 'buckley+' % Buckley-Leverett, CFL_max: 0.20 & tEnd: 0.40
        c = 0.5; % constant parameter
        flux = @(w) -(w.^2)./(w.^2 + c*(1-w).^2);
        dflux = @(w) -2*c*w.*(1-w)./(w.^2 + c*(1-w).^2).^2;
        ICcase=2; CFL=0.2; IC=11; tEnd=3.6;
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add', S = @(w) 0.1*w.^2;
    case 'dont',S = @(w) zeros(size(w));
end

% Build discrete domain
a=-4*pi; b=4*pi; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
%ICcase=2;  % overide for: {1}Testing, {2}Costum ICs
switch ICcase
    case 1, u0=TestingIC(x);  tEnd=2.0; % Jiang and Shu Testing IC
    case 2, u0=CommonIC(x,IC); % scalar cases 1-10 <- check them out!
    otherwise, error('IC file not listed');
end

% Build Exact Solution using the quasilinear solver
RunCoulouvratSolver(a,b,IC,fluxfun,tEnd);
load('quasiAnalitical.mat');%<-- Build solution with the quasilinear solver

switch scheme
    case 1, AdvecRes1d = @WENO5_AdvecRes1d_FluxSplitting; name='WENO';
    case 2, AdvecRes1d = @WENO5_THINC_AdvecRes1d; name='WENO5-THINC';
    case 3, AdvecRes1d = @MUSCL_AdvecRes1d_FV; name='MUSCL';
    case 4, AdvecRes1d = @MUSCL_THINC_AdvecRes1d; name='MUSCL-THINC';
    case 5, AdvecRes1d = @THINC_AdvecRes1d; name='THINC-BVD';
end

% Plot range
dl=0.1; plotrange=[a,b,min(u0)-dl,max(u0)+2*dl];

%% Solver Loop

% load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    % Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end

    % Update time and iteration counter
    t=t+dt; it=it+1;

    % RK Initial step
    uo = u;

    % 1st stage
    L = AdvecRes1d(u,flux,dflux,S,dx,limit);
    u = uo-dt*L;

    % 2nd Stage
    L = AdvecRes1d(u,flux,dflux,S,dx,limit);
    u = 0.75*uo+0.25*(u-dt*L);

    % 3rd stage
    L = AdvecRes1d(u,flux,dflux,S,dx,limit);
    u = (uo+2*(u-dt*L))/3;

    % Plot solution
    if rem(it,10) == 0
        plot(x,u0,'-',xe,ue,'-',x,u,'o'); axis(plotrange); shg; drawnow;
    end
end

%% Final Plot
plot(x,u0,'-',xe,ue,'-',x,u,'o'); axis(plotrange);
title('Cell averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);
legend({'IC','Quasi-Analytical',name},...
    'orientation','horizontal','location','north'); legend boxoff;