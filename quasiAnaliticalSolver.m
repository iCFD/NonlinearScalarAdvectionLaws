function [X,U,PHI]=quasiAnaliticalSolver(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Quasi-analytical solutions for scalar nonlinear conservation models
%           Here these models are assumed to be of the form:
%
%                      u_t - f(u)_x = 0,    {x>R, t>0}
%
%                Coded by Manuel A. Diaz, ISNA, 2018.07.20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref.: 
% [1] Coulouvrat, François. "A quasi-analytical shock solution for general
%     nonlinear progressive waves." Wave Motion 46.2 (2009): 97-107. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: We follow the sign convention of the flux function as in Ref.[1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; close all; clc;

switch nargin
    case 0 % Manually set base parameters
        fluxFunc = 'buckley';
        t = 2.4; % output time
        debug = true;
    case 4
        x=varargin{1};
        u0=varargin{2};
        t=varargin{3};
        fluxFunc=varargin{4};
        debug = false;
    case 5
        x=varargin{1};
        u0=varargin{2};
        t=varargin{3};
        fluxFunc=varargin{4};
        debug=varargin{5};
    otherwise
        error('Incorrect input arguments')
end

%% 0. Flux functions and their exact derivatives
switch fluxFunc
    case 'burgers' % inviscid burgers' equation (traditional sign)
        f = @(u) -u.^2/2; % the flux function
        df= @(u) -u; % the derivative f'(u)
        Df= @(u) -u.^2/2; % and Df: u*f'(u)-f(u)
        I = @(u) -u.^3/3; % I(u): energy flux
	case 'burgers+' % inviscid burgers' equation (Couluovrat convention)
        f = @(u) u.^2/2; % the flux function
        df= @(u) u; % the derivative f'(u)
        Df= @(u) u.^2/2; % and Df: u*f'(u)-f(u)
        I = @(u) u.^3/3; % I(u): energy flux
    case 'cubic' % inviscid burgers' equation (traditional sign)
        f = @(u) -u.^3/3; % the flux function
        df= @(u) -u.^2; % the derivative f'(u)
        Df= @(u) -2*u.^3/3; % and Df: u*f'(u)-f(u)
        I = @(u) -u.^4/4; % I(u): energy flux
	case 'cubic+' % inviscid burgers' equation (Couluovrat convention)
        f = @(u) u.^3/3; % the flux function
        df= @(u) u.^2; % the derivative f'(u)
        Df= @(u) 2*u.^3/3; % and Df: u*f'(u)-f(u)
        I = @(u) u.^4/4; % I(u): energy flux
    case 'buckley' % buckley-leverett equation (traditional sign, but not working yet!)
        a = 0.5; % constant parameter 
        f = @(u) -(u.^2)./(u.^2 + a*(1-u).^2); % the flux function
        df= @(u) -2*a*u.*(1-u)./(u.^2 + a*(1-u).^2).^2; % the derivative f'(u)
        Df= @(u) (a*(u.^2-4*u+3).*u.^2+u.^4)./(u.^2 + a*(1-u).^2).^2; % and Df: u*f'(u)-f(u)
        w = @(u) (a+1)*u/sqrt(a)-sqrt(a); % dummy function
        I = @(u) (sqrt(a)/(1+a)^2)*(w(u)+(a-1)*atan(w(u))+sqrt(a)*log(1+w(u).^2)); % I(u): energy flux
    case 'buckley+' % buckley-leverett equation (Couluovrat convention)
        a = 0.5; % constant parameter 
        f = @(u) (u.^2)./(u.^2 + a*(1-u).^2); % the flux function
        df= @(u) 2*a*u.*(1-u)./(u.^2 + a*(1-u).^2).^2; % the derivative f'(u)
        Df= @(u) (u.^2).*(a-(1+a)*u.^2)./(u.^2 + a*(1-u).^2).^2; % and Df: u*f'(u)-f(u)
        w = @(u) (a+1)*u/sqrt(a)-sqrt(a); % dummy function
        I = @(u) (sqrt(a)/(1+a)^2)*(w(u)+(a-1)*atan(w(u))+sqrt(a)*log(1+w(u).^2)); % I(u): energy flux
    otherwise
        error('case not set');
end

%% 1. Set initial conditions for the scalar (u) and potential (Ø) equation
if nargin==0
    dx=pi/100; x=(-5*pi:dx:5*pi); nx=numel(x);
    u0 = zeros(1,nx) + sin(x).*exp(-x.^2/50).*(abs(x) < 4*pi);
else
    nx=numel(x); dx=max(x(2:nx)-x(1:nx-1));
end
phi0=zeros(1,nx); for j=1:nx, phi0(j)=phi0(j)+sum(dx*u0(1:j)); end;
if debug, figure(1); end
if debug, subplot(321); plot(x,u0,'.k'); ylabel('u(x,t)'); end
if debug, subplot(322); plot(x,phi0,'.k'); ylabel('\phi(x,t)'); end

%% 2. The Poisson solution
% The Poisson problem tell us that the solution is left similar and it can
% obtained by simply shifting the abscisas of the original profile. 
% Therefore, the new y-coordinates and the phi-operator are computed as
y = x-t*df(u0); phi=phi0-t*Df(u0);
if debug; subplot(323); plot(y,u0,':k'); hold on; end
if debug; subplot(324); plot(y,phi,':k'); hold on; end

%% 3. Identify positive branches 
% 3.1 Check the monotonicity of y-absissas
TVD=y(2:nx)-y(1:nx-1)<0; y(TVD)=NaN; % If i-cell is TVD, set y(i)=NaN

% 3.2 Find indexes of the positive branches
idx=find(diff(isnan(y(1:nx-1)))); nShocks=numel(idx)/2; nBranches=nShocks+1;
idxBranches=[0,idx,nx]; Branches=reshape(idxBranches,[2,nBranches]);
Branches(1,:)=Branches(1,:)+1; % Correct for the initial index of each branch
if isnan(prod(y(Branches(:))))
    nBranches=numel(idx)/2; nShocks=nBranches-1;
    idxBranches=idx; Branches=reshape(idxBranches,[2,nBranches]);
    Branches(1,:)=Branches(1,:)+1;
end
if debug; subplot(323); plot(y(Branches),u0(Branches),'o'); hold on; end
if debug; subplot(324); plot(y(Branches),phi(Branches),'o'); hold on; end

%% 4. Interpolation of positive branches 
Branch{nBranches}.x=[];
Branch{nBranches}.u=[];
Branch{nBranches}.phi=[];
Branch{nBranches}.x_sh=[];
Branch{nBranches}.EC=[];
for j = 1:nBranches
    Branch{j}.x = x(x>=y(Branches(2*j-1)) & x<=y(Branches(2*j)));
    Branch{j}.u = interp1(y(Branches(2*j-1):Branches(2*j)),u0(Branches(2*j-1):Branches(2*j)),Branch{j}.x,'linear');
    Branch{j}.phi=interp1(y(Branches(2*j-1):Branches(2*j)),phi(Branches(2*j-1):Branches(2*j)),Branch{j}.x,'linear');
    if debug; subplot(323); plot(Branch{j}.x,Branch{j}.u); hold on; end
    if debug; subplot(324); plot(Branch{j}.x,Branch{j}.phi); hold on; end
end

%% 5. Identify intersections of positive branches
for j=1:nBranches-1
    intersection = NaN*ones(1,nBranches);
    up_sh = NaN*ones(1,nBranches);
    un_sh = NaN*ones(1,nBranches);
    entropyCondition = zeros(1,nBranches);
    
    for k=j+1:nBranches
        [x,jB,kB]=intersect(Branch{j}.x,Branch{k}.x); % Evaluate if branches x-range intersect
        % If phi-branches intersect in the solution grid
        if numel(x)>0 
            % Ø-Ø1 = (Ø2 -Ø1 )/(x2-x1).(x-x1) : line Eq.(1),
            % Ø-Ø1*= (Ø2*-Ø1*)/(x2-x1).(x-x1) : line Eq.(2),
            x1=x(1:end-1); phi1=Branch{j}.phi(jB(1:end-1)); phi1s=Branch{k}.phi(kB(1:end-1));
            x2=x( 2:end ); phi2=Branch{j}.phi(jB( 2:end )); phi2s=Branch{k}.phi(kB( 2:end ));
            % The intersection is computed as:
            % Ø1*-Ø1= ((Ø2-Ø1)-(Ø2*-Ø1*))/(x2-x1).(x-x1) : Eq.(1) - Eq.(2).
            x = x1 + (phi1s-phi1).*(x2-x1)./((phi2-phi1)-(phi2s-phi1s)); 
            if numel(x(x>x1 & x<x2))>0 
                intersection(k) = x(x>x1 & x<x2); % then
                up_sh(k)= interp1(y(Branches(2*k-1):Branches(2*k)),u0(Branches(2*k-1):Branches(2*k)),intersection(k),'linear');
                un_sh(k)= interp1(y(Branches(2*j-1):Branches(2*j)),u0(Branches(2*j-1):Branches(2*j)),intersection(k),'linear');
                phi_sh = interp1(y(Branches(2*k-1):Branches(2*k)),phi(Branches(2*k-1):Branches(2*k)),intersection(k),'linear');
                entropyCondition(k) = (I(up_sh(k))-I(un_sh(k))) > ((f(up_sh(k))-f(un_sh(k))).*(up_sh(k)+un_sh(k))/2);
            if debug; subplot(323); plot([intersection(k),intersection(k)],[up_sh(k),un_sh(k)],'.-k'); hold on; end
            if debug; subplot(324); plot(intersection(k),phi_sh,'.k'); hold on; end
            end
        end
        % If phi-branches do not explicitly in the solution grid
        if j==k-1 && Branch{j}.x(end)-Branch{k}.x(1)<dx
            intersection(k) = (Branch{j}.x(end)+Branch{k}.x(1))/2; % then
            up_sh(k)= interp1(y(Branches(2*k-1):Branches(2*k)),u0(Branches(2*k-1):Branches(2*k)),intersection(k),'linear');
            un_sh(k)= interp1(y(Branches(2*j-1):Branches(2*j)),u0(Branches(2*j-1):Branches(2*j)),intersection(k),'linear');
            phi_sh = interp1(y(Branches(2*k-1):Branches(2*k)),phi(Branches(2*k-1):Branches(2*k)),intersection(k),'linear');
            entropyCondition(k) = (I(up_sh(k))-I(un_sh(k))) > ((f(up_sh(k))-f(un_sh(k))).*(up_sh(k)+un_sh(k))/2);
            if debug; subplot(323); plot([intersection(k),intersection(k)],[up_sh(k),un_sh(k)],'.-k'); hold on; end
            if debug; subplot(324); plot(intersection(k),phi_sh,'.k'); hold on; end
        end
    end
    Branch{j}.x_sh = intersection;
    Branch{j}.up_sh= up_sh;
    Branch{j}.un_sh= un_sh;
    Branch{j}.EC=entropyCondition;
end
if debug; subplot(323); hold off; ylabel('u(x,t)'); end
if debug; subplot(324); hold off; ylabel('\phi(x,t)'); end

%% 6. Selection of admisible shocks
% Among all possible intersections jumping from a branch j to a branch k>j,
% rule out:
% (a) those which do not satisfy the entropy condition Eq.(20) in Ref.[1];
% (b) those which do not cross a branch of higher number;
% (c) similarly, those which do come from a branch that does not cross a 
%     previously lower branch.
% 

% Examing the full path (The simplest solution first!)
path = nchoosek(1:nBranches,nBranches); 
shocks = zeros(1,nShocks);
for p = 1:size(path,1)
    for b = 1:nShocks
        from=path(p,b); to=path(p,b+1);
        shocks(p,b) = Branch{from}.x_sh(to);
    end
end
if issorted(shocks) % if shocks are in order: the full path is the solution
    solutionPath=path;
    intersections=shocks;
else % if the full solution path is not the solution then,
    i=0; idx=[];
    while numel(idx)==0
        i = i+1; 
        path = nchoosek(1:nBranches,nBranches-i); % List all possible paths
        shocks= NaN*ones(size(path,1),nShocks-i);
        index = ones(size(path,1),1);
        for p = 1:size(path,1)
            if path(p,1)==1 && path(p,nBranches-i)==nBranches
                for b = 1:nShocks-i
                    from=path(p,b); to=path(p,b+1);
                    shocks(p,b) = Branch{from}.x_sh(to);
                end
            end
        end 
        for b=1:nShocks-1-i
            index = index.*(shocks(:,b) < shocks(:,b+1));
        end
        idx=find(index);
    end
    solutionPath=path(idx,:);
    intersections=shocks(idx,:);
end

%% 7. Build (and plot) the Solution Profile
X=[]; U=[]; PHI=[];
if numel(Branch)==1
    X=[X, Branch{1}.x];
    U=[U, Branch{1}.u];
    PHI=[PHI, Branch{1}.phi];
else
    for i = 1:length(solutionPath)
        if i ==1 
            X=[X, Branch{1}.x(Branch{1}.x<=intersections(1))]; %#ok<*AGROW>
            U=[U, Branch{1}.u(Branch{1}.x<=intersections(1))];
            PHI=[PHI, Branch{1}.phi(Branch{1}.x<=intersections(1))];
        elseif i==length(solutionPath)
            X=[X, Branch{nBranches}.x(Branch{ nBranches }.x>=intersections(i-1))];
            U=[U, Branch{nBranches}.u(Branch{ nBranches }.x>=intersections(i-1))];
            PHI=[PHI, Branch{nBranches}.phi(Branch{ nBranches }.x>=intersections(i-1))];
        else
            X=[X, Branch{solutionPath(i)}.x(Branch{solutionPath(i)}.x>intersections(i-1) & Branch{solutionPath(i)}.x<intersections(i))];
            U=[U, Branch{solutionPath(i)}.u(Branch{solutionPath(i)}.x>intersections(i-1) & Branch{solutionPath(i)}.x<intersections(i))];
            PHI=[PHI, Branch{solutionPath(i)}.phi(Branch{solutionPath(i)}.x>intersections(i-1) & Branch{solutionPath(i)}.x<intersections(i))];
        end
    end
end
if debug; subplot(325); plot(X,U,'-k'); ylabel('u(x,t)'); xlabel('x'); end
if debug; subplot(326); plot(X,PHI,'-k'); ylabel('\phi(x,t)'); xlabel('x'); end

% Manuel A. Diaz, ISNA 2018, Santa Fe, NM.
