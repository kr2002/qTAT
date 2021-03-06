%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qTAT: qTAT image reconstruction with least-square minimization method 
%       based on the Helmholtz model with Dirichlet boundary conditions
%
% Author:      Kui Ren  ren@math.utexas.edu
% Institution: Department of Mathematics, UT Austin
% Last update: 2010-10-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathematical model:
% \Delta u + k^2 (1+n) u +i k \sigma u =0
% u = f
% 
% Internal data: H = \sigma |u|^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  close all;

tic; tb=toc;

% Load geometrical information on the domain
load 'geo-2b2';

MaxIT=50;

Ns=36;

% Set wave number
wnum=1;

% Decide which parameter to be reconstructed
MinVar='Sigma'; % 'Ref','Sigma' or 'All'

% Create a Cartesian grid for inversion
dx=0.025; x=0:dx:2;
dy=0.025; y=0:dy:2;
Nx=length(x);
Ny=length(y);
% [X,Y]=meshgrid(x,y);

% Generate regular finite element mesh on rectangular geometry
[P,E,T]=poimesh(geo,Nx-1,Ny-1); 
%figure;
%pdemesh(P,E,T);
%axis tight; axis square; box on; axis off;
%title('Finite element mesh');

M=Nx*Ny; % total number of nodes in the mesh

% Setup information on sources and boundary edges
SrcInfo=SetSources(Ns);
BdaryInfo=SetBdaryInfo(P,E);

% Set true parameters for wave propagation
rec1=[0.5 0.8; 0.5 0.8]; 
rec2=[1.4 1.7; 1.4 1.7];
rec3=[1.0 1.8; 0.3 0.6];
circ1=[1.0 1.5 0.2];
circ2=[1.5 1.0 0.3];

%reft=0.3*ones(M,1); % true refractive index
%reft=(0.1+0.2*exp(-(P(1,:)-1).^2-(P(2,:)-1).^2))';
reft=0.2+0.2*ind_rec(P,rec1)+0.4*ind_rec(P,rec2)+0.6*ind_rec(P,rec3);

sigmat=0.02*ones(M,1); % true absorption
%sigmat=(0.1+0.2*exp(-(P(1,:)-1).^2-(P(2,:)-1).^2))';
%sigmat=0.2+0.2*ind_rec(P,rec1)+0.4*ind_rec(P,rec2)+0.6*ind_rec(P,rec3);
for k1=1:10
    for k2=1:10
        dom=[0.5+(k1-1)*0.1 0.5+k1*0.1;0.5+(k2-1)*0.1 0.5+k2*0.1];
        sigmat=sigmat+0.02*(sign(rand-0.5)+1)*ind_rec(P,dom);
    end
end

if strcmp(MinVar,'Ref')||strcmp(MinVar,'All')
    reftg=tri2grid(P,T,reft,x,y);
    figure;
    pcolor(x,y,reftg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    title('true n');
    drawnow;
end
if strcmp(MinVar,'Sigma')||strcmp(MinVar,'All')
    sigmatg=tri2grid(P,T,sigmat,x,y);
    figure;
    pcolor(x,y,sigmatg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    title('true \sigma');
    drawnow;
end

disp('Finished setting simulation geometry and parameters .......');

% Generating synthetic data
disp(' ');
disp(' ');
disp('Generating synthetic data .......');
disp(' ');

noiselevel=0.0; % set noise level
srczero=zeros(M,1);
Hm=zeros(M,Ns);
for ks=1:Ns
    
    % Solve the Helmholtz equation
    ut=HelmholtzSolve('Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,reft,sigmat,srczero);

    Ht=sigmat.*abs(ut).^2;
    
    % Plot data/solutions
    %utg=tri2grid(P,T,Ht,x,y);
    %figure;
    %pcolor(x,y,real(utg)); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    %pause
   
    % Add noise to data
	Hm(:,ks)=Ht.*(1+noiselevel*2*(rand(M,1)-0.5));
    
    disp(['Synthetic data generated for source #: ' num2str(ks)]);
    disp('  ');

    clear ut Ht;
end
disp('Finished generating synthetic data .......');

% Setup initial guess
disp(' ');
disp(' ');
disp('Setting initial guess .......');
disp(' ');

ref0=0.2*ones(M,1);
if ~strcmp(MinVar,'Ref')
    ref0=reft;
end

sigma0=0.02*ones(M,1);
if ~strcmp(MinVar,'Sigma')
    sigma0=sigmat;
end

if strcmp(MinVar,'Ref')||strcmp(MinVar,'All')
    ref0g=tri2grid(P,T,ref0,x,y);
    figure;
    pcolor(x,y,ref0g); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    %caxis([0.05 0.25]);
    title('initial guess of n');
    drawnow;
end
if strcmp(MinVar,'Sigma')||strcmp(MinVar,'All')
    sigma0g=tri2grid(P,T,sigma0,x,y);
    figure;
    pcolor(x,y,sigma0g); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    title('initial guess of \sigma');
    drawnow;
end

X0=[ref0' sigma0']';

disp('Finished setting initial guess .......');

% This short part is only for debugging
%[f0 g0]=qTATObj(X0,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,Ns,Hm,SrcInfo,BdaryInfo,wnum);
%g0g=tri2grid(P,T,g0(1:M),x,y);
%g0g=tri2grid(P,T,g0(M+1:2*M),x,y);
%figure;
%pcolor(x,y,g0g); axis tight; colorbar('SouthOutside');
%axis square; shading interp;
%title('Gradient');
%drawnow;

OptimMethod='UNCON';

% Setup the minimization algorithm
disp(' ');
disp(' ');
disp('Minimizing objective function .......');
disp(' ');

f=@(X) qTATObj(X,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,Ns,Hm,SrcInfo,BdaryInfo,wnum);

if strcmp(OptimMethod,'UNCON')
    options=optimoptions(@fminunc,'Algorithm','quasi-newton', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    [X,fval,exitflag,output,grad]=fminunc(f,X0,options);
else
    % Set inequality constraint
    Aieq=zeros(1,2*M);
    Bieq=0;
    % Set equality constraint
    Aeq=zeros(1,2*M);
    Beq=0;
    % Set upper and lower bounds
    LB=[0.1*ones(1,M) 0.1*ones(1,M)]';
    UB=[0.9*ones(1,M) 0.4*ones(1,M)]';

    options=optimoptions(@fmincon,'Algorithm','trust-region', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    %options=optimset('Display','iter-detailed','GradObj','on','TolFun',1e-12,'MaxIter',MaxIT);
    %options = optimset('algorithm','sqp','maxfunevals',5000,'maxiter',100);
    %options = optimset(options,'tolx',1e-9,'tolcon',1e-9,'tolfun',1e-6);
    %options = optimset(options,'GradObj','on','GradConstr','off');
    
    [X,fval,exitflag,output,lambda]=fmincon(f,X0,Aieq,Bieq,Aeq,Beq,LB,UB,[],options);
    %[X,fval,exitflag,output]=fmincon(f,X0,zeros(M,M),zeros(M,1),[],[],LB,UB);
end

disp(' ');
disp(' ');
disp('Finished minimizing objective function .......');

disp(' ');
disp(' ');
disp('Plotting final results .......');
disp(' ');

refr=X(1:M);
sigmar=X(M+1:2*M);
% Plot reconstruction results
if strcmp(MinVar,'Ref')||strcmp(MinVar,'All')
    refrg=tri2grid(P,T,refr,x,y);
    figure;
    pcolor(x,y,refrg); axis tight; colorbar('SouthOutside');
    %caxis([0.40 1.10]);
    axis square; axis off; shading interp;
    title('recovered n');
    drawnow;
end
if strcmp(MinVar,'Sigma')||strcmp(MinVar,'All')
    sigmarg=tri2grid(P,T,sigmar,x,y);
    figure;
    pcolor(x,y,sigmarg); axis tight; colorbar('SouthOutside');
    %caxis([0.10 0.50]);
    axis square; axis off; shading interp;
    title('recovered \sigma');
    drawnow;
end

disp('Finished plotting final results .......');

save Exp01-Info geo P E T SrcInfo BdaryInfo wnum Ns MaxIT ...
                  OptimMethod noiselevel dx dy -ASCII
save Exp01-Results reft ref0 refr sigmat sigma0 sigmar -ASCII

te=toc;
disp(' ');
disp(' ');
disp(['The code run for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

% This last line is used to close MATLAB after the computation. It is 
% only used when runing the code in background.

%exit; % to exit MATLAB