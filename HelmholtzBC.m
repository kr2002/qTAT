%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the forward and the adjoint 
% Helmholtz problems.
%
% The forward Helmholtz model:
%
% \Delta u + k^2(1+n)u+ ik\sigma u = S  in \Omega
% u = f, on \partial\Omega
%
% with S=0, f given by the boundary source
%
% The adjoint Helmholtz model is:
%
% \Delta w + k^2(1+n)u + ik\sigma w = S  in \Omega
% w=f, on \partial\Omega
%
% with S=-(\sigma |u|^2-H)\sigma \conj(u), f=0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qmatrix,gmatrix,hmatrix,rmatrix] = ...
         HelmholtzBC(Type,SrcInfo,BdaryInfo,ks,p,e,u,time)

ne = size(e,2); % number of edges on the domain boundary

% This two lines indicate that the BC is NOT Robin
qmatrix = zeros(1,ne);
gmatrix = zeros(1,ne);

% The following two lines set the BC to homogeneous Dirichlet type
hmatrix = ones(1,2*ne);
rmatrix = zeros(1,2*ne);

% Set the rmatrix
if strcmp(Type,'Adjoint') % boundary source for adjoint: zero
    for k = 1:ne
        rmatrix(k)=0.0;
        rmatrix(ne+k)=0.0;
    end
elseif strcmp(Type,'Forward') % boundary ources for forward: Gaussians
    xs=SrcInfo(1,ks);
    ys=SrcInfo(2,ks);
    srcseg=SrcInfo(3,ks);
    for k = 1:ne
        x1 = p(1,e(1,k)); % x at first point in segment
        y1 = p(2,e(1,k)); % y at first point in segment
        x2 = p(1,e(2,k)); % x at second point in segment
        y2 = p(2,e(2,k)); % y at second point in segment
        rmatrix(k)=1.0;
        if BdaryInfo(2,k)==srcseg % if the edge lives on the same side with the source
            rmatrix(k) = 2*exp(-((x1-xs)^2+(y1-ys)^2)/0.1);
            rmatrix(k+ne) = 2*exp(-((x2-xs)^2+(y2-ys)^2)/0.1);
        end        
    end
else
    disp('Must specific problem type (Forward or Adjoint) to fix BC!');
end