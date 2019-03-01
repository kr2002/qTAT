%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves the the forward and adjoint Helmholtz 
% equations with a P_1 finite element method
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=HelmholtzSolve(Type,SrcInfo,BdaryInfo,ks,P,E,T,wnum,ref,sigma,S)

% interpolation to triangle middle point
refm=pdeintrp(P,T,ref);
sigmam=pdeintrp(P,T,sigma);
Sm=pdeintrp(P,T,S);

% construct mass matrices
[K,M,F]=assema(P,T,-1,wnum^2*(1+refm)+i*wnum*sigmam,Sm);

% construct boundary conditions
pdebound =@(p,e,u,time)HelmholtzBC(Type,SrcInfo,BdaryInfo,ks,p,e,[],[]);
[Q,G,H,R] = assemb(pdebound,P,E);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);