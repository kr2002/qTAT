%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the objective function and its gradients with 
% respect to the optimization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f g]=qTATObj(X,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,...
    Ns,Hm,SrcInfo,BdaryInfo,wnum)

M=Nx*Ny; % total number of nodes in the mesh
ne = size(SrcInfo,2); % number of edges/nodes on the domain boundary

refc=X(1:M);% current value of n
sigmac=X(M+1:2*M); % current value of sigma

f=0.0;
g=zeros(2*M,1);
for ks=1:Ns
    
    Hc=zeros(M,1); % predicted data
    rz=zeros(M,1); % residual on measurement locations
    srczero=zeros(M,1); % zero volume source for forward problems
 
    uc=HelmholtzSolve('Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,srczero);
    Hc=sigmac.*abs(uc).^2;
    
    %Hcg=tri2grid(P,T,Hc,x,y);
    %figure;
    %pcolor(x,y,Hcg); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    
    HmL=Hm(:,ks);
    rz=(Hc-HmL); % for unnormalized objective function
    %rz=(Hc-HmL)./HmL; % for normalized objective function
    
    % the contribution to the objective function from source ks
    f=f+0.5*sum(rz.^2)*dx*dy;
    
    % the contribution to the gradient from source ks
    if nargout > 1         

        % solve the adjoint equation
        srcadj=-sigmac.*rz.*conj(uc);        
        wc=HelmholtzSolve('Adjoint',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,srcadj);
        
        %wcg=tri2grid(P,T,wc,x,y);
        %figure;
        %pcolor(x,y,real(wcg)); axis tight; colorbar('SouthOutside');
        %axis square; axis off; shading interp;
        %drawnow;
        %pause;
    
        % the gradient w.r.t n            
        if strcmp(MinVar,'Ref')||strcmp(MinVar,'All')
            g(1:M)=g(1:M)+2*wnum^2*real(uc.*wc)*dx*dy;
        end
        % the gradient w.r.t sigma
        if strcmp(MinVar,'Sigma')||strcmp(MinVar,'All')
            g(M+1:2*M)=g(M+1:2*M)+(rz.*abs(uc).^2+2*wnum*real(i*uc.*wc))*dx*dy;
        end
        
    end
    
end

% Add regularization terms to both the objective function and its gradients
betan=0e-16; betaS=1*betan; % regularization parameters

if strcmp(MinVar,'Ref')||strcmp(MinVar,'All')
    [Rx,Ry] = pdegrad(P,T,refc);
    Rx1=pdeprtni(P,T,Rx); Ry1=pdeprtni(P,T,Ry);
    f=f+0.5*betan*sum(Rx1.^2+Ry1.^2)*dx*dy;
    if nargout >1
        [Rxx, Rxy]=pdegrad(P,T,Rx1); [Ryx, Gyy]=pdegrad(P,T,Ry1);
        Rx2=pdeprtni(P,T,Rxx); Ry2=pdeprtni(P,T,Gyy);
        Deltan=Rx2+Ry2;
        g(1:M)=g(1:M)-betan*Deltan*dx*dy;
        for j=1:ne
            nd=BdaryInfo(1,j);
            g(nd)=g(nd)-betan*BdaryInfo(3,j)*Rx1(nd)+BdaryInfo(4,j)*Ry1(nd)*BdaryInfo(5,j);
        end
    end
end
if strcmp(MinVar,'Sigma')||strcmp(MinVar,'All')
    [Sx,Sy] = pdegrad(P,T,sigmac);
    Sx1=pdeprtni(P,T,Sx); Sy1=pdeprtni(P,T,Sy);
    f=f+0.5*betaS*sum(Sx1.^2+Sy1.^2)*dx*dy;
    if nargout >1
        [Sxx, Sxy]=pdegrad(P,T,Sx1); [Syx, Syy]=pdegrad(P,T,Sy1);
        Sx2=pdeprtni(P,T,Sxx); Sy2=pdeprtni(P,T,Syy);
        DeltaSigma=Sx2+Sy2;
        g(M+1:2*M)=g(M+1:2*M)-betaS*DeltaSigma*dx*dy;
        for j=1:ne
            nd=BdaryInfo(1,j);
            g(M+nd)=g(M+nd)-betaS*BdaryInfo(3,j)*Sx1(nd)+BdaryInfo(4,j)*Sy1(nd)*BdaryInfo(5,j);
        end
    end
end