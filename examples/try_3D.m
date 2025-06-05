% @ -0,0 +1,47 @@
clc;
close all;
clear all;

%% add E I l A
syms E G L A Iz Iy Jk
Ke=[E*A/L 0 0 0 0 0 -E*A/L 0 0 0 0 0;
    0 12*E*Iz/L^3 0 0 0 6*E*Iz/L^2 0 -12*E*Iz/L^3 0 0 0 6*E*Iz/L^2;
    0 0 12*E*Iy/L^3 0 6*E*Iy/L^2 0 0 0 -12*E*Iy/L^3 0 6*E*Iy/L^2 0;
    0 0 0 G*Jk/L 0 0  0 0 0 -G*Jk/L 0 0 ;
    0 0 6*E*Iy/L^2 0 4*E*Iy/L 0 0 0 -6*E*Iy/L^2 0 2*E*Iy/L 0;
    0 6*E*Iz/L^2 0 0 0 4*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 2*E*Iz/L;
    -E*A/L 0 0 0 0 0 E*A/L 0 0 0 0 0;
    0 -12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2 0 12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2;
    0 0 -12*E*Iy/L^3 0 -6*E*Iy/L^2 0 0 0 12*E*Iy/L^3 0 -6*E*Iy/L^2 0;
    0 0 0 -G*Jk/L 0 0  0 0 0 G*Jk/L 0 0 ;
    0 0 6*E*Iy/L^2 0 2*E*Iy/L 0 0 0 -6*E*Iy/L^2 0 4*E*Iy/L 0;
    0 6*E*Iz/L^2 0 0 0 2*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 4*E*Iz/L]
[X,D] = eig(Ke) ;
r=rank(Ke)                      % rank of (A_1g)
d=diag(D)
X1=X(:,1:r)         % Mechanical mode 
X2=X(:,r+1:end)     % Pre stress mode 
 
mcs=rref(X1')       % Mechanical mode 
pre=rref(X2')       % Pre stress mode
%% beam
if 0
K=[1 0 0 -1 0 0
     0 12 6 0 -12 6;
     0 6 4 0 -6 2;
    -1 0 0 1 0 0
    0 -12 -6 0 12 -6;
    0 6 2 0 -6 4]
[U,S,V] = svd(K);
r=rank(K)                      % rank of (A_1g)
U1=U(:,1:r)
U2=U(:,r+1:end)        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r)
V2=V(:,r+1:end)       % V1 is C(A_1g'); V2 is N(A_1g) self stress mode

rref(U2')
rref(V1')
%%  Use eigenvalue analysis 
[X,D] = eig(K) ;
X1=X(:,1:r)         % Mechanical mode 
X2=X(:,r+1:end)     % Pre stress mode 
 
mcs=rref(X1')       % Mechanical mode 
pre=rref(X2')       % Pre stress mode 

end
%% plot  Pre stress mode 

N=[0 0 0;1 0 0]';
C_b=[-1 1];
C_s=[];
ne=1;
C_bar={[1 0;0 1]};
T_ei={eye(3)}
strut_s.T_ei=T_ei;
strut_s.C_bar=C_bar;


if isfield(strut_s,'displs')
strut_s=rmfield(strut_s,'displs');
end
for i=1:size(X2,2)          
strut_s.stress=kron(eye(ne),[kron(eye(2),[0 0 1])])*round(X2(:,i),3);     % Moment.
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['moment-',num2str(i)]);
axis([0 1 -0.8 0.8]);
strut_s.stress=kron(eye(ne),[kron(eye(2),[1 0 0])])*round(X2(:,i),3);     % axial force
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['axial force-',num2str(i)]);
axis([0 1 -0.8 0.8]);
strut_s.stress=kron(eye(ne),[kron(eye(2),[0 1 0])])*round(X2(:,i),3);     % shear force
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['shear force-',num2str(i)]);
axis([0 1 -0.8 0.8]);
end
%% Plot mechanism mode(use countor plot)
% Plot the structure to make sure it looks right
if isfield(strut_s,'stress')
strut_s=rmfield(strut_s,'stress');
end

for i=1:size(X1,2)  
    
fig=figure
tenseg_plot_dash(N,C_b,C_s,fig);
title('scissor hinge with cables');
N_d=[reshape(X1([1 2 4 5],i),2,[]);zeros(1,2)];
strut_s.displs=sqrt(sum(N_d.^2)');

N_motion=N+N_d;
tenseg_plot_stress(N_motion,C_b,C_s,fig,[],[],[],[],strut_s);

end
%% N C of the structure
% Manually specify node positions of double layer prism.
N=[0 0 0;1 1 0;2 0 0;1 -1 0;1.2 0.2 0]';  
n=N(:);
N2=N(1:2,:);                 %3D to 2D
n2=N2(:);
% Manually specify connectivity indices.
C_s_in = [1 2;2 3;3 4;4 1];  % This is indicating that string connection
C_b_in = [1 5;2 5;3 5;4 5];  % Similarly, this is saying bar 1 connects node 1 to node 2,

% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
C_abs=abs(C);
H=N*C';                     % element's direction matrix
l=sqrt(diag(H'*H));         % elements' length
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
n_m=sum(abs(C));        % n_m: No. of element in a node
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
title('scissor hinge with cables');
tenseg_plot(N,C,[]);
%% rigid/pin connection of members
cnct={0;0;0;0;{[1 3];[2 4]}};       %connection of rigid for 1; pin for 0; otherwise, a cell

% generate the E_nri:relation between rotation angle \theta, and reduced angle
% \theta_r,E_nr: collection of E_nri;E_eri: relation between element DOF
% and all DOF vector
% [E_nri,E_nr,E_eri]=tsgb_trans_dof(cnct,nn,n_m);

% generate the relation between rotation angle \theta, and reduced angle
% \theta_r
E_nri=cell(nn,1);
for i=1:nn
    if iscell(cnct{i})
         temp=zeros(n_m(i),numel(cnct{i}));
        for j=1:numel(cnct{i})
            I_temp=eye(n_m(i));
            temp(:,j)=sum(I_temp(:,cnct{i}{j}),2);            
        end
        E_nri{i}=temp;
    else if cnct{i}==0
        E_nri{i}=eye(n_m(i));
    else cnct{i}==1
            E_nri{i}=ones(n_m(i),1);    
    end
    end
end

E_nr=blkdiag(E_nri{:});         % r, for \theta rotation angle

n_cr=size(E_nr,2);              % num of \theta reduced

E_eri=cell(ne,1);               % \theta in a member
for i=1:ne
  E_eri{i}=zeros(2,n_cr) ; 

C_temp=C_abs;
C_temp(i,:)=C_temp(i,:)*2;
[~,~,v]=find(C_temp);
row_num=find(v==2);
E_eri{i}=E_nr(row_num,:);
end


%% C_bar
C_bar=cell(ne,1);
for i=1:ne
C_bar{i}=zeros(2,nn);
start_num=find(C(i,:)==-1);     %start node num
end_num=find(C(i,:)==1);       % end node num
C_bar{i}(1,start_num)=1;
C_bar{i}(2,end_num)=1;
end

E_ei=cell(nn,1);
for i=1:ne                      % relation between member DOF and minimal coordinate
E_ei{i}=blkdiag(kron(C_bar{i},eye(2)),E_eri{i});
end
E_e=cell2mat(E_ei);             % 
%% Boundary constraints
% node constraints
pinned_X=[]; pinned_Y=[]; 
[E_na,E_nb,a,b]=tenseg_boundary_2D(pinned_X,pinned_Y,nn);
% rotation constraints
ro_const={0;0;0;0;0};       %all constraint for -1; all free for 0; otherwise, a cell vector



% generate E_ra,E_rb,
E_rai=cell(nn,1);
E_rbi=cell(nn,1);

for i=1:nn

 if ro_const{i}==0
        E_rai{i}=eye(size(E_nri{i},2));
        E_rbi{i}=[];
    else if ro_const{i}==-1
             E_rai{i}=[];
            E_rbi{i}=ones(size(E_nri{i},2),1);   
    
 else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         I_temp=eye(size(E_nri{i},2));
         a_temp=ro_const{i};
 b_temp=setdiff((1:size(E_nri{i},2)),a_temp);
 E_rai{i}=I_temp(:,a_temp);
 E_rbi{i}=I_temp(:,b_temp);
    end
 end 
end

E_qa=blkdiag(E_na,blkdiag(E_rai{:}));           % relation between free DOF and minimal DOF
E_qb=blkdiag(E_nb,blkdiag(E_rbi{:}));           % relation between constraint DOF and minimal DOF
%% rotation matrix 
T_i=cell(ne,1);
T_ei=cell(ne,1);
dxi=cell(ne,1);
dyi=cell(ne,1);
for i=1:ne
dxi{i}=[1 0]*kron(C(i,:),eye(2))*n2;
dyi{i}=[0 1]*kron(C(i,:),eye(2))*n2;
T_i{i}=l(i)^-1*[dxi{i} -dyi{i};dyi{i} dxi{i}];
T_ei{i}=blkdiag(T_i{i},T_i{i},eye(2));
end
%% stiffness matrix
E=1e5*ones(ne,1);
r=0.1;
t=0.01;
I=pi/4*(r^4-(r-t)^4)*ones(ne,1);
A=pi*(r^2-(r-t)^2)*ones(ne,1);



k_i=cell(ne,1);
I_temp=eye(6);
seq_chg=I_temp([1 2 4 5 3 6],:);
for i=1:ne
    k_i{i}=T_ei{i}*seq_chg*[E(i)*A(i)/l(i) 0 0 -E(i)*A(i)/l(i) 0 0
     0 12*E(i)*I(i)*l(i)^-3 6*E(i)*I(i)*l(i)^-2 0 -12*E(i)*I(i)*l(i)^-3 6*E(i)*I(i)*l(i)^-2;
     0 6*E(i)*I(i)*l(i)^-2 4*E(i)*I(i)*l(i)^-1 0 -6*E(i)*I(i)*l(i)^-2 2*E(i)*I(i)*l(i)^-1;
    -E(i)*A(i)/l(i) 0 0 E(i)*A(i)/l(i) 0 0
    0 -12*E(i)*I(i)*l(i)^-3 -6*E(i)*I(i)*l(i)^-2 0 12*E(i)*I(i)*l(i)^-3 -6*E(i)*I(i)*l(i)^-2;
    0 6*E(i)*I(i)*l(i)^-2 2*E(i)*I(i)*l(i)^-1 0 -6*E(i)*I(i)*l(i)^-2 4*E(i)*I(i)*l(i)^-1]*seq_chg'*T_ei{i}';


end

%% equilibrium matrix 1
A_tsgb1=E_qa'*E_e';%*blkdiag(k_i{:});


%% equilibrium matrix 2
A_tsgb2=zeros(3*ne,6*ne);
for i=1:ne
A_tsgb2(3*i-2:3*i,6*i-5:6*i)=[1 0 1 0 0 0;...
                      0 1 0 1 0 0;...
                      0 0 -dyi{i} dxi{i} 1 1];
end
%SVD of equilibrium matrix
A_tsgb=[A_tsgb1;A_tsgb2];
[U1,U2,V1,V2,S1]=tenseg_svd(A_tsgb);


V2_loc=blkdiag(T_ei{:})'*V2;