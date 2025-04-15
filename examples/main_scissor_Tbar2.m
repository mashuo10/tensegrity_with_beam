% @ -0,0 +1,47 @@
clc;
close all;


%% N C of the structure
% Manually specify node positions of double layer prism.
N=[0 0 0;1 1 0;2 0 0;1 -1 0;1.2 0.2 0]';  
n=N(:);
N2=N(1:2,:);                 %3D to 2D
n2=N2(:);
% Manually specify connectivity indices.
C_s_in = [1 2;2 3;3 4;4 1];  % This is indicating that string connection
C_b_in = [1 5;2 5;5 3;5 4];  % Similarly, this is saying bar 1 connects node 1 to node 2,

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
%% rigid/pin connection of members ，rotation  Relationship
cnct={0;0;0;0;{[1 3];[2 4]}};       %connection of rigid for 1; pin for 0; otherwise, a cell

% generate the E_nri:relation between rotation angle \theta, and reduced angle
% \theta_r,E_nr: collection of E_nri;E_eri: relation between element DOF
% and all DOF vector
% [E_nri,E_nr,E_eri]=tsgb_trans_dof(cnct,nn,n_m);

% generate the relation between rotation angle \theta, and reduced angle
% \theta_r
E_ri=cell(nn,1);
for i=1:nn
    if iscell(cnct{i})
         temp=zeros(n_m(i),numel(cnct{i}));
        for j=1:numel(cnct{i})
            I_temp=eye(n_m(i));
            temp(:,j)=sum(I_temp(:,cnct{i}{j}),2);            
        end
        E_ri{i}=temp;
    else if cnct{i}==0
        E_ri{i}=eye(n_m(i));
    else cnct{i}==1
            E_ri{i}=ones(n_m(i),1);    
    end
    end
end

E_r=blkdiag(E_ri{:});         % r, for \theta rotation angle

n_cr=size(E_r,2);              % num of \theta reduced

E_eri=cell(ne,1);               % \theta in a member
for i=1:ne
  E_eri{i}=zeros(2,n_cr) ; 

C_temp=C;
C_temp(i,:)=C_temp(i,:)*2;
[~,~,v]=find(C_temp);
row_num1=find(v==-2);
row_num2=find(v==2);
E_eri{i}=E_r([row_num1,row_num2],:);
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
%% E_ei   Relationship metrics for an element
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
        E_rai{i}=eye(size(E_ri{i},2));
        E_rbi{i}=[];
    else if ro_const{i}==-1
             E_rai{i}=[];
            E_rbi{i}=ones(size(E_ri{i},2),1);   
    
 else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         I_temp=eye(size(E_ri{i},2));
         a_temp=ro_const{i};
 b_temp=setdiff((1:size(E_ri{i},2)),a_temp);
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

%% SVD of equilibrium matrix
A_tsgb=[A_tsgb1;A_tsgb2];
[U1,U2,V1,V2,S1]=tenseg_svd(A_tsgb);        %V2 Is the   Prestress. mode. In global coordinate.

V2_loc=blkdiag(T_ei{:})'*V2;                % This is the stress in local coordinated.

%%  Plot the stress in self equilibrium tsgb.
strut_s.T_ei=T_ei;
strut_s.C_bar=C_bar;
for i=1:size(V2_loc,2)          
strut_s.stress=kron(eye(ne),[zeros(2,4),eye(2)])*round(V2_loc(:,i),3);     % Moment.
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['moment-',num2str(i)]);
strut_s.stress=kron(eye(ne),[1 0 0 0 0 0;0 0 1 0 0 0])*round(V2_loc(:,i),3);     % axial force
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['axial force-',num2str(i)]);
strut_s.stress=kron(eye(ne),[0 1 0 0 0 0;0 0 0 1 0 0])*round(V2_loc(:,i),3);     % shear force
tenseg_plot_stress(N,C_b,C_s,[],[],[],[],[],strut_s);
title(['shear force-',num2str(i)]);
end