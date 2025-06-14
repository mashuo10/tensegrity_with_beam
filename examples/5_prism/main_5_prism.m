% @ -0,0 +1,47 @@
%%%% this code calculate the equilibrium matrix, prestress mode of beam-tsg
%%%% 
clc;
clear all;
close all;


%% N C of the structure
% Manually specify node positions of double layer prism.
R=1;
h=2;
q=3;
N=zeros(3,2*q)
N(:,1)=[R,0,0]';
beta=2*pi/q;
T_n=[cos(beta), -sin(beta) 0; sin(beta) cos(beta) 0; 0 0 1];    % To generate bottom node
for i=2:q
N(:,i)=T_n*N(:,i-1);
end
beta=pi*(0.5-1/q); 	% rotation angle
T_n=[cos(beta), -sin(beta) 0; sin(beta) cos(beta) 0; 0 0 1];    % To the Top node
N(:,q+1:end)=T_n*N(:,1:q);
N(3,q+1:end)=h;
n=N(:);

% Manually specify connectivity indices.
C_s_in = [[1:q],[q+1:2*q],[1:q];...
    [2:q,1],[q+2:2*q,q+1],[q+1:2*q]]';  % This is indicating that string connection; Bottom string Top string the diagonal string

C_b_in = [[1:q];[q+2:2*q,q+1]]';  % Similarly, this is saying bar 1 connects node 1 to node 2,

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
title('beam tensegrity prism');
tenseg_plot(N,C,[]);

%% revise C and N
% add middle node in beam
rate=0.8;

B=N*C_b';
N_add=N(:,1:q)+rate*B;
N=[N,N_add];               % New nodal coordinat
n=N(:);

% Manually specify connectivity indices.
C_s_in = [[1:q],[q+1:2*q],[1:q];...
    [2:q,1],[q+2:2*q,q+1],[3*q,2*q+1:3*q-1]]';  % This is indicating that string connection; Bottom string Top string the diagonal string

C_b_in = [[1:q,2*q+1:3*q];...
    [2*q+1:3*q],[q+2:2*q,q+1]]';  % Similarly, this is saying bar 1 connects node 1 to node 2,

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
title('beam tensegrity prism');

%% rigid/pin connection of members ，rotation  Relationship
   

cnct=cell(nn,1); %connection of rigid for 1; pin for 0; otherwise, a cell
for i=[1:2*q]
    cnct{i}=0;
end
for i=2*q+1:3*q
    cnct{i}={[1 2],3};
end


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
        E_ri{i}=kron(temp,eye(1));  
    else if cnct{i}==0
        E_ri{i}=eye(n_m(i));
    else cnct{i}==1
            E_ri{i}=kron(ones(n_m(i),1),eye(1));    
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
E_ei{i}=blkdiag(kron(C_bar{i},eye(3)),kron(E_eri{i},eye(3)));
end
E_e=cell2mat(E_ei);             % 
%% Boundary constraints
% node constraints
pinned_X=[]; pinned_Y=[]; pinned_Z=[]; 
[E_na,E_nb,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
% 
% [E_na,E_nb,a,b]=tenseg_boundary_2D(pinned_X,pinned_Y,nn);
% rotation constraints
ro_const=cell(nn,1);       %all constraint for -1; all free for 0; otherwise, a cell vector
 
for i=[1:3*q]
    ro_const{i}=0;
end

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

E_qa=blkdiag(E_na,kron(blkdiag(E_rai{:}),eye(3)));           % relation between free DOF and minimal DOF
E_qb=blkdiag(E_nb,kron(blkdiag(E_rbi{:}),eye(3)));          % relation between constraint DOF and minimal DOF
%% rotation matrix 
T_i=cell(ne,1);
T_ei=cell(ne,1);

for i=1:ne
xx=[1 0 0]*kron(C(i,:),eye(3))*n;
yy=[0 1 0]*kron(C(i,:),eye(3))*n;
zz=[0 0 1]*kron(C(i,:),eye(3))*n;
if (abs(xx)<1e-5)&&(abs(yy)<1e-5)
T_i{i}=[0 0 1;1 0 0;0 1 0]';
else
temp=[xx yy zz; -yy xx 0;-xx*zz,-yy*zz,xx^2+yy^2]';
T_i{i}=temp/diag(sqrt(sum(temp.^2)));
end
T_ei{i}=kron(eye(4),T_i{i});          
end
%% stiffness matrix
E=1e5*ones(ne,1);
mue=0.3;
G=E/2/(1+mue);
r=0.1;
t=0.01;
Iy=pi/4*(r^4-(r-t)^4)*ones(ne,1);
Iz=pi/4*(r^4-(r-t)^4)*ones(ne,1);
Jk=pi/2*(r^4-(r-t)^4)*ones(ne,1);
A=pi*(r^2-(r-t)^2)*ones(ne,1);



k_i=cell(ne,1);         % Stiffness metrics in global frame
k_e_i=cell(ne,1);         % Stiffness metrics in local frame
I_temp=eye(4);
seq_chg=kron(I_temp(:,[1 3 2 4]),eye(3));
for i=1:ne
    k_e_i{i}=[E(i)*A(i)/l(i) 0 0 0 0 0 -E(i)*A(i)/l(i) 0 0 0 0 0;
    0 12*E(i)*Iz(i)/l(i)^3 0 0 0 6*E(i)*Iz(i)/l(i)^2 0 -12*E(i)*Iz(i)/l(i)^3 0 0 0 6*E(i)*Iz(i)/l(i)^2;
    0 0 12*E(i)*Iy(i)/l(i)^3 0 6*E(i)*Iy(i)/l(i)^2 0 0 0 -12*E(i)*Iy(i)/l(i)^3 0 6*E(i)*Iy(i)/l(i)^2 0;
    0 0 0 G(i)*Jk(i)/l(i) 0 0  0 0 0 -G(i)*Jk(i)/l(i) 0 0 ;
    0 0 6*E(i)*Iy(i)/l(i)^2 0 4*E(i)*Iy(i)/l(i) 0 0 0 -6*E(i)*Iy(i)/l(i)^2 0 2*E(i)*Iy(i)/l(i) 0;
    0 6*E(i)*Iz(i)/l(i)^2 0 0 0 4*E(i)*Iz(i)/l(i) 0 -6*E(i)*Iz(i)/l(i)^2 0 0 0 2*E(i)*Iz(i)/l(i);
    -E(i)*A(i)/l(i) 0 0 0 0 0 E(i)*A(i)/l(i) 0 0 0 0 0;
    0 -12*E(i)*Iz(i)/l(i)^3 0 0 0 -6*E(i)*Iz(i)/l(i)^2 0 12*E(i)*Iz(i)/l(i)^3 0 0 0 -6*E(i)*Iz(i)/l(i)^2;
    0 0 -12*E(i)*Iy(i)/l(i)^3 0 -6*E(i)*Iy(i)/l(i)^2 0 0 0 12*E(i)*Iy(i)/l(i)^3 0 -6*E(i)*Iy(i)/l(i)^2 0;
    0 0 0 -G(i)*Jk(i)/l(i) 0 0  0 0 0 G(i)*Jk(i)/l(i) 0 0 ;
    0 0 6*E(i)*Iy(i)/l(i)^2 0 2*E(i)*Iy(i)/l(i) 0 0 0 -6*E(i)*Iy(i)/l(i)^2 0 4*E(i)*Iy(i)/l(i) 0;
    0 6*E(i)*Iz(i)/l(i)^2 0 0 0 2*E(i)*Iz(i)/l(i) 0 -6*E(i)*Iz(i)/l(i)^2 0 0 0 4*E(i)*Iz(i)/l(i)];
    k_i{i}=T_ei{i}*k_e_i{i}*T_ei{i}';
end

%% equilibrium matrix 1
A_tsgb1=E_qa'*E_e'*kron(eye(ne),seq_chg')*blkdiag(T_ei{:});%*blkdiag(k_i{:});


%% equilibrium matrix 2


A_tsgb2=zeros(6*ne,12*ne);

for i=1:ne
A_tsgb2(6*i-5:6*i,12*i-11:12*i)=[1 0 0 0 0 0 1 0 0 0 0 0;...
                      0 1 0 0 0 0 0 1 0 0 0 0;...
                      0 0 1 0 0 0 0 0 1 0 0 0
                      0 0 0 1 0 0 0 0 0 1 0 0
                      0 0 0 0 1 0 0 0 l(i) 0 1 0
                      0 0 0 0 0 1 0 l(i) 0 0 0 1];
end

%% SVD of equilibrium matrix
A_tsgb=[A_tsgb1;A_tsgb2];
[U1,U2,V1,V2,S1]=tenseg_svd(A_tsgb);        %V2 Is the   Prestress. mode. In global coordinate.

V2_loc=V2;                % This is the stress in local coordinated.
% V2_loc=V1(:,end);                % This is the stress in local coordinated.

%%  Plot the stress in self equilibrium tsgb.

strut_data.T_i=T_i;
strut_data.C_bar=C_bar;

if isfield(strut_data,'displs')
strut_data=rmfield(strut_data,'displs');
end

mod_c=0.5*max(l);           % Change the plot according to member length 
for i=1:size(V2_loc,2)          
strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[0 0 0 1 0 0])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment x-',num2str(i)]);

strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[0 0 0 0 1 0])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 3];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment y-',num2str(i)]);

strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[0 0 0 0 0 1])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment z-',num2str(i)]);

strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[1 0 0 0 0 0])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force x-',num2str(i)]);

strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[0 1 0 0 0 0])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force y-',num2str(i)]);

strut_data.stress=mod_c*kron(eye(ne),[kron(eye(2),[0 0 1 0 0 0])])*round(V2_loc(:,i),3);     % Torque X.
strut_data.seq_plot=[1 3];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force z-',num2str(i)]);
end





%% Group/Clustered information 
%generate group index

gr={};     % number of elements in one group

Gp=kron(tenseg_str_gp(gr,C),eye(12));    %generate group matrix
%% prestress design
%external force in equilibrium design
w0=zeros(size(A_tsgb,1),1);

%prestress design
index_gp=[2*q*12+7];                   % number of groups with designed force; represent the bottom string axial force
fd=1e2;                        % force in bar is given as -1000

I=eye(size(Gp,2));
e_d=I(:,index_gp);        % e_d is the matrix to select group of member with designed force
z=(e_d'*V2)\(fd-e_d'*(pinv(A_tsgb*Gp)*w0));   %self-stress coefficient

fe_c=pinv(A_tsgb*Gp)*w0+V2*z;

fe=Gp*fe_c;



%%  Plot the stress in self equilibrium tsgb.



for i=1:size(fe,2)          
strut_data.stress=kron(eye(ne),[kron(eye(2),[0 0 0 1 0 0])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment x-',num2str(i)]);

strut_data.stress=kron(eye(ne),[kron(eye(2),[0 0 0 0 1 0])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 3];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment y-',num2str(i)]);

strut_data.stress=kron(eye(ne),[kron(eye(2),[0 0 0 0 0 1])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Moment z-',num2str(i)]);

strut_data.stress=kron(eye(ne),[kron(eye(2),[1 0 0 0 0 0])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force x-',num2str(i)]);

strut_data.stress=kron(eye(ne),[kron(eye(2),[0 1 0 0 0 0])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 2];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force y-',num2str(i)]);

strut_data.stress=kron(eye(ne),[kron(eye(2),[0 0 1 0 0 0])])*round(fe(:,i),3);     % Torque X.
strut_data.seq_plot=[1 3];  %plot in xy plane
tenseg_plot_stress_3d(N,C_b,C_s,[],[],[],[],[],strut_data);
title(['Force z-',num2str(i)]);
end













return;

%% mechanism mode
B_tb1=A_tsgb1';
B_tb2=zeros(3*ne,6*ne);

for i=1:ne
% B_tb2(3*i-2:3*i,6*i-5:6*i)=[1 0 0 -1 0 0;...
%                       0 1 0 0 -1 1/l(i) ;...
%                       0 0 1 0 0 -1];

B_tb2(3*i-2:3*i,6*i-5:6*i)=[1 0 0 -1 0 0;...
                      0 -2 -1/l(i) 0 2 -1/l(i) ;...
                      0 0 1 0 0 -1];
end

B_tb=B_tb2*B_tb1;

[U1,U2,V1,V2,S1]=tenseg_svd(B_tb);        %V2 Is the   Prestress. mode. In global coordinate.

V2_loc=V2;  

% V2=[V1(:,end),V2]
% B_tb*V2
%% Plot mechanism mode(use countor plot)
% Plot the structure to make sure it looks right
if isfield(strut_data,'stress')
strut_data=rmfield(strut_data,'stress');
end

for i=1:size(V2,2)  
    
fig=figure
tenseg_plot_dash(N,C_b,C_s,fig);
title('scissor hinge with cables');
N_d=[reshape(E_na*V2(1:size(E_na,2),i),2,[]);zeros(1,nn)];
strut_data.displs=sqrt(sum(N_d.^2)');

N_motion=N+N_d;
tenseg_plot_stress(N_motion,C_b,C_s,fig,[],[],[],[],strut_data);

end
%% plot zero state configuration %%%%% To be finished

displacemt=blkdiag(k_e_i{:})\round(V2_loc(:,1),3)
k_e_i{1}* lsqminnorm(k_e_i{1},V2_loc(1:6,1))-V2_loc(1:6,1)
dis=lsqminnorm(k_e_i{1},V2_loc(1:6,1));

%% Tangent stiffness matrix(only material stiffness)
Ktaa_e=E_qa'*E_e'*kron(eye(ne),seq_chg)*blkdiag(k_i{:})*kron(eye(ne),seq_chg')*E_e*E_qa;
Ktaa_e-Ktaa_e'
[K_mode,D1] = eig(Ktaa_e);         % eigenvalue of tangent stiffness matrix
k=diag(D1);  
[k_sort,I]=sort(k);
K_mode_sort=K_mode(:,I);
% plot the mode eigenvalue
k_sort=round(k_sort,6)
figure
% plot(1:size(D1,2),k_sort,'k-o','linewidth',1.5);
semilogy(1:size(D1,2),k_sort,'k-o','linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Order of Vibration Mode','fontsize',18,'Interpreter','latex');
ylabel('Frequency (Hz)','fontsize',18,'Interpreter','latex');
% grid on;

%plot mode shapes
if isfield(strut_data,'stress')
strut_data=rmfield(strut_data,'stress');
end

for i=1:8
    
fig=figure
title=({['Mode',num2str(i)];['\lambda=',num2str(k_sort(i),'%.2f'),'N/m']});
    
tenseg_plot_dash(N,C_b,C_s,fig);
% title('scissor hinge with cables');
N_d=[reshape(E_na*K_mode_sort(1:size(E_na,2),i),2,[]);zeros(1,nn)];
strut_data.displs=sqrt(sum(N_d.^2)');

N_motion=N+N_d;
tenseg_plot_stress(N_motion,C_b,C_s,fig,[],[],title,[],strut_data);
end

%% statics analysis To be finished
substep=30;

ind_w=[];w=[];

ind_dqb=[]; dqb0=[];
ind_dl0=[5,7]'; dl0=[]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dqb_t,l0_t,E_qa_new,E_qb_new]=tenseg_load_static_RDT(substep,ind_w,w,ind_dqb,dqb0,ind_dl0,dl0,l0,E_qa,E_qb,gravity,[0;9.8;0],C,mass);
[nq,nqa]=size(E_qa_new);
[~,nqb]=size(E_qb_new);
% modify external force(optional)
% w_t(:,1:substep/2)=w_t(:,2:2:end); w_t(:,substep/2+1:end)=w_t(:,end)*ones(1,substep/2);
% dqb_t(:,1:substep/2)=dqb_t(:,2:2:end); dqb_t(:,substep/2+1:end)=dqb_t(:,end)*ones(1,substep/2);
% l0_t(:,1:substep/2)=l0_t(:,2:2:end); l0_t(:,substep/2+1:end)=l0_t(:,end)*ones(1,substep/2);

% input data
rdm_q=[0.01*rand(size(E_na,1),1);zeros(size(E_sa,1),1)];  rdm_q([1:4,5*level+[1:4]])=0;% initial form with disturbance
q_disturb=q+rdm_q;
data.q=q; data.C=C; data.ne=ne; data.nn=nn; data.E_qa=E_qa_new; data.E_qb=E_qb_new;%data.S=S;
data.E=E; data.A=A; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dqb_t=dqb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.substep=substep;    % substep


data_out=static_solver_RDT(data);
t_t=data_out.t_t;          %member force in every step
q_t=data_out.q_t;          %nodal coordinate in every step
l_t=data_out.l_t;          %member length in every step
K_t_t= data_out.Kt_t; %tangent stiffness of whole struct.

n_t=q_t(1:3*nn,:);          % nodal coordinate n
sld_t=q_t(3*nn+1:end,:);    % sliding distance


