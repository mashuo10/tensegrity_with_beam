function volume = stability_test_func(z)
%% function of optimization_model5_min_csa_all.m
global T_diag V2 A ne w L d_min1 d_min3 d1   
%% 考虑外力F=10000
% w=T_diag*(V2*z+pinv(A)*[0 0 0 0 10000 0 0 10000 0 0 10000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'); % w：局部坐标系下的节点力
% z=pinv(V2)*(pinv(T_diag)*w-pinv(A)*[0 0 0 0 10000 0 0 10000 0 0 10000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]')
[ex_force,~,~,~,~,~,~] = input_var;
w=T_diag*(V2*z+pinv(A)*ex_force); % w：局部坐标系下的节点力
z=pinv(V2)*(pinv(T_diag)*w-pinv(A)*ex_force);
% ww=pinv(A)*[0 0 10000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
x=1;y=6;
% for i=1:ne
%     u{i,:}=w(x:y,:);
%     x=x+6;y=y+6;
% end
%% 截面抵抗矩W:圆：W=pi*d^3/32;d：圆形截面直径;c：材料承受的最大应力
for i=1:ne
    %% 强度计算
    % u_i=cell2mat(u(i));
    u_i=w(x:y,:);x=x+6;y=y+6;
    F_max=[u_i(1,:);u_i(4,:)];
    M_max=[u_i(3,:);u_i(6,:)];
    % 用符号表达式，推荐
    syms d;     % 环径0.8d
    sigma=235e6;   % d：变量；sigma：应力
    W=pi*(d^4-(0.8*d)^4)/(32*d);   % W：圆环截面抗弯抵抗矩
    A_=0.25*pi*(d^2-(0.8*d)^2);     % A_：圆环截面积
    % eqn(1)=sigma*pi*d^3-4*abs(F_max(1,:))*d-32*abs(M_max(1,:))==0;    % eqn1：微分方程：杆端1
    eqn(1)=abs(F_max(1,:)/A_)+abs(M_max(1,:)/W)-sigma;
    d1=double(solve(eqn(1),d));    % d1：微分方程的解
    d1_p=d1(find(d1>0));     % d1_p：找出d1中大于0的数
    % eqn(2)=sigma*pi*d^3-4*abs(F_max(2,:))*d-32*abs(M_max(2,:))==0;    % eqn2：微分方程：杆端2
    eqn(2)=abs(F_max(2,:)/A_)+abs(M_max(2,:)/W)-sigma;
    d2=double(solve(eqn(2),d));    % d2：微分方程的解
    d2_p=d2(find(d2>0));     % d2_p：找出d2中大于0的数
    d_tuple1{i,:}=max(min(d1_p),min(d2_p));    % d_tuple： 各个杆件满足条件圆环外径最小值组成的元组 
end
%%%%%%%%%%
d_min1=real(cell2mat(d_tuple1));   % d_min2：行：表示各个杆件圆环外径最小值
%%%%%%%%%%
%% 稳定性验算
d_tuple2=root(d_min1);   % d3_tuple：各个杆件满足条件圆环外径最小值组成的元组 
d_min2=cell2mat(d_tuple2);   % d_min2：行：表示各个杆件圆环外径最小值
xy=[d_min1,d_min2]';d1=max(xy);    % d1:环形圆杆最终的外径
d_min3=0.25*pi*(d1.^2-(0.8*d1).^2)*1e6;   % d_min3：杆件截面积最小值(mm2)
volume=d_min3*L*1e3; % 体积单位转换为立方毫米
end

