function d_tuple2=root(xyz) % xyz:直径
global ne w L f E 
%% 考虑外力F=1
x=1;y=6;
for i=1:ne
    u{i,:}=w(x:y,:);
    x=x+6;y=y+6;
end
for i=1:ne
    %% 稳定性计算
    % 各种参数
    u_i=cell2mat(u(i));
    M_max=[u_i(3,:);u_i(6,:)];
    N=max(u_i(1,:),u_i(4,:)); % 所计算构件范围内轴心压力设计值,压力
    f=235e6; % 材料抗拉强度设计值
    E=210e9; % 弹性模量 
    M1=max(abs(M_max(1,:)),abs(M_max(2,:)));M2=min(abs(M_max(1,:)),abs(M_max(2,:)));
    beta_mx=0.6+0.4*(M2/M1);% 无侧移框架柱和两端支承的构件，无横向荷载作用时计算公式，端弯矩abs(M1)>=abs(M2)
    M_x=M1; % M_x所计算构件范围内的最大弯矩设计值
    gamma_x=1.15; % 截面塑形发展系数
        % eta=[0.7]; % 截面影响系数，闭口截面0.7，其它截面1.0
        % beta_tx=0.65+0.35*(M2/M1); % 等效弯矩系数，两端支承的构件段取其中央 1/3 范围内的最大弯矩与全段最大弯矩之比，但不小于0.5;悬臂段取1.0
        % phi_b=[1.0]; % 均匀弯曲的受弯构件整体稳定系数，对闭口截面取1.0   
    % 筛选出受压构件
    press=u_i(4,:);
    if press<0
        % 平面内稳定计算公式
        L_i=L(i);
        d3=fsolve(@(xyzx)(N/(phi_x(xyzx,L_i)*0.25*pi*(xyzx^2-(0.8*xyzx)^2)*f)+beta_mx*M_x/(gamma_x*W_1x(xyzx)*(1-0.8*N/N_Ex(xyzx,L_i))*f)-1),double(xyz(i)));
        %%%%%%%%%%%
        d3_p=d3(find(d3>0));d3_p=real(d3_p);%d3_p=d3_p(abs(imag(d3_p))<eps(d3_p));
        %%%%%%%%%%%
        if isempty(d3_p)
            d_tuple2{i,:}=0; 
        else
            d_tuple2{i,:}=min(d3_p);  
        end
    else
        d_tuple2{i,:}=0; 
    end
    % d_mi=cell2mat(d_tuple2)
end
end