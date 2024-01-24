clc; clear;

%Q = 1e5;%throughput
Qk = 1000;


%UAV参数
Wei = 20;%UAV重量
p_air = 1.225;%air density kg/m^3
R_UAV = 0.4;%UAV扇叶半径
A_UAV = 0.503;%UAV扇叶面积
V_BLADE = 300;%叶片角速度
U_tip = 120;%叶片尖端速度
b = 4;%叶片数量
c = 0.0157;%叶片或翼型弦长度
s = 0.05;%Rotor solidity
S_PF = 0.0151;%Fuselage euivalent flat plate area in m^2
d0 = 0.6;%Fuselage drag ratio
k = 0.1;%Incremental correction factor to induced power
v0 = 4.03;%Mean rotor induced velocity in hover
f_p = 0.012;%Profile drag coefficient


%位置信息
v_max = 30;%UAV最大飞行速度
qi = [0;
      0];%UAV initial location
qf = [800;
      800];%UAV final location

%通信相关参数
H = 100;%UAV的高度
B = 1e6;%带宽
gama_0 = 52.5;%信噪比
Pc = 5;%UAV通信功率
P = 20;%发射功率

%GN参数
K_GN = 3;%GN numbers GN1 = [500;500];%GN1的位置GN2 = [390;900];%GN2的位置GN3 = [390;100];%GN3的位置
wk = [390,500,400;
      100,500,900];%GN location

wk_t=wk';

%LOS
C = 10;%环境参数
D = 0.6;
k_los = 0.2;%附加因子
a = 2.3;%路径损失的指数



P0 = (f_p/8)*p_air*s*A_UAV*V_BLADE^3*R_UAV^3;%P0 constant
Pi = ((1+k)*Wei^(3/2))/sqrt(2*p_air*A_UAV);%Pi constant
Ph = P0+Pi;%Ph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A1 Loc
qk1=[0,390,500,400,800;
      0,100,500,900,800];


%pre-define M and delta_max

delta_max=10;

M=402; %分4段 M+2需要可以整除


%init var

%init v
V_ini=18;

%init gama_k
pklos=1/(1+C*exp(-D*(90-C)));
gama_k = (pklos+(1-pklos)*k_los)*gama_0; %cal gama_K

%init waypoint
qi=[0;0];
qf=[800;800];
                        %有4段
q1_x=linspace(qk1(1,1),qk1(1,2),(M+2)/4);
q1_y=linspace(qk1(2,1),qk1(2,2),(M+2)/4);
q1=[q1_x;q1_y]';
q2_x=linspace(qk1(1,2),qk1(1,3),(M+2)/4);
q2_y=linspace(qk1(2,2),qk1(2,3),(M+2)/4);
q2=[q2_x;q2_y]';
q3_x=linspace(qk1(1,3),qk1(1,4),(M+2)/4);
q3_y=linspace(qk1(2,3),qk1(2,4),(M+2)/4);
q3=[q3_x;q3_y]';
q4_x=linspace(qk1(1,4),qk1(1,5),(M+2)/4);
q4_y=linspace(qk1(2,4),qk1(2,5),(M+2)/4);
q4=[q4_x;q4_y]';
q_last=[q1;q2;q3;q4];

% hold on;
% plot(q_last(:,1),q_last(:,2),'m--');

% initialize the time
T_last = zeros(M+1,1);
for m = 1:M+1
    T_last(m) = norm(q_last(m+1,:)-q_last(m,:))/V_ini;
end




% initialize ym (52)
y_last = zeros(M+1,1);
for m = 1:M+1
    y_last(m) = sqrt(...
        sqrt(T_last(m)^4 + norm(q_last(m+1,:)-q_last(m,:))^4/(4*v0^4))...
        - norm(q_last(m+1,:)-q_last(m,:))^2/(2*v0^2));
end

% initialize tao
tao_last = zeros(M+1,K_GN);
for i=1:K_GN
    for m = 1:M+1
        tao_last(m,i)=T_last(m)/K_GN;
    end
end

% initialize A
A_last = zeros(M+1,K_GN);
for i=1:K_GN
    for m = 1:M+1
        A_last(m,i) = sqrt(tao_last(m,i)*...
            log2(1+gama_k/(H^2 + norm(q_last(m,:) - wk_t(i,:))^2)));
    end
end

%%%%%%%%optimizing by cvx
 l=1;
 max_ite=10;
 last_opt=0;

 cvx_clear;

for ite=1:2
    fprintf('==================== liter = %d ====================\n',l);
    beta_mk=zeros(M+1,K_GN);
    current_rate=zeros(M+1,K_GN);
    for i=1:K_GN
        for m=1:(M+1)
            beta_mk(m,i)=(log2(exp(1))*gama_k*(a/2))/...
                ((H^2+norm(q_last(m,:)-wk_t(i,:))^2)*...
                ((H^2+norm(q_last(m,:)-wk_t(i,:))+gama_k)));
            current_rate(m,i)=log2(1+gama_k/(H^2+norm(q_last(m,:)-wk_t(i,:))));
        end
    end

    cvx_begin
    cvx_solver mosek

    variable q(M+2,2)
    variable T(M+1,1) nonnegative
    variable Ym(M+1,1) nonnegative
    variable A(M+1,K_GN) nonnegative
    variable Y(M+1,1) % delta^3/T^2
    variable tao(M+1,K_GN) nonnegative%tao
    variable S(M+1,1) %Y>=delta^3/T^2
    variable S2(M+1,1) nonnegative %slack variable for convex term T^4/Ym^2
    variable delta_m(M+1,1) %delta_m=||q{m+1}-q{m}||

    expression energy(M+1,K_GN)

    for i=1:K_GN
        for m=1:M+1
            energy(m,i)=P0*(T(m)+(3/U_tip^2)*quad_over_lin(delta_m(m),T(m)))...
                +Pi*Ym(m)...
                +0.5*d0*p_air*s*A_UAV*Y(m)...
                +Pc*tao(m,i);
        end
    end

    minimize sum(sum(energy))

    subject to

    for i=1:K_GN
        sum(A_last(:,i)).^2+2*sum(A_last(:,i)).*(sum(A(:,i))-sum(A_last(:,i)))>=Qk;
        %sum(A_last(:,i).^2+2*A_last(:,i).*(A(:,i)-A_last(:,i)))>=Qk;
        %sum_Al(1,i).^2+2*sum_Al(1,i).*(sum_A(1,i)-sum_Al(1,i))>=Qk;
        %temp_Am(1,i)>=Qk;
    end

    % T_m^4/z_m^2 <= {z_m^l}^2 + 2*z_m^l*(z_m-z_m^l)
    % - || q_{m+1}^l - q_{m}^l ||^2/v_0^2
    % + (2/v_0^2)*(q_{m+1}^l - q_{m}^l)^T(q_{m+1} - q_{m})
    for m=1:M+1
        pow_p(S2(m),2)<=...
            y_last(m)^2+2*y_last(m)*(Ym(m)-y_last(m))...
            - norm(q_last(m+1,:)-q_last(m,:))^2/(v0^2) ...
            + (2/v0^2).*(q_last(m+1,:)-q_last(m,:))'*(q(m+1,:)-q(m,:));
    end

    for m=1:M+1
        S2(m)>=quad_over_lin(T(m),Ym(m));
    end

    % Amk^2/tao_mk <= R(q(m))^(l)
    for i=1:K_GN
        for m=1:M+1
            quad_over_lin(A(m,i),tao(m,i))<=current_rate(m,i)-beta_mk(m,i)*...
                (sum_square_abs(q(m,:)-wk_t(i,:))-norm(q_last(m,:)-wk_t(i,:))^2);
        end
    end


    % %ym>=0
    % for m=1:M+1
    %     Ym(m)>=0;
    % end

    %q1=qi q_m+2=qf
    q(1,:)==qi';
    q(M+2,:)==qf';

    %||q_{m+1}-q_{m}||<=min(delta_max,T_m*V_max)
    for m=1:M+1
        norm(q(m+1,:)-q(m,:))<=delta_max;
        norm(q(m+1,:)-q(m,:))<=T(m)*v_max;
        norm(q(m+1,:)-q(m,:))<=delta_m(m);
    end

    %sum(taomk)<=Tm
    for m=1:M+1
        sum(tao(m,:))<=T(m);
    end

    % %taomk>=0
    % for i=1:K_GN
    %     for m=1:M+1
    %         tao(m,i)>=0;
    %     end
    % end

    % Y(m) >= ||q(m+1,:)-q(m,:)||^3/T(m)^2
    for m = 1:M+1
        {delta_m(m),S(m),Y(m)} == rotated_lorentz(1);
        {S(m),Y(m),delta_m(m)} == rotated_lorentz(1);
    end


    cvx_end

    %update
    q_last=q;
    T_last=T;
    y_last=Ym;
    A_last=A;
    tao_last=tao;
    
    l=l+1;
    %stopping condition
    if norm(cvx_optval-last_opt)/cvx_optval <= 0.1
        break
    else
        last_opt = cvx_optval;
    end

end

hold on;
plot(qk1(1,:),qk1(2,:),'m--');

hold on
plot(wk_t(:,1),wk_t(:,2),'ks');

hold on
plot(q(:,1),q(:,2),...
        'Color','r','LineWidth',1),title('Q = 1000Mbits');
hold off
plot_v=zeros(M+1,1);
for m=1:M+1
    plot_v(m)=norm(q(m+1,:)-q(m,:))/T(m);
end
save("1000v.mat","plot_v","T");