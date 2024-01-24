clc; clear;

Q = 1e5;%throughput
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
wk=wk';

%LOS
C = 10;%环境参数
D = 0.6;
k_los = 0.2;%附加因子
a = 2.3;%路径损失的指数

P0 = (f_p/8)*p_air*s*A_UAV*V_BLADE^3*R_UAV^3;%P0 constant
Pi = ((1+k)*Wei^(3/2))/sqrt(2*p_air*A_UAV);%Pi constant
Ph = P0+Pi;%Ph

%%%%%%%%%%%%%%%%%%%%%%%Algorith1

qk = [0,390,500,400,800;
      0,100,500,900,800];%initial location

pklos=1/(1+C*exp(-D*(90-C)));
gama_k = (pklos+(1-pklos)*k_los)*gama_0; %cal gama_K

v_mr = 18;%Vmr
E0_new = P0*(1/v_mr+3*v_mr/U_tip^2)+Pi*(sqrt(sqrt(v_mr^(-4)+1/(4*v0^4))-1/(2*v0^2)))+1/2*d0*p_air*s*A_UAV*v_mr^2;%E0*

qk_last=qk';

beta_k=zeros(K_GN);
current_rate=zeros(K_GN);
% pklos=zeros(K_GN);
% gama_k=zeros(K_GN);
last_opt=0;
cvx_clear;
for cy=1:10

% for i=1:K_GN
%     pklos(i)=1/(1+C*exp(-D*((180/pi)*asin(H/sqrt(H^2+norm(qk_last(i+1,:)-wk(i,:))^2))-C)));
%     gama_k(i)= (pklos(i)+(1-pklos(i))*k_los)*gama_0; %cal gama_K
% end

% for i=1:K_GN
%     beta_k(i)=(gama_k(i)*(a/2)*(log2(exp(1))))/((H^2+norm(qk_last(i+1,:)-wk(i,:))^2)*(((H^2+norm(qk_last(i+1,:)-wk(i,:))^2))^(a/2)+gama_k(i)));
%     current_rate(i)=log2(1+gama_k(i)/((H^2+norm(qk_last(i+1,:)-wk(i,:))^2)^(a/2)));
% 
% end

for i=1:K_GN
    beta_k(i)=(gama_k*(a/2)*(log2(exp(1))))/((H^2+norm(qk_last(i+1,:)-wk(i,:))^2)*(((H^2+norm(qk_last(i+1,:)-wk(i,:))^2))^(a/2)+gama_k));
    current_rate(i)=log2(1+gama_k/((H^2+norm(qk_last(i+1,:)-wk(i,:))^2)^(a/2)));
    
end
     cvx_begin

        variable Dtr(1)
        variable yita_k(K_GN) nonnegative
        variable qk(K_GN+2,2)
        expression temp(K_GN) 
            for i=1:K_GN
                temp(i) = Qk.*inv_pos(yita_k(i));
            end

        OBJ = E0_new*Dtr+(Ph+Pc)*sum(temp);

        minimize(OBJ);

        subject to
        
        qk(1,:)==[0,0];
        qk(K_GN+2,:)==[800,800];

        sum_square_abs(qk(2,:)-qk(1,:))+sum_square_abs(qk(3,:)-qk(2,:))+sum_square_abs(qk(4,:)-qk(3,:))+sum_square_abs(qk(5,:)-qk(4,:))<=Dtr;
        

        for i=1 :K_GN
            yita_k(i)<=current_rate(i)-beta_k(i)*(sum_square_abs(qk(i+1,:)-wk(i,:))-norm(qk_last(i+1,:)-wk(i,:))^2);  
        end

    cvx_end

    %update
    qk_last=qk;

     %stopping condition
    if norm(cvx_optval-last_opt)/cvx_optval <= 0.1
        break
    else
        last_opt = cvx_optval;
    end

end
hold on
plot(wk(:,1),wk(:,2),'ks');
hold on
plot(qk(:,1),qk(:,2),...
        'Color','r','LineWidth',1),title('Q=1000Mbits');

save("A1_1000v.mat","qk");