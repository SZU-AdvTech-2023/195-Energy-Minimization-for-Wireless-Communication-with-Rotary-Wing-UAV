clc;
clear;

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

P0 = (f_p/8)*p_air*s*A_UAV*V_BLADE^3*R_UAV^3;%P0 constant
Pi = ((1+k)*Wei^(3/2))/sqrt(2*p_air*A_UAV);%Pi constant
            
P=zeros(1,31);

Parasite=zeros(1,31);

Blade_p=zeros(1,31);

induce=zeros(1,31);

 for i=0:30
     P(:,i+1)=P0*(1+3*i^2/U_tip^2)+Pi*sqrt(sqrt(1+i^4/(4*v0^4))-i^2/(2*v0^2))+1/2*d0*p_air*s*A_UAV*i^3;
 end

for i=0:30
    Blade_p(:,i+1)=P0*(1+3*i^2/U_tip^2);
end

for i=0:30
    induce(:,i+1)=Pi*sqrt(sqrt(1+i^4/(4*v0^4))-i^2/(2*v0^2));
end

for i=0:30
    Parasite(:,i+1)=1/2*d0*p_air*s*A_UAV*i^3;
end


%draw
plot(0:30,P(:,1:31),'red-');
hold on;
plot(0:30,Blade_p(:,1:31),'blue-.');
plot(0:30,induce(:,1:31),'k--');
plot(0:30,Parasite(:,1:31),'k:*');
hold off;


% function Pv = cal_Pv(V)
% 
%     Pv=P0*(1+3*V^2/U_tip^2)+Pi*(sqrt(1+V^4/(4*v0^2))-V^2/(2*v0^2))^(1/2)+1/2*d0*p_air*A_UAV*V^3;
% 
% end