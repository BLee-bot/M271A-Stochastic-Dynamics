% MAE 271A Calibration of an Accelerometer Using GPS
clear all
%% Task1 True Model

w=0.1; % Frequency rad/s
td1=0.02; % Sampling time 50Hz 
V0=100; % Initial Velocity
P0=0; % Initial Position
t=0:td1:30; % Discrete time step
a=10;

A=a*sin(w*t); % Actual acceleration
V=V0+(a/w)-(a/w)*cos(w*t); % True velocity
P=P0+(V0+(a/w))*t-(a/w^2)*sin(w*t); % True Position


% figure(1)
% subplot(3,1,1);
% plot(t,A,'r');
% title('True Acceleration');
% grid on
% xlabel('time(sec)');
% ylabel('acceleration(m/s^{2})')
% subplot(3,1,2);
% plot(t,V,'b')
% title('True Velocity');
% grid on
% xlabel('time(sec)');
% ylabel('velocity(m/s)')
% subplot(3,1,3);
% plot(t,P,'k')
% title('True Position');
% grid on
% xlabel('time(sec)');
% ylabel('positon(m)')
% sgtitle('Truth model')

%% Task2 50Hz acc Model
% Roll
b=0+sqrt(0.01)*randn; % bias roll once
V0=V0+sqrt(1)*randn; % V0
P0=P0+sqrt(100)*randn; % P0
W=0+sqrt(0.0001)*randn; % Gaussia white noise
Ac=[];
Vc=[];
Pc=[];
Vc_temp=V0;
Pc_temp=P0;
for k=1:1501
    Ac_temp=A(1,k)+b+0+sqrt(0.0001)*randn;
    Vc_temp=Vc_temp+Ac_temp*td1;
    Pc_temp=Pc_temp+Vc_temp*td1+Ac_temp*(td1^2)/2;
    Ac=[Ac Ac_temp];
    Vc=[Vc Vc_temp];
    Pc=[Pc Pc_temp];
end

ErrA=A-Ac;
ErrV=V-Vc;
ErrP=P-Pc;



% figure(2)
% subplot(3,1,1);
% plot(t,Ac,'r'); hold on;
% title('Sampling Acceleration');
% grid on
% xlabel('time(sec)'); 
% ylabel('acceleration(m/s^{2})')
% subplot(3,1,2);
% plot(t,Vc,'b'); hold on;
% title('Integraled Velocity');
% grid on
% xlabel('time(sec)'); 
% ylabel('velocity(m/s)')
% subplot(3,1,3);
% plot(t,Pc,'k'); hold on;
% title('Integraled Position');
% grid on
% xlabel('time(sec)');
% ylabel('positon(m)')
% sgtitle('Sampling model X5')
% 
% figure(3)
% subplot(3,1,1);
% plot(t,ErrA,'r'); hold on;
% title('Acceleration Error');
% grid on
% xlabel('time(sec)');
% ylabel('acceleration(m/s^{2})')
% subplot(3,1,2);
% plot(t,ErrV,'b'); hold on;
% title('Velocity Error');
% grid on
% xlabel('time(sec)');
% ylabel('velocity(m/s)')
% subplot(3,1,3);
% plot(t,ErrP,'k'); hold on;
% title('Position Error');
% grid on
% xlabel('time(sec)');
% ylabel('positon(m)')
% sgtitle('Errors between task1 and task2 X5')


e=[];
for c=1:100
    % Task 3 Kalman filter model
    % Setting
    td2=0.5; %2 Hz
    T=0:td2:30;
    Pi=[1 td2 -(td2^2)/2;0 1 -td2; 0 0 1];
    Ga=[-(td2^2)/2;-td2;0];
    H=[1 0 0;0 1 0];
    % Eta
    EtaP=0+sqrt(1)*randn; EtaZ=0+0.04*randn;
    X0=[0;0;0]; P0temp=[sqrt(1)*randn sqrt(0.0016)*randn sqrt(0.01)*randn];
    P0=diag(P0temp);
    Vk=[1 0;0 0.04^2]; W=0.0001;

    Xhat=X0; X0=[X0+P0temp'];
    X_stack=X0;
    Pk=P0;
    Rk_stack=[0;0];
    P_stack=[Pk];
    X_bar_stack=[0;0;0];

    % Kalman Filter 2Hz
    for i=1:length(T)-1 % 63 Steps
        Xbar=Pi*Xhat;
        M=Pi*Pk*Pi'+Ga*W*Ga';
        Kk=Pk*H'*inv(Vk); % Kalman gain
        Pk=inv(inv(M)+H'*inv(Vk)*H);
        Rk=[P(1,25*i)+sqrt(1)*randn-Pc(1,25*i)-Xbar(1,1);V(1,25*i)+0.04*randn-Vc(1,25*i)-Xbar(2,1)]; % Residual
        Xhat=Xbar+Kk*Rk;

        % Stacking data
        Rk_stack=[Rk_stack Rk]; % Residual rk
        P_stack=[P_stack M]; % Covariance P
        X_stack=[X_stack Xhat]; % State X
        X_bar_stack=[X_bar_stack Xbar];
    end

    Pn=[P(1,1)];
    Vn=[V(1,1)];
    Pcn=[Pc(1,1)];
    Vcn=[Vc(1,1)];
    for i=1:60
    Pn=[Pn P(1,25*i)];
    Pcn=[Pcn Pc(1,25*i)];
    Vn=[Vn V(1,25*i)];
    Vcn=[Vcn Vc(1,25*i)];
    end

e=[e;[Pn-Pcn;Vn-Vcn;b*ones(1,61)]-[X_bar_stack(1,:);X_bar_stack(2,:);X_bar_stack(3,:)]];
   
end

e1=sum(e(4:3:end,10))/99
e2=sum(e(5:3:end,10))/99
e3=sum(e(6:3:end,10))/99

Pave=[0 0 0;0 0 0;0 0 0];
for i=1:99
    Pave=Pave+[e(3*i+1,10)-e1;e(3*i+2,10)-e2;e(3*i,10)-e2]*[e(3*i+3,10)-e1;e(3*i-1,10)-e2;e(3*i,10)-e2]';
end

Pave=Pave/98;
L=Pave-P_stack(1:3,28:30);
P_stack(1:3,28:30)
L
    
    
    
    
    
    
    
    
    
    
% figure(4)
% subplot(241)
% plot(T,X_stack(3,:)); hold on;
% title('bias')
% grid on
% subplot(242)
% plot(T,X_stack(3,:)); hold on;
% title('bias')
% grid on
% 
% subplot(243)
% plot(T,P_stack(1,1:3:end)); hold on;
% plot(T,P_stack(1,2:3:end)); 
% plot(T,P_stack(1,3:3:end)); 
% plot(T,P_stack(2,1:3:end)); 
% plot(T,P_stack(2,2:3:end)); 
% plot(T,P_stack(2,3:3:end)); 
% plot(T,P_stack(3,1:3:end)); 
% plot(T,P_stack(3,2:3:end));
% plot(T,P_stack(3,3:3:end)); 
% title('Covariance P')
% subplot(244)
% plot(T,P_stack(1,1:3:end)); hold on;
% plot(T,P_stack(1,2:3:end)); 
% plot(T,P_stack(1,3:3:end)); 
% plot(T,P_stack(2,1:3:end)); 
% plot(T,P_stack(2,2:3:end)); 
% plot(T,P_stack(2,3:3:end)); 
% plot(T,P_stack(3,1:3:end)); 
% plot(T,P_stack(3,2:3:end)); 
% plot(T,P_stack(3,3:3:end)); 
% title('Covariance P')
% subplot(245)
% plot(T,Rk_stack(1,:)); hold on;
% title('Residual P')
% subplot(246)
% plot(T,Rk_stack(1,:)); hold on;
% title('Residual P')
% 
% subplot(247)
% plot(T,Rk_stack(2,:)); hold on;
% title('Residual V')
% subplot(248)
% plot(T,Rk_stack(2,:)); hold on;
% title('Residual V')







% 

% 
% CompP=X_stack(1,:)+Pn;
% CompV=X_stack(2,:)+Vn;
% 
% figure(5)
% plot(t,P,'r'); hold on;
% plot(T,CompP,'ob');
% grid on
% 
% axis([0 30 0 10000]);
% 
% figure(6)
% plot(t,V,'r'); hold on;
% plot(T,CompV,'ob');
% grid on
% 
% axis([0 30 0 1000]);
% 




%Monte Carelo Check Residual
% percent1=0;
% percent2=0;
% 
% Rkmc1=Rkmc(1,:);
% Rkmc2=Rkmc(2,:);
% mean=[0;0];
% 
% for i=1:60
% 
% Rke=H*Pmc(1:3,(3*i-2):(3*i))*H'+Vk;
% 
% if -sqrt(Rke(1,1))<Rkmc1(1,i)<sqrt(Rke(1,1))
%     percent1=percent1+1;
% end
% if -sqrt(Rke(2,2))<Rkmc2(1,i)<sqrt(Rke(2,2))
%     percent2=percent2+1;
% end
% end
% disp(percent1/60)
% disp(percent2/60)



    
    
    

    
    
    
    
    
    
    
    
    







