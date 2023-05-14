close all; clear all; clc; clearvars;
N=100; L=N;    % The length of x coordinate
x1=linspace(-L,0,N+1); x2=linspace(0,L,N+1); dx=x1(2)-x1(1);

dt=0.2*dx.^2; % Time step
J=500; % The number of points on t coordinate
% J=1;
t0=0; for j=1:J+1 t(j)=(j-1)*dt; end % Determining t coordinate
aa=1i*dt/(2*dx^2); % Constant coefficient

%% Nonlinearity coefficients
% %%%%% Fulfilling case
% beta_m1=46/46; beta_p1=46/40; 
% beta_p2=46/24; beta_m2=46/21;
% beta_p3=46/22; beta_m3=46/19;
%%%%%% Broken case
beta_m1=46/70; beta_p1=46/58; 
beta_p2=46/22; beta_m2=46/17;
beta_p3=46/16; beta_m3=46/15;

constr=1/beta_m1+1/beta_m2+1/beta_m3-1/beta_p1-1/beta_p2-1/beta_p3


alpha_m1=sqrt(beta_m1); alpha_p1=sqrt(beta_p1);
alpha_m2=sqrt(beta_m2); alpha_p2=sqrt(beta_p2);
alpha_m3=sqrt(beta_m3); alpha_p3=sqrt(beta_p3);

%% The initial conditions
z01=-L/2; z02=L/2;
k0=1;
sigma=8; A=0.1;
Qm1(:,1)=A/alpha_m1*exp(-1i*k0*x1-(x1-z01).^2/(2*sigma)^2)/sqrt(sqrt(2*pi)*sigma);
Qp1(:,1)=A/alpha_p1*exp(1i*k0*x2-(x2-z02).^2/(2*sigma)^2)/sqrt(sqrt(2*pi)*sigma);
% Qm1(1:L+1,1)=0;

q=zeros(6*N+6,1);
q(1:N+1,1)=Qm1(:,1); q(N+2:2*N+2,1)=Qp1(:,1);
%% Defining the matrices
left=ones(1,6*N+5); center=ones(1,6*N+6);
AA=aa*diag(left,-1)+aa*diag(left,1)+(1-2*aa)*diag(center,0);
BB=-aa*diag(left,-1)-aa*diag(left,1)+(1+2*aa)*diag(center,0);
% -1 and +1 bonds
AA(N+1,N  )=0; AA(N+1,N+1)=0; AA(N+1,N+2)=0;
AA(N+2,N+1)=0; AA(N+2,N+2)=0; AA(N+2,N+3)=0; 
BB(N+1,N  )=0; BB(N+1,N+1)=0; BB(N+1,N+2)=0;
BB(N+2,N+1)=0; BB(N+2,N+2)=0; BB(N+2,N+3)=0;
% -2 and +2 bonds
AA(3*N+3,3*N+2)=0; AA(3*N+3,3*N+3)=0; AA(3*N+3,3*N+4)=0;
AA(3*N+4,3*N+3)=0; AA(3*N+4,3*N+4)=0; AA(3*N+4,3*N+5)=0; 
BB(3*N+3,3*N+2)=0; BB(3*N+3,3*N+3)=0; BB(3*N+3,3*N+4)=0;
BB(3*N+4,3*N+3)=0; BB(3*N+4,3*N+4)=0; BB(3*N+4,3*N+5)=0; 
% -3 and +3 bonds
AA(5*N+5,5*N+4)=0; AA(5*N+5,5*N+5)=0; AA(5*N+5,5*N+6)=0;
AA(5*N+6,5*N+5)=0; AA(5*N+6,5*N+6)=0; AA(5*N+6,5*N+7)=0; 
BB(5*N+5,5*N+4)=0; BB(5*N+5,5*N+5)=0; BB(5*N+5,5*N+6)=0;
BB(5*N+6,5*N+5)=0; BB(5*N+6,5*N+6)=0; BB(5*N+6,5*N+7)=0; 

AA(2*N+3,2*N+2)=0; AA(2*N+2,2*N+3)=0; AA(4*N+5,4*N+4)=0; AA(4*N+4,4*N+5)=0; 
BB(2*N+3,2*N+2)=0; BB(2*N+2,2*N+3)=0; BB(4*N+5,4*N+4)=0; BB(4*N+4,4*N+5)=0; 

% % % VBC  eski
% % AA(N+1,N+1)=alpha_m1; AA(N+1,3*N+3)=-alpha_m2;
% % AA(N+2,N+1)=alpha_m1; AA(N+2,5*N+5)=-alpha_m3;
% % AA(3*N+3,N+1)=1/alpha_m1; AA(3*N+3,3*N+3)=1/alpha_m2; AA(3*N+3,5*N+5)=1/alpha_m3;
% % BB(3*N+3,N+3)=1/alpha_p1; BB(3*N+3,3*N+5)=1/alpha_p2; BB(3*N+3,5*N+7)=1/alpha_p3;
% % 
% % AA(3*N+4,N+2)=alpha_p1; AA(3*N+4,3*N+4)=-alpha_p2;
% % AA(5*N+5,N+2)=alpha_p1; AA(5*N+5,5*N+6)=-alpha_p3;
% % AA(5*N+6,N+2)=1/alpha_p1; AA(5*N+6,3*N+4)=1/alpha_p2; AA(5*N+6,5*N+6)=1/alpha_p3;
% % BB(5*N+6,N  )=1/alpha_m1; BB(5*N+6,3*N+2)=1/alpha_m2; BB(5*N+6,5*N+4)=1/alpha_m3;

%% VBC
AA(N+1,N+1)=alpha_m1; AA(N+1,3*N+3)=-alpha_m2;
AA(N+2,N+1)=alpha_m1; AA(N+2,5*N+5)=-alpha_m3;
AA(3*N+3,N+1)=alpha_m1; AA(3*N+3,N+2)=-alpha_p1;
AA(3*N+4,N+1)=alpha_m1; AA(3*N+4,3*N+4)=-alpha_p2;
AA(5*N+5,N+1)=alpha_m1; AA(5*N+5,5*N+6)=-alpha_p3;
AA(5*N+6,N+1)=1/alpha_m1; AA(5*N+6,3*N+3)=1/alpha_m2; AA(5*N+6,5*N+5)=1/alpha_m3;
AA(5*N+6,N+2)=1/alpha_p1; AA(5*N+6,3*N+4)=1/alpha_p2; AA(5*N+6,5*N+6)=1/alpha_p3;
BB(5*N+6,N+3)=1/alpha_p1; BB(5*N+6,3*N+5)=1/alpha_p2; BB(5*N+6,5*N+7)=1/alpha_p3;
BB(5*N+6,N  )=1/alpha_m1; BB(5*N+6,3*N+2)=1/alpha_m2; BB(5*N+6,5*N+4)=1/alpha_m3;


AA=sparse(AA); BB=sparse(BB);

time=1;
%% The main cycle
for jj=1:J+1
%     %%%%%%% b_{-1} va b_{1}
%     for n=2:2*N+1
%         BB(n,n)=BB(n,n)-1i*dt*conj(q(2*N+2-n,1))*(q(n+1,1)+q(n-1,1));
%     end
%     BB(N+1,N+1)=0; BB(N+2,N+2)=0; 
%     %%%%%%% b_{-2} va b_{2}
%     for n=2*N+4:4*N+3
%         BB(n,n)=BB(n,n)-1i*dt*conj(q(6*N+6-n,1))*(q(n+1,1)+q(n-1,1));
%     end
%     BB(3*N+3,3*N+3)=0; BB(3*N+4,3*N+4)=0;
%     %%%%%%% b_{-3} va b_{3}
%     for n=4*N+6:6*N+5
%         BB(n,n)=BB(n,n)-1i*dt*conj(q(10*N+10-n,1))*(q(n+1,1)+q(n-1,1));
%     end
%     BB(5*N+5,5*N+5)=0; BB(5*N+6,5*N+6)=0;
    
    
    nnn=1;
    %%%%%%% b_{-1} va b_{1}
    for n=2:2*N+1
        q(n,jj)=q(n,jj)*exp(-1i*sqrt(beta_m1*beta_p1)*dt*conj(q(2*N+2-n,jj))*(q(n+nnn,jj)+q(n-nnn,jj))/2);
    end
    %%%%%%% b_{-2} va b_{2}
    for n=2*N+3:4*N+4
        q(n,jj)=q(n,jj)*exp(-1i*sqrt(beta_m2*beta_p2)*dt*conj(q(6*N+6-n,jj))*(q(n+nnn,jj)+q(n-nnn,jj))/2);
    end
    %%%%%%% b_{-3} va b_{3}
    for n=4*N+5:6*N+5
        q(n,jj)=q(n,jj)*exp(-1i*sqrt(beta_m3*beta_p3)*dt*conj(q(10*N+10-n,jj))*(q(n+nnn,jj)+q(n-nnn,jj))/2);
    end

    
    q(:,jj+1)=AA\BB*q(:,jj);
    
    
    %%%%%%% b_{-1} va b_{1}
    for n=2:2*N+1
        q(n,jj+1)=q(n,jj+1)*exp(-1i*sqrt(beta_m1*beta_p1)*dt*conj(q(2*N+2-n,jj+1))*(q(n+nnn,jj+1)+q(n-nnn,jj+1))/2);
    end
    %%%%%%% b_{-2} va b_{2}
    for n=2*N+3:4*N+4
        q(n,jj+1)=q(n,jj+1)*exp(-1i*sqrt(beta_m2*beta_p2)*dt*conj(q(6*N+6-n,jj+1))*(q(n+nnn,jj+1)+q(n-nnn,jj+1))/2);
    end
    %%%%%%% b_{-3} va b_{3}
    for n=4*N+5:6*N+5
        q(n,jj+1)=q(n,jj+1)*exp(-1i*dt*sqrt(beta_m3*beta_p3)*conj(q(10*N+10-n,jj+1))*(q(n+nnn,jj+1)+q(n-nnn,jj+1))/2);
    end
    
    
    
% set(gcf,'Position',[200 0 900 600]);
% subplot(3,2,1)
% plot(x1,abs(q(1:N+1,jj+1)).^2,'r','Linewidth',2); title('Bond 11'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% % ylim([0 0.002])
% 
% subplot(3,2,2)
% plot(x2,abs(q(N+2:2*N+2,jj+1)).^2,'r','Linewidth',2); title('Bond 12'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% % ylim([0 0.002])
% 
% subplot(3,2,3)
% plot(x1,abs(q(2*N+3:3*N+3,jj+1)).^2,'r','Linewidth',2); title('Bond 22'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.06])
% 
% subplot(3,2,4)
% plot(x2,abs(q(3*N+4:4*N+4,jj+1)).^2,'r','Linewidth',2);  title('Bond 21'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.6])
% 
% subplot(3,2,5)
% plot(x1,abs(q(4*N+5:5*N+5,jj+1)).^2,'r','Linewidth',2);  title('Bond 32'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.6])
% 
% subplot(3,2,6)
% plot(x2,abs(q(5*N+6:6*N+6,jj+1)).^2,'r','Linewidth',2);  title('Bond 31'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.6])
% drawnow

if mod(jj,10)==0
        Q(:,time)=q(:,jj+1);
        time=time+1;
        fprintf('%s%i\n','time=',time)
end

fprintf('%s%i\n','j=',jj)
end


% jj=510;
% for jj=1:880
% set(gcf,'Position',[200 0 900 600]);
% subplot(3,2,1)
% plot(x1,abs(q(1:N+1,jj+1)).^2,'r','Linewidth',2); title('Bond 11'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% % ylim([0 0.0015])
% 
% subplot(3,2,2)
% plot(x2,abs(q(N+2:2*N+2,jj+1)).^2,'r','Linewidth',2); title('Bond 12'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% % ylim([0 0.0015])
% 
% subplot(3,2,3)
% plot(x1,abs(q(2*N+3:3*N+3,jj+1)).^2,'r','Linewidth',2); title('Bond 22'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.0015])
% 
% subplot(3,2,4)
% plot(x2,abs(q(3*N+4:4*N+4,jj+1)).^2,'r','Linewidth',2);  title('Bond 21'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.0015])
% 
% subplot(3,2,5)
% plot(x1,abs(q(4*N+5:5*N+5,jj+1)).^2,'r','Linewidth',2);  title('Bond 32'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.0015])
% 
% subplot(3,2,6)
% plot(x2,abs(q(5*N+6:6*N+6,jj+1)).^2,'r','Linewidth',2);  title('Bond 31'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.0015])
% drawnow
% % hold on
% end


% for TT=1:200
% set(gcf,'Position',[200 0 900 600]);
% subplot(3,2,1)
% plot(x1,abs(Q(1:N+1,TT)).^2,'r','Linewidth',2); title('Bond 11'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% ylim([0 1])
% 
% subplot(3,2,2)
% plot(x2,abs(Q(N+2:2*N+2,TT)).^2,'r','Linewidth',2); title('Bond 12'); xlabel('x','Fontsize',14); ylabel('|q_{1}|^2','Fontsize',14);
% ylim([0 1])
% 
% subplot(3,2,3)
% plot(x1,abs(Q(2*N+3:3*N+3,TT)).^2,'r','Linewidth',2); title('Bond 22'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.06])
% 
% subplot(3,2,4)
% plot(x2,abs(Q(3*N+4:4*N+4,TT)).^2,'r','Linewidth',2);  title('Bond 21'); xlabel('x','Fontsize',14); ylabel('|q_{2}|^2','Fontsize',14);
% % ylim([0 0.6])
% 
% subplot(3,2,5)
% plot(x1,abs(Q(4*N+5:5*N+5,TT)).^2,'r','Linewidth',2);  title('Bond 32'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.6])
% 
% subplot(3,2,6)
% plot(x2,abs(Q(5*N+6:6*N+6,TT)).^2,'r','Linewidth',2);  title('Bond 31'); xlabel('x','Fontsize',14); ylabel('|q_{3}|^2','Fontsize',14);
% % ylim([0 0.6])
% drawnow
% end




for hhh=1:1
    %% Plotting
middle=130; final=300;
font=16;
figure('DefaultAxesFontSize',14); set(gcf,'Position',[200 0 600 1600]);
subplot(3,2,1)
plot(x1,abs(q(1:N+1,1)).^2,'r','Linewidth',2); title('b_{-1}','Fontsize',font);
ylabel('|q_{-1}|^2','Fontsize',font); set(gca,'XTickLabel',[]);
ylim([0 0.0008])
hold on
plot(x1,abs(q(1:N+1,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q(1:N+1,final)).^2,'b','Linewidth',2);

subplot(3,2,2)
plot(x2(2:N+1),abs(q(N+3:2*N+2,1)).^2,'r','Linewidth',2); title('b_{1}','Fontsize',font); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
ylim([0 0.0008])
hold on
plot(x2(2:N+1),abs(q(N+3:2*N+2,middle)).^2,'g','Linewidth',2);
hold on
plot(x2(2:N+1),abs(q(N+3:2*N+2,final)).^2,'b','Linewidth',2);

subplot(3,2,3)
plot(x1,abs(q(2*N+3:3*N+3,1)).^2,'r','Linewidth',2); title('b_{-2}','Fontsize',font);
ylabel('|q_{-2}|^2','Fontsize',font); set(gca,'XTickLabel',[]);
ylim([0 0.00025])
hold on
plot(x1,abs(q(2*N+3:3*N+3,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q(2*N+3:3*N+3,final)).^2,'b','Linewidth',2);

subplot(3,2,4)
plot(x2(2:N+1),abs(q(3*N+5:4*N+4,1)).^2,'r','Linewidth',2);  title('b_{2}','Fontsize',font); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
ylim([0 0.00025])
hold on
plot(x2(2:N+1),abs(q(3*N+5:4*N+4,middle)).^2,'g','Linewidth',2);
hold on
plot(x2(2:N+1),abs(q(3*N+5:4*N+4,final)).^2,'b','Linewidth',2);

subplot(3,2,5)
plot(x1,abs(q(4*N+5:5*N+5,1)).^2,'r','Linewidth',2);  title('b_{-3}','Fontsize',font);
ylabel('|q_{-3}|^2','Fontsize',font); xlabel('x','Fontsize',font);
ylim([0 0.0002])
hold on
plot(x1,abs(q(4*N+5:5*N+5,middle)).^2,'g','Linewidth',2);
hold on
plot(x1,abs(q(4*N+5:5*N+5,final)).^2,'b','Linewidth',2); 

subplot(3,2,6)
plot(x2(2:N+1),abs(q(5*N+7:6*N+6,1)).^2,'r','Linewidth',2);  title('b_{3}','Fontsize',font);
xlabel('x','Fontsize',font); set(gca,'YTickLabel',[]);
ylim([0 0.0002])
hold on
plot(x2(2:N+1),abs(q(5*N+7:6*N+6,middle)).^2,'g','Linewidth',2); 
hold on
plot(x2(2:N+1),abs(q(5*N+7:6*N+6,final)).^2,'b','Linewidth',2);
legend('t=0','t=45','t=100')
end


%% Plotting contour plot figure('DefaultAxesFontSize',14); set(gcf,'Position',[200 0 600 1600]);
final=300;
for gg=1:1
figure('DefaultAxesFontSize',14); set(gcf,'Position',[200 0 700 1200]);%set(gcf,'Position',[100 10 600 1200]);
subplot(3,2,1)
mesh(x1,t(1:final),abs(q(1:N+1,1:final)').^2); ylim([0 dt*final]); title('b_{-1}','Fontsize',16); ylabel('t','Fontsize',16);
colorbar; set(gca,'XTickLabel',[]);
view(0,90); 

subplot(3,2,2)
mesh(x2(2:N+1),t(1:final),abs(q(N+3:2*N+2,1:final)').^2); ylim([0 dt*final]); title('b_{1}','Fontsize',16);
colorbar; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
view(0,90)

subplot(3,2,3)
mesh(x1,t(1:final),abs(q(2*N+3:3*N+3,1:final)').^2); ylim([0 dt*final]); title('b_{-2}','Fontsize',16); ylabel('t','Fontsize',16);
colorbar; set(gca,'XTickLabel',[]);
view(0,90)

subplot(3,2,4)
mesh(x2(2:N+1),t(1:final),abs(q(3*N+5:4*N+4,1:final)').^2); ylim([0 dt*final]); title('b_{2}','Fontsize',16); 
colorbar; set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
view(0,90)

subplot(3,2,5)
mesh(x1,t(1:final),abs(q(4*N+5:5*N+5,1:final)').^2); ylim([0 dt*final]); title('b_{-3}','Fontsize',16); xlabel('n','Fontsize',16); ylabel('t','Fontsize',16);
colorbar;
view(0,90)

subplot(3,2,6)
mesh(x2(2:N+1),t(1:final),abs(q(5*N+7:6*N+6,1:final)').^2); ylim([0 dt*final]); title('b_3','Fontsize',16); xlabel('n','Fontsize',16);
colormap(jet(256)); set(gca,'YTickLabel',[]);
colorbar;
view(0,90)
end




    
    