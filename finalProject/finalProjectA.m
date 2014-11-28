close all

object_sunverlet
figure()
subplot(4,1,1)
plot3(A(:,1),A(:,2),A(:,3),'-b')
hold on

object_earthverlet
plot3(A(:,1),A(:,2),A(:,3),'-r')
xlabel('x')
ylabel('y')
zlabel('z')

title('Verlet: Position of earth-sun-like system in three dimensions in 80 years, dt = 1 month')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)



% rk43D
% 
% figure()
% subplot(4,1,1)
% plot3(A(:,2),A(:,3),A(:,4),'-r')
% hold on
% plot3(A(:,7),A(:,8),A(:,9),'-b')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend('sun','earth')
% title('RK4: Position of earth-sun-like system in three dimensions in 4 years, dt = 1 week')
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 100)
% 
Verletenergy

subplot(4,1,2)
plot(A(:,1),A(:,2),'-r')
xlabel('time[yr]')
ylabel('energy')
title('Verlet: Energy of earth-sun-like system for 80 years, dt = 1 month')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)


subplot(4,1,3)

object_sunrk4
plot3(A(:,1),A(:,2),A(:,3),'-b')
hold on

object_earthrk4
plot3(A(:,1),A(:,2),A(:,3),'-r')
xlabel('x')
ylabel('y')
zlabel('z')

title('RK4: Position of earth-sun-like system in three dimensions in 80 years, dt = 1 month')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)

RK4energy
subplot(4,1,4)
plot(A(:,1),A(:,2),'-r')
xlabel('time[yr]')
ylabel('energy')
title('RK4: Energy of earth-sun-like system for 80 years, dt = 1 month')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)

