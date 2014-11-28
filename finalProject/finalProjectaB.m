close all

figure()
object_sunverlet_
subplot(2,1,1)
plot3(A(:,1),A(:,2),A(:,3),'ob')
hold on
% subplot(3,1,2)
object_earthverlet_
plot3(A(:,1),A(:,2),A(:,3),'-r')
% 
hold on
% subplot(3,1,3)
object_moonverlet_
plot3(A(:,1),A(:,2),A(:,3),'-g')
% 
xlabel('x[AU]')
ylabel('y[AU]')
zlabel('z[AU]')
legend('Sun', 'Earth', 'Moon')

title('Verlet: Position of earth-sun-moon like system with constant timesteps of one hour')
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 100)

object_sunrk4_
subplot(2,1,2)
plot3(A(:,1),A(:,2),A(:,3),'ob')
hold on
% subplot(3,1,2)
object_earthrk4_
plot3(A(:,1),A(:,2),A(:,3),'-r')
% 
hold on
% subplot(3,1,3)
object_moonrk4_
plot3(A(:,1),A(:,2),A(:,3),'-g')
% 
xlabel('x[AU]')
ylabel('y[AU]')
zlabel('z[AU]')
legend('Sun', 'Earth', 'Moon')
title('RK4: Position of earth-sun-moon like system with constant timesteps of one hour')

figure()
object_sunverletAdapt_
subplot(2,1,1)
plot3(A(:,1),A(:,2),A(:,3),'ob')
hold on
% subplot(3,1,2)
object_earthverletAdapt_
plot3(A(:,1),A(:,2),A(:,3),'-r')
% 
hold on
% subplot(3,1,3)
object_moonverletAdapt_
plot3(A(:,1),A(:,2),A(:,3),'-g')
% 
xlabel('x')
ylabel('y')
zlabel('z')
legend('Sun', 'Earth', 'Moon')
title('Verlet: Position of earth-sun-moon like system with adaptive timesteps')

object_sunrk4Adapt_
subplot(2,1,2)
plot3(A(:,1),A(:,2),A(:,3),'ob')
hold on
% subplot(3,1,2)
object_earthrk4Adapt_
plot3(A(:,1),A(:,2),A(:,3),'-r')
% 
hold on
% subplot(3,1,3)
object_moonrk4Adapt_
plot3(A(:,1),A(:,2),A(:,3),'-g')
% 
xlabel('x')
ylabel('y')
zlabel('z')
legend('Sun', 'Earth', 'Moon')
title('RK4: Position of earth-sun-moon like system with adaptive timesteps')