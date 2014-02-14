
close all;
length = 6000;
%Gravity free fall
grav = [0,-9.82,0];
%Step size
delta_t = 0.01;
%Time vector for plotting
t_vector = 0:1:length;
%Position vector with time component [px,py,pz,t]
p_matrix = zeros(3, length+1);
p_matrix(:,1) = [2 , 20, 2];
p_matrix = cat(1, p_matrix , t_vector);
%Velocity vector with time component [vx,vy,vz,t]
v_matrix = zeros(3, length+1);
v_matrix(:,1) = [10, 0, 15];
v_matrix = cat(1, v_matrix , t_vector);
count = 0;
for dt = 1:1:length
    index_new = dt+1;
    index_old = dt;
    %If ball hits ground y-direction
    if p_matrix(2, index_old) <= 0 %&& v_matrix(2, index_old) <= 0 
        %x-pos
        v_matrix(1, index_new) = v_matrix(1, index_old); 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        %y-pos
        p_matrix(2, index_old) = 0;
        v_matrix(2, index_new) = -(0.90*v_matrix(2, index_old)) + grav(2)*delta_t;
        p_matrix(2, index_new) = v_matrix(2, index_new)*delta_t;
        %z-pos
        v_matrix(3, index_new) = v_matrix(3, index_old); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;
        count = count +1;
    %If ball hits wall x = 0 in x-direction
    elseif p_matrix(1, index_old) <= 0
        %x-pos
        p_matrix(1, index_old) = 0;
        v_matrix(1, index_new) = -(0.90*v_matrix(1, index_old)) 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        % y-pos
        v_matrix(2, index_new) = v_matrix(2, index_old) + grav(2)*delta_t; 
        p_matrix(2, index_new) = p_matrix(2, index_old) + v_matrix(2, index_new)*delta_t;
        % z-pos
        v_matrix(3, index_new) = v_matrix(3, index_old); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;
        count = count +1;
        %If ball hits wall x = 300 in x-direction
    elseif p_matrix(1, index_old) >= 300 
        %x-pos
        p_matrix(1, index_old) = 300;
        v_matrix(1, index_new) = -(0.90*v_matrix(1, index_old)); 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        % y-pos
        v_matrix(2, index_new) = v_matrix(2, index_old) + grav(2)*delta_t; 
        p_matrix(2, index_new) = p_matrix(2, index_old) + v_matrix(2, index_new)*delta_t;
        % z-pos
        v_matrix(3, index_new) = v_matrix(3, index_old); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;
        count = count +1;
    %If ball hits wall z = 0 in z-direction
    elseif p_matrix(3, index_old) <= 0
        %x-pos
        v_matrix(1, index_new) = v_matrix(1, index_old); 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        % y-pos
        v_matrix(2, index_new) = v_matrix(2, index_old) + grav(2)*delta_t; 
        p_matrix(2, index_new) = p_matrix(2, index_old) + v_matrix(2, index_new)*delta_t;
        % z-pos
        p_matrix(3, index_old) = 0;
        v_matrix(3, index_new) = -(0.90*v_matrix(3, index_old)); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;
        count = count +1;
         %If ball hits wall z = 300 in z-direction
    elseif p_matrix(3, index_old) >= 300 
        %x-pos
        v_matrix(1, index_new) = v_matrix(1, index_old); 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        % y-pos
        v_matrix(2, index_new) = v_matrix(2, index_old) + grav(2)*delta_t; 
        p_matrix(2, index_new) = p_matrix(2, index_old) + v_matrix(2, index_new)*delta_t;
        % z-pos
        p_matrix(3, index_old) = 300;
        v_matrix(3, index_new) = -(0.90*v_matrix(3, index_old)); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;
        count = count +1;
    %If no collision between ball and wall
    else 
        %x-pos
        v_matrix(1, index_new) = v_matrix(1, index_old); 
        p_matrix(1, index_new) = p_matrix(1, index_old) + v_matrix(1, index_new)*delta_t;
        % y-pos
        v_matrix(2, index_new) = v_matrix(2, index_old) + grav(2)*delta_t;
        p_matrix(2, index_new) = p_matrix(2, index_old) + v_matrix(2, index_new)*delta_t;
        % z-pos
        v_matrix(3, index_new) = v_matrix(3, index_old); 
        p_matrix(3, index_new) = p_matrix(3, index_old) + v_matrix(3, index_new)*delta_t;

    end
    
    if count>50
        v_matrix(1:3,dt:length+1)=0;
        p_matrix(1:3,dt:length+1)=0;

    end
end
%Plot a box works in 2013
%plot(plot::Box(0..300, 0..300, 0..300, Filled = FALSE,
%               LineColor = RGB::Black))

% Plot for position in (x, y)
% plot3(p_matrix(4,:), p_matrix(1,:), p_matrix(2,:));
% title('Position of Bouncing Ball 3D')
% xlabel('Time (t)')
% ylabel('Position (x)')
% zlabel('Position (y)')
% grid on
% % Plot for position in (y, z)
% figure;
% plot3(p_matrix(4,:), p_matrix(2,:), p_matrix(3,:), 'g'); 
% axis square
% title('Position of Bouncing Ball 3D')
% xlabel('Time (t)')
% ylabel('Position (y)')
% zlabel('Position (z)')
% grid on
% % Plot for velocity in (x,y)
% figure;
% plot3(p_matrix(4,:), v_matrix(1,:), v_matrix(2,:), 'r'); 
% axis square
% title('Velocity of Bouncing Ball 3D')
% xlabel('Time (t)')
% ylabel('Velocity (x)')
% zlabel('Velocity (y)')
% grid on

% Plot for position in (x,y,z) % plottar y uppåt pga gravitationen
figure;
plotcube([300 300 30],[0 0 0],.3,[1.0 .73 1.0]); 
hold on
plot3(p_matrix(1,:), p_matrix(3,:), p_matrix(2,:), 'b'); 
axis square
title('Position of Bouncing Ball 3D')

xlabel('Position (x)')
ylabel('Position (z)')
zlabel('Position (y)')
grid on

% Plot for velocity in (x,y,z)
% figure;
% figure;
% plotcube([300 300 30],[0 0 0],.3,[1.0 .73 1.0]); 
% hold on
% plot3(v_matrix(1,:), v_matrix(3,:), v_matrix(2,:), 'm'); 
% axis square
% title('Velocity of Bouncing Ball 3D')
% 
% xlabel('Velocity (x)')
% ylabel('Velocity (z)')
% zlabel('Velocity (y)')
% grid on

% Plot for velocity in (x,t)
figure;
subplot(1,3,1), plot(v_matrix(4,:), v_matrix(1,:), 'g'); 
axis square
title('Velocity of Bouncing Ball 3D along x-axis')
xlabel('Time (t)')
ylabel('Velocity (x)')
grid on

% Plot for velocity in (y,t)

subplot(1,3,2), plot(v_matrix(4,:), v_matrix(2,:), 'b'); 
axis square
title('Velocity of Bouncing Ball 3D along y-axis')
xlabel('Time (t)')
ylabel('Velocity (y)')
grid on

% Plot for velocity in (z,t)

subplot(1,3,3), plot(v_matrix(4,:), v_matrix(3,:), 'black'); 
axis square
title('Velocity of Bouncing Ball 3D along z-axis')
xlabel('Time (t)')
ylabel('Velocity (z)')
grid on
