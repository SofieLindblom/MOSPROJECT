close all;
%Gravity free fall
grav = -0.0982;
%Step size
delta_t = 1;
%Time vector for plotting
t_vector = 0:1:2000;
%Store ball's position for corresponding t
p_vector = [2000];
%Stroe ball's velocity for corresponding t
v_vector = [0];
%counter
count=0;


for dt = 1:1:2000
    index_new = dt+1;
    index_old = dt;
    %If ball hits the floor
    if p_vector(index_old) <= 0 && v_vector(index_old) <= 0
        p_vector(index_old) = 0;
        v_vector(index_new) = -(0.6*v_vector(index_old));
        p_vector(index_new) = v_vector(index_new)*delta_t;
        count = count+1;
    else
        v_vector(index_new) = v_vector(index_old) + grav*delta_t; 
        p_vector(index_new) = p_vector(index_old) + v_vector(index_old)*delta_t;

    end
    %Avoid infinite amount of bouncing
    if count>10
        p_vector(dt:2001) = 0;
        v_vector(dt:2001) = 0;
        break;
    end
end

plot(t_vector, p_vector);
title('Position of Bouncing Ball 1D')
xlabel('Time (t)')
ylabel('Position (y)')
%%hold on;
figure;
plot(t_vector, v_vector, 'r');
title('Velocity of Bouncing Ball 1D')
xlabel('Time (t)')
ylabel('Velocity (v)')