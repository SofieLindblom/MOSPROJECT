
close all;
%length of simulation
length = 6000;
%ball radius
r = 1;
%Time vector for plotting
t_vector = 0:1:length;
%nr of bounces
count = 0;
%nr of balls
nr_of_balls = 2;

%Ball 1
    %Position vector with time component [px,py,pz,t]
    b1_p_matrix = zeros(3, length+1);
    b1_p_matrix(:,1) = [2 , 20, 2];
    b1_p_matrix = cat(1, b1_p_matrix , t_vector);
    %Velocity vector with time component [vx,vy,vz,t]
    b1_v_matrix = zeros(3, length+1);
    b1_v_matrix(:,1) = [10, 0, 15];
    b1_v_matrix = cat(1, b1_v_matrix , t_vector);
    
%Ball 2
    %Position vector with time component [px,py,pz,t]
    b2_p_matrix = zeros(3, length+1);
    b2_p_matrix(:,1) = [2 , 20, 298];
    b2_p_matrix = cat(1, b2_p_matrix , t_vector);
    %Velocity vector with time component [vx,vy,vz,t]
    b2_v_matrix = zeros(3, length+1);
    b2_v_matrix(:,1) = [10, 0, -15];
    b2_v_matrix = cat(1, b2_v_matrix , t_vector);
    
% Ball vector [bn_p_matrix , bn_v_matrix]
ball_vector = zeros((8*nr_of_balls),(length+1));
ball_vector = [b1_p_matrix; b1_v_matrix; b2_p_matrix; b2_v_matrix];
    
% Main loop
for dt = 1:1:length

for i = 1:1:nr_of_balls % P_temp = ball_vector(((i-1*4)+1):i*4,dt)
%     p_temp = ball_vector(((i-1*4)+1):i*4,1:length);
%     v_temp = ball_vector(((i*4)+1):i*8,1:length);

    %Extract p and v for ball number n
    p_temp = ball_vector((((i-1)*8)+1):4+(i-1)*8,dt);
    v_temp = ball_vector(((i-1)*8)+5:i*8,dt);
    
    %Calculate new values for p and v
    [p_temp,v_temp,count] = propagation(p_temp, v_temp, count, dt);

    %Insert new values at next time frame
    ball_vector((((i-1)*8)+1):4+(i-1)*8,dt+1) = p_temp;
    ball_vector(((i-1)*8)+5:i*8,dt+1) = v_temp;
    
end  
    %Extract newly calculated p and v
    p1_temp = ball_vector(1:3,dt+1);
    p2_temp = ball_vector(9:11,dt+1);
    v1_temp = ball_vector(5:7,dt+1);
    v2_temp = ball_vector(13:15,dt+1);
    
    %Calculate distance between balls
    d = (p1_temp(1)-p2_temp(1))^2 + (p1_temp(2)-p2_temp(2))^2 + (p1_temp(3)-p2_temp(3))^2;
    dist = sqrt(d);
    
    %tangent = [p1_temp(1)-p2_temp(1) , p1_temp(2)-p2_temp(2), p1_temp(3)-p2_temp(3)];
    %normal = 
    %här borde nog något med tangent och normalriktningarna räknas ut, men
    %vi ger bollarna varandras hastighetsvektorer vid kollision istället,
    %om dist <2r kolliderar bollarna:
    
    %If collision occur
    if dist < 2*r
        ball_vector(5:7,dt+1) = v2_temp*0.95;
        ball_vector(13:15,dt+1) = v1_temp*0.95;
        count = count+1;
    end
    
    %Avoid infinite bouncing
    if count>100
        ball_vector(5:7,dt:length+1)=0;
        ball_vector(13:15,dt:length+1)=0;

    end
end

%Plot Bonding box and balls
 figure;
 plotcube([300 300 30],[0 0 0],.3,[1.0 .73 1.0]); 
 hold on
 plot3(ball_vector(1,:), ball_vector(3,:), ball_vector(2,:), 'b'); 
 hold on 
 plot3(ball_vector(9,:), ball_vector(11,:), ball_vector(10,:), 'g'); 
 axis square
 title('Position of Bouncing Ball 3D')
 
 xlabel('Position (x)')
 ylabel('Position (z)')
 zlabel('Position (y)')
 grid on
