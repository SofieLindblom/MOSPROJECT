function [ p_out, v_out, count ] = propagation( p_in, v_in, count, dt )

%Gravity free fall
grav = [0,-9.82,0];
%Step size
delta_t = 0.01;
%p_out = p_in;
%v_out = v_in;

 %If ball hits ground y-direction
    if p_in(2) <= 0 %&& v_in(2) <= 0 
        %x-pos
        v_out(1) = v_in(1);
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        %y-pos
        p_in(2) = 0;
        v_out(2) = -(0.90*v_in(2))+ grav(2)*delta_t;
        p_out(2) = v_out(2)*delta_t;
        %z-pos
        v_out(3) = v_in(3); 
        p_out(3) = p_in(3) + v_out(3)*delta_t;
        count = count +1;
    %If ball hits wall x = 0 in x-direction
    elseif p_in(1) <= 0
        %x-pos
        p_in(1) = 0;
        v_out(1) = -(0.90*v_in(1));
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        % y-pos
        v_out(2) = v_in(2) + grav(2)*delta_t; 
        p_out(2) = p_in(2) + v_out(2)*delta_t;
        % z-pos
        v_out(3) = v_in(3); 
        p_out(3) = p_in(3) + v_out(3)*delta_t;
        count = count +1;
        %If ball hits wall x = 300 in x-direction
    elseif p_in(1) >= 300 
        %x-pos
        p_in(1) = 300;
        v_out(1) = -(0.90*v_in(1));
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        % y-pos
        v_out(2) = v_in(2) + grav(2)*delta_t; 
        p_out(2) = p_in(2) + v_out(2)*delta_t;
        % z-pos
        v_out(3) = v_in(3); 
        p_out(3) = p_in(3) + v_out(3)*delta_t;
        count = count +1;
    %If ball hits wall z = 0 in z-direction
    elseif p_in(3) <= 0
        %x-pos
        v_out(1) = v_in(1); 
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        % y-pos
        v_out(2) = v_in(2) + grav(2)*delta_t; 
        p_out(2) = p_in(2) + v_out(2)*delta_t;
        % z-pos
        p_in(3) = 0;
        v_out(3) = -(0.90*v_in(3)); 
        p_out(3) = p_in(3) + v_out(3)*delta_t;
        count = count +1;
         %If ball hits wall z = 300 in z-direction
    elseif p_in(3) >= 300 
        %x-pos
        v_out(1) = v_in(1); 
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        % y-pos
        v_out(2) = v_in(2) + grav(2)*delta_t; 
        p_out(2) = p_in(2) + v_out(2)*delta_t;
        % z-pos
        p_in(3) = 300;
        v_out(3) = -(0.90*v_in(3)) 
        p_out(3) = p_in(3) + v_out(3)*delta_t;
        count = count +1;
    %If no collision between ball and wall
    else 
        %x-pos
        v_out(1) = v_in(1); 
        p_out(1) = p_in(1) + v_out(1)*delta_t;
        % y-pos
        v_out(2) = v_in(2) + grav(2)*delta_t;
        p_out(2) = p_in(2) + v_out(2)*delta_t;
        % z-pos
        v_out(3) = v_in(3); 
        p_out(3) = p_in(3) + v_out(3)*delta_t;

    end
    
    p_out(4) = dt;
    v_out(4) = dt;
    
%     if count>50 - utanför
%         v_matrix(1:3,dt:2001)=0;
%         p_matrix(1:3,dt:2001)=0;
% 
%     end

end

