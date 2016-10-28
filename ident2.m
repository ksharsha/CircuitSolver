function [ I ] = ident2( I_plot,matrix_index, centers, radii, a )
[C_wi,e] = hough_cap(a,10,10);
[start_X,end_X,start_Y,end_Y] = getNodes(e,15);


BW = edge(a,'sobel');
%figure,  imshow(imcomplement(BW));
BW = imcomplement(BW);
N = [start_X; end_X; start_Y; end_Y];
w_x = abs(start_X - end_X);
w_y = abs(start_Y - end_Y);
[m,n] = size(N); % m=4 always, n=no. of nodes
N_m = [(start_X+end_X)/2; (start_Y+end_Y)/2];
N_wi = zeros(4,n);
r_w = [start_X; end_X;end_Y; end_Y+w_y];
l_w = [start_X; end_X;start_Y-w_y; start_Y];        
d_w = [end_X; end_X+w_x ;start_Y; end_Y];
u_w = [start_X-w_x; start_X; start_Y; end_Y];




for i = 1:n
    c_r = 0; % right
    for j = max(1,r_w(1,i)) : min(size(BW,1),r_w(2,i))
        for k = max(1,r_w(3,i)) :  min(size(BW,2),r_w(4,i)) 
            if (BW(j,k)== 0)
                c_r = 1;
            end 
        end
    end    
    if (c_r == 1)
        N_wi(1,i) = 1;
    end   
    
    c_l = 0; % left
    for j = max(1,l_w(1,i)) :  min(size(BW,1),l_w(2,i))
        for k = max(1,l_w(3,i)) :  min(size(BW,2),l_w(4,i)) 
            if (BW(j,k)== 0)
                c_l = 1;
            end 
        end
    end    
    if (c_l == 1)
        N_wi(2,i) = 1;
    end   
    
    c_u = 0; % up
    for j = max(1,u_w(1,i)) :  min(size(BW,1),u_w(2,i))
        for k = max(1,u_w(3,i)) :  min(size(BW,2),u_w(4,i)) 
            if (BW(j,k)== 0)
                c_u = 1;
            end 
        end
    end    
    if (c_u == 1)
        N_wi(3,i) = 1;
    end   
    
    c_d = 0; % down
    for j = max(1,d_w(1,i)) :  min(size(BW,1),d_w(2,i))
        for k = max(1,d_w(3,i)) :  min(size(BW,2),d_w(4,i)) 
            if (BW(j,k)== 0)
                c_d = 1;
            end 
        end
    end    
    if (c_d == 1)
        N_wi(4,i) = 1;
    end   
end

%T = [N_wi(4,:);N_wi(3,:);N_wi(2,:);N_wi(1,:)];
N_wi = N_wi';

N_m = [N_m(2,:); N_m(1,:)]; 

I = comp2mesh(I_plot, N_m, N_wi, C_wi, matrix_index, centers, radii);

end

