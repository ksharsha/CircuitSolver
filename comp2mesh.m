function [ I ] = comp2mesh( I_plot, N_m, N_wi, C_wi, matrix_index, centers, radii )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Generate test values

%clear all;
%%inputs


f = 60;     %frequency
clear N;
N = N_m';
Nf = N;
%N = [5 6;4 95;42 7;95 5;42 52; 94 54; 44 92; 96 94]; % node matrix containg (x,y) cordinate of each node in each row 

%N_wi = [1 0 0 1;1 0 1 0;1 1 0 1; 0 1 0 1;1 0 1 1; 0 1 1 1; 1 1 1 0; 0 1 1 0];
% N_wi(i,j) is 1 if there is wire going from node i in j direction
% j = 1-> right
%     2-> left
%     3-> up
%     4-> down
R = [matrix_index(:,3) matrix_index(:,2) matrix_index(:,5) matrix_index(:,4)];
%R = [41 20 45 22; 93 21 97 24;39 68 45 72; 93 66 97 70; 63 50 68 55 ];
%[R(i,:) = [x1 y1 x2 y2] for resistor i

% while( time ==1)
%     
% disp('Enter value of ');
% disp(size(R,1));
% val_R = input('resistors');
% if(size(val_R,2)== size(R,1))
%     break;
% end
% end
%val_R = [3 12 8 2 10];
%val_R(i) = value of resistor i
% Y_end ,Y_start, X_end, X_start
C = [];
val_C = [];
%V = [6 47 6];
V = [centers(:,1) centers(:,2) radii];
%V(i,:) = [x y radii]
% voltages are positive on upper side and also on the left side ie. [+ -]
% or [+
%     -]


L = [];
val_L = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1-> right
%2-> left
%3-> up
%4-> down
% N -> contains (only orthogonal intersctions)nodes of the form (x,y)
% R,C,L -> contains diagonal corners (x1,y1,x2,y2)
% V -> contains centre of the voltage sources and radius (x,y,r)
% N_wi -> contains 4 columns correspondng to each direction and they will
% be 1 if there is a wire in that direction and will be 0 otherwise (right left up down)
% val_R,val_C,val_L,val_V -> contains values of corresponding R,C,L,V
% outputs
% du_N -> a matrix containing duplicate nodes clubed together,number of
% different nodes = size(du_N,1). The node numbers corresponding to the
% input matrix are given from the 2nd colum of each row. The first column
% of each row gives the number of duplicate nodes in each cluster
% Z_N,V_N -> and n*n matrix where element (i,j) gives the impedence and voltage source connected
% between node i & j
% adj_N -> element (i,j) is 1 if there is a connection between node i and j

% finding reactance of C and L
for i = 1:size(val_C,2)
    val_C(i) = -(1/(2*pi*f*val_C(i)))*(i);
end

for i = 1:size(val_L,2)
    val_L(i) = (2*pi*f*val_L(i))*(i);
end

if(size(R,1)~=0)
Rp = [(R(:,1)+R(:,3))/2 (R(:,2)+R(:,4))/2];
end

%Read values of resistors
val_R = [];
for i = 1: size(R,1)
    t_I = I_plot;
    %R
    t_I(R(i,4):R(i,2),R(i,3):R(i,1),1) = 0;
    imshow(t_I);
    
    disp('Enter value of highlighted resistor');
    
    k= input('resistors');
    val_R = [val_R k];
    
end

%
if(size(C,1)~=0)
Cp = [(C(:,1)+C(:,3))/2 (C(:,2)+C(:,4))/2];
end
if(size(L,1)~=0)
Lp = [(L(:,1)+L(:,3))/2 (L(:,2)+L(:,4))/2];
end

Vp = V(:,1:2);
r_m = mean(V(:,3));

Z_N = zeros(size(N,1));
V_N = zeros(size(N,1));
%making window for nodes 
for i=1:size(N,1)
    w_N(i,:) = [N(i,1)+r_m N(i,2)+r_m N(i,1)-r_m N(i,2)-r_m ];
end

%making window for voltage source 
for i=1:size(V,1)
    w_V(i,:) = [V(i,1)+V(i,3) V(i,2)+V(i,3) V(i,1)-V(i,3) V(i,2)-V(i,3) ];
end

val_V = [];

for i = 1: size(V,1)
    t_I = I_plot;
    %R
    t_I(w_V(i,4):w_V(i,2),w_V(i,3):w_V(i,1),1) = 0;
    imshow(t_I);
    
    disp('Enter value of highlighted resistor');
    
    k= input('resistors');
    val_V = [val_V k];
    
end


clear V;
V = w_V;
%matrix with centre points of everything in the circuit
 p_CO = N;
 if(size(R,1)~=0)
 p_CO(size(p_CO,1)+1:size(p_CO,1)+size(Rp,1),:) = Rp;
 end
 if(size(C,1)~=0)
 p_CO(size(p_CO,1)+1:size(p_CO,1)+size(Cp,1),:) = Cp;
 end
 if(size(L,1)~=0)
 p_CO(size(p_CO,1)+1:size(p_CO,1)+size(Lp,1),:) = Lp;
 end
 p_CO(size(p_CO,1)+1:size(p_CO,1)+size(V,1),:) = Vp;
 
 %matrix with window of everything in the circuit
 w_CO = w_N;
 if(size(R,1)~=0)
 w_CO(size(w_CO,1)+1:size(w_CO,1)+size(R,1),:) = R;
 end
 if(size(C,1)~=0)
 w_CO(size(w_CO,1)+1:size(w_CO,1)+size(C,1),:) = C;
 end
 if(size(L,1)~=0)
 w_CO(size(w_CO,1)+1:size(w_CO,1)+size(L,1),:) = L;
 end
 w_CO(size(w_CO,1)+1:size(w_CO,1)+size(V,1),:) = V;

 %generate matrix with all the values
 val_CO = zeros(size(p_CO,1));
 val_CO(size(N,1)+1:size(N,1)+size(R,1)) = val_R;
 val_CO(size(N,1)+size(R,1)+1:size(N,1)+size(R,1)+size(C,1)) = val_C;
 val_CO(size(N,1)+size(R,1)+size(C,1)+1:size(N,1)+size(R,1)+size(L,1)+size(C,1)) = val_L;
 val_CO(size(N,1)+size(R,1)+size(C,1)+size(L,1)+1:size(N,1)+size(R,1)+size(L,1)+size(C,1)+size(V,1)) = val_V;
 
 %% compute distance b/w each components and nodes
 %d_CO -> distance matrix
 d_CO = 10000*ones(size(p_CO,1),size(p_CO,1));
 o_CO = zeros(size(p_CO,1),size(p_CO,1));
 for i =1:size(p_CO,1)
    for j =i+1:size(p_CO,1)
        
        d1 = p_CO(i,1)-p_CO(j,1);
        d2 = p_CO(i,2)-p_CO(j,2);
        %define distance matrix ->distance between everything in the matrix
        d_CO(i,j) = sqrt(d1.^2 + d2.^2);
        d_CO(j,i) = d_CO(i,j);
        
        %obtaining orientation
        %1-> right
        %2-> left
        %3-> up
        %4-> down
        %o_CO -> orientation matrix
        if((p_CO(i,1)<max(w_CO(j,1),w_CO(j,3)) && (p_CO(i,1)>min(w_CO(j,1),w_CO(j,3)))))
            if(d2>0) %left
                o_CO(i,j) = 3;
                o_CO(j,i) = 4;
            else
                o_CO(i,j) = 4;
                o_CO(j,i) = 3;
            end
        elseif((p_CO(i,2)<max(w_CO(j,2),w_CO(j,4))&& (p_CO(i,2)>min(w_CO(j,4),w_CO(j,2)))))
            if(d1>0) %left
                o_CO(i,j) = 2;
                o_CO(j,i) = 1;
            else
                o_CO(i,j) = 1;
                o_CO(j,i) = 2;
            end
        end
    end
 end
%% 
 %we have distance and orientation between each components
 
 N_wi_d = ones(size(N,1),4);
 N2N = zeros(size(N,1),4);
 CO_d = zeros(size(p_CO,1),1);
 adj_N = zeros(size(N,1),size(N,1));
 for i = 1:size(N,1) %checking each node
    for j = 1:4 %checking each direction
        
        if(N_wi(i,j)==1 && N_wi_d(i,j)==1)
            Z = 0+0i;
            V_s = 0;
            flag_n = 0;
            i1 = i;
            while (flag_n ==0)
                t_min = 100000;
                for k = 1:size(p_CO,1)
                    if(o_CO(i1,k)== j && d_CO(i1,k)<t_min && k~=i1 && CO_d(k)== 0)
                        t_min = d_CO(i1,k);
                        min_ind = k;
                    end
                end
                if(min_ind > size(N,1))
                    CO_d(min_ind) = 1;
                    if(min_ind > size(p_CO,1)-size(V,1))
                        V_s = V_s + val_CO(min_ind);
                    else
                        Z = Z+val_CO(min_ind);
                    end
                    i1 = min_ind;
                else
                    flag_n = 1;
                    adj_N(i,min_ind) = 1;
					adj_N(min_ind,i) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%is it necessary
                    N2N(i,j) = min_ind;
                    if(j ==1 || j == 3)
                        N2N(min_ind,j+1) = i;
                    else
                        N2N(min_ind,j-1) = i;
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(j ==1 || j==4)
                            V_N(i,min_ind) = -V_s;
                        else
                            V_N(i,min_ind) = V_s;
                        end
                        V_N(min_ind,i) = -V_N(i,min_ind);
                        Z_N(i,min_ind) = Z;
                        Z_N(min_ind,i) = Z;
                    N_wi_d(i,j) = 0;
                    N_wi_d(min_ind,j+rem(j,2)-rem(rem(j,2)+1,2)) = 0;
                end
            end
        end
    end
    CO_d(i) = 1;
 end
 
 %% Finding cycles
 C = f_cycle(adj_N,[1],1);
 for i = 2:size(adj_N,1)
    C = [C ;f_cycle(adj_N,[i],i)];
 end
 num = zeros(size(C,1),1);
 for i = 1:size(C,1)
    num(i,1) = -1;
    for j=1:size(C,2)
        if(C(i,j)~=0)
            num(i,1) = num(i,1)+1;
        end
    end
 end
 
cy_N = [num C(:,1:(max(num)+1))];
%remove redundant loops in cy_N
  
 for i = 1: size(cy_N,1)
     for j = i+1:size(cy_N,1)
         if(j<=size(cy_N,1) && i<=size(cy_N,1))
         if(cy_N(i,1)==cy_N(j,1) && i~=j)
             k = 0;
             while (k<=cy_N(i,1)+1)
               if(j<=size(cy_N,1) && i<=size(cy_N,1))  
                if(cy_N(i,1)==cy_N(j,1) && i~=j)
                 if(cy_N(i,2:(cy_N(i,1)+1)) == circshift(cy_N(j,2:(cy_N(j,1)+1))',k)' | cy_N(i,2:(cy_N(i,1)+1)) == fliplr(circshift(cy_N(j,2:(cy_N(j,1)+1))',k)'))
                    cy_N(j,:) = [];
                    %j = j-1;
                    k = -1;
                 end
                end
               end
               k = k+1;
            end
         end
         end
         
     end
 end
 t_cyN = cy_N;
 tmp_cy = cy_N(:,2:end);
 index = [];
 
 %mark the shortest loops for each 'edge'
 
 for i = 1:size(adj_N,1)-1
     for j = i+1:size(adj_N,1)
         if(adj_N(i,j)~=0)
            l = size(adj_N,1)+2;
            ind =0;
            for k = 1:size(cy_N,1)
                i1 = find(tmp_cy(k,1:cy_N(k,1))==i);
                if(size(i1,2)~=0)
                    if(cy_N(k,i1+2)==j || (cy_N(k,i1)==j && i1~=1) || (cy_N(k,cy_N(k,1)+1)==j && i1==1))
                        if(cy_N(k,1)<l)
                            l = cy_N(k,1);
                            ind = k;
                        end
                    end
                end
                
            end
            t=find(index ==ind);
                if(size(t,2)==0 & ind~=0)
                    if(ind == 6)
                        check = [i,j];
                    end
                    index = [index ind];
                end
         end
     end
 end
 % remove the rest of the loops
 new_cy = [];
 for i = 1:size(index,2)
    new_cy = [new_cy; cy_N(index(i),:)];
 end
 
 clear cy_N;
 cy_N = new_cy;



%% Cycle2Mesh

%cy_N -> matrix with each row containing one cycle eg: cy_N(1,:) = [3 1 3 5
%1] 1st element is the number of nodes in the cycle
%i_cy ->impedence matrix
%cy_N = [5 1 3 5 7 2 1;4 3 4 6 5  3 0;4 5 6 8 7 5 0];
%%
i_cy = zeros(size(cy_N,1),size(cy_N,1));
v_cy = zeros(1,size(cy_N,1));

for i = 1:size(cy_N,1)
    for j = 2:cy_N(i,1)+ 1
        v_cy(i) = v_cy(i) + V_N(cy_N(i,j),cy_N(i,j+1));
        i_cy(i,i) = i_cy(i,i)+Z_N(cy_N(i,j),cy_N(i,j+1)); % generates all (i,i) terms
        for k = i+1:size(cy_N,1)
            %flag = 0;
            
                i1 = find(cy_N(k,2:cy_N(k,1)+1)==cy_N(i,j));
                if(size(i1,2)~=0)
                    %tmp = cy_N(k,2:end);
                    if(cy_N(k,i1+2) == cy_N(i,j+1) || (cy_N(k,i1) == cy_N(i,j+1) && i1~=1) || (i1==1 && cy_N(k,cy_N(k,1)+1) == cy_N(i,j+1)))
                        i_cy(i,k) = i_cy(i,k)-Z_N(cy_N(i,j),cy_N(i,j+1));
                        i_cy(k,i) = i_cy(k,i)-Z_N(cy_N(i,j),cy_N(i,j+1));
                    end
                end
            
        end
    end
end
%cycle current
%%
cy_I = inv(i_cy)*v_cy';

I = zeros(size(Z_N));

for i = 1:size(cy_N,1)
    for j = 2:cy_N(i,1)+1
        I(cy_N(i,j),cy_N(i,j+1)) = I(cy_N(i,j),cy_N(i,j+1))+cy_I(i);
        I(cy_N(i,j+1),cy_N(i,j)) = I(cy_N(i,j+1),cy_N(i,j))-cy_I(i);
    end
end

Im = I_plot(:,:,1);
In = imread('img.jpg');
In = In(1:220,1:160,1);
for i =1:size(Nf,1)
    %if(i<=10)
    i  
    %if(rem(i,10)~=0)
        Im(Nf(i,2)-10:Nf(i,2)+11,Nf(i,1)-7:Nf(i,1)+8) = In(floor(((i-1)/10))*22+1:floor(((i-1)/10+1))*22,rem((i-1),10)*16+1:(rem((i-1),10)+1)*16);
        imshow(Im);
    %end 
    %end
end
imshow(Im);

Nf

end