function [start_X,end_X,start_Y,end_Y] = getNodes(R,N)
% R is binary 
K = size(R);
Y = zeros(size(R));
% cross of size (2N + 1) X (2N + 1)
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);
d = zeros(1,N);
for i = 1:K(1)
    for j = 2:K(2)
        if R(i,j) == 1
            for k = 1:N
                if (i + k) <= K(1) 
                    a(k) = R(i+k,j);
                end
                if (i-k) >= 1
                    b(k) = R(i-k,j);
                end
                if (j+k)<= K(2)
                    c(k) = R(i,j+k); 
                end
                if (j-k)>= 1
                    d(k) = R(i,j-k);
                end
            end
            a_min = min(a);
            b_min = min(b);
            c_min = min(c);
            d_min = min(d);
            if (a_min == 1) || (b_min == 1)
               if (c_min == 1) || (d_min == 1)
                   % i,j is a node
                   Y(i,j) = 1;
               else
                   Y(i,j) = 0;
               end
            else
                Y(i,j) = 0;
            end
        end
    end
end
figure,imshow(Y)
J =Y;
k = 1;
RR =zeros(size(R));
    for i = 2:K(1)
        for j = 2:K(2)
            if J(i,j) == 1
                if J(i,j-1) == 0 && J(i-1,j) == 0
                    RR(i,j) = k;
                    k = k + 1;
                elseif J(i,j-1) == 0 && J(i-1,j) == 1
                    RR(i,j) = RR(i-1,j);
                elseif J(i,j-1) == 1 && J(i-1,j) == 0
                    RR(i,j) = RR(i,j-1);
                elseif J(i,j-1) == 1 && J(i-1,j) == 1
                    RR(i,j) = RR(i-1,j);
                    if RR(i,j-1) ~= RR(i-1,j)
                        for l = 1:i
                            for m = 1:j
                                if RR(l,m) == RR(i,j-1)
                                    RR(l,m) = RR(i-1,j);%tracing back
                                end
                            end
                        end
                    end
                end
            end
        end
    end
nodes = max(max(RR));
j = 1;
for i = 1:nodes
    [indx,indy] = find(RR == i);
    if indx ~= 0
        start_X(j) = indx(1) - 10;
        if start_X(j) <= 0
            start_X(j) = 1;
        end
        end_X(j) = indx(1) + 10;
        if end_X(j) > K(1)
            end_X(j) = K(1);
        end
        start_Y(j) = indy(1) -10;
        if start_Y(j) <= 0
            start_Y(j) = 1;
        end
        end_Y(j) = indy(1) + 10;
        if end_Y(j) > K(2)
            end_Y(j) = K(2);
        end
        j = j+1;
    end
end

start_X
end_X
start_Y
end_Y

    
   
    
        