function [start_Y_cap,end_Y_cap,X_start_cap,X_end_cap] =  region_I(IM_H,Q);
% takes greyscale image
if Q == 1
    J = im2bw(IM_H);
else
    J = im2bw(IM_H');
end
kk = size(J);   
% Apply connected components
k = 1;
R =zeros(kk);
    for i = 2:kk(1)
        for j = 2:kk(2)
            if J(i,j) == 1
                if J(i,j-1) == 0 && J(i-1,j) == 0
                    R(i,j) = k;
                    k = k + 1;
                elseif J(i,j-1) == 0 && J(i-1,j) == 1
                    R(i,j) = R(i-1,j);
                elseif J(i,j-1) == 1 && J(i-1,j) == 0
                    R(i,j) = R(i,j-1);
                elseif J(i,j-1) == 1 && J(i-1,j) == 1
                    R(i,j) = R(i-1,j);
                    if R(i,j-1) ~= R(i-1,j)
                        for l = 1:i
                            for m = 1:j
                                if R(l,m) == R(i,j-1)
                                    R(l,m) = R(i-1,j);%tracing back
                                end
                            end
                        end
                    end
                end
            end
        end
    end
reg = max(max(R));
reg
if reg == 0
    start_Y_cap =[0];
    end_Y_cap = [0];
    X_start_cap = [0];
    X_end_cap = [0];
else
    cap_num = reg/2;
    start_Y = zeros(1,reg); % store starting Y points of regions
    end_Y = zeros(1,reg); % store starting Y points of regions
    X = zeros(1,reg); % store X points of regions
    for h = 1:reg
        for i = 1:kk(1)
            for j = 1:kk(2)
                if R(i,j) == h
                    X(h) = i;
                    start_Y(h) = j;
                    while R(i,j)==h
                        j = j+1;
                    end
                    end_Y(h) = j-1;    
                    break;
                end
            end
        end
    end
% XX = unique(X);
% X_start = XX(1:2:length(XX));
% X_end = XX(2:2:length(XX));
% l = 1;
% g = length(X_start);
% for i = 1:g
%     ind = find(X == X_start(i));
%     % length of ind denotes how many cap in a row
%     hh = length(ind);
%     for d = 1:hh
%         sY(l) = start_Y(ind(d));
%         eY(l) = end_Y(ind(d));
%         Xe(l) = X_end(i);
%         Xs(l) = X_start(i);
%         l = l+1;
%     end
% end
    l = 1;
    YY = unique(start_Y);
    for i = 1:length(YY)
        ind = find(start_Y == YY(i));
        N = length(ind);
        if N ~= 0
            X_start = X(ind(1:2:N));
            X_end = X(ind(2:2:N));
            for d = 1:(N/2)
                sY(l) = YY(i);
                eY(l) = end_Y(ind(d));
                Xs(l) = X_start(d);
                Xe(l) = X_end(d);
                l = l+1;
            end
        end
    end
    if Q == 1
        start_Y_cap = sY;
        end_Y_cap = eY;
        X_start_cap = Xs;
        X_end_cap = Xe;
    else
        start_Y_cap = Xs;
        end_Y_cap = Xe;
        X_start_cap = sY;
        X_end_cap = eY;
    end
end


    




            



