
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