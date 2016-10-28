I = imread('rlcinput.png');
%figure('name','Image'),imshow(I);
%figure('name','Histogram'),imhist(I);
level = graythresh(I);
J = im2bw(I,level);
%figure('name','Otsus segmentation'),imshow(J);
dimJ = size(J);
R = zeros(dimJ);
k = 1;
scount = 1;
J(1,:) = 1;
J(:,1) = 1;
for i = 2:dimJ(1)
    for j = 2:dimJ(2)
        if J(i,j) == 0
            if J(i,j-1) == 1 && J(i-1,j) == 1
                R(i,j) = k;
                k = k + 1;
            elseif J(i,j-1) == 1 && J(i-1,j) == 0
                R(i,j) = R(i-1,j);
            elseif J(i,j-1) == 0 && J(i-1,j) == 1
                R(i,j) = R(i,j-1);
            elseif J(i,j-1) == 0 && J(i-1,j) == 0
                R(i,j) = R(i-1,j);
                if R(i,j-1) ~= R(i-1,j)
                    for l = 1:dimJ(1)
                        for m = 1:dimJ(2)
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
figure('name','segmented'),imshow(J);
sign = [1,1,0];
flag = 0;
for i = 2:dimJ(1)
    for j = 2:dimJ(1)
        if R(i,j) ~= 0
            dimS = size(sign);
            for k = 1:dimS(1)%checks which region the current pixel belongs to and increments that particular value
                if R(i,j) == sign(k,2)
                    sign(k,3) = sign(k,3) + 1;%gives the total number of pixels in that region
                    flag = 1;
                    break
                end
            end
            if flag == 0 %new region
                sign = [sign;dimS(1)+1,R(i,j),1];%adding a new roww to the matrix
            end
            flag = 0;
        end
    end
end
%identifying the region with maximum pixels%
dimS = size(sign);
mx = max(sign(:,3));

for i = 1:dimS(1)
    if sign(i,3) == mx
        obj = sign(i,2);
        break
    end
end
%%modifying the program
%%making the mx go to zero
for i = 1:dimS(1)
    if sign(i,3) == mx
        sign(i,3)=1;
        break
    end
end
mx1 = max(sign(:,3));%%second maxima
for i = 1:dimS(1)
    if sign(i,3) == mx1
        obj1 = sign(i,2);
        break
    end
end
if mx1<(mx/2)
    obj1=obj;
end
for i = 1:dimS(1)
    if sign(i,3) == mx1
        sign(i,3)=1;
        break
    end
end
mx2 = max(sign(:,3));%%second maxima
S = zeros(dimJ);
for i = 1:dimJ(1)
    for j = 1:dimJ(2)
        if R(i,j) == obj||R(i,j)== obj1
            S(i,j) = 0;
        else
            S(i,j) = 1;
        end
    end
end
figure('name','Object Plot'),imshow(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%done obtaining only the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%circuit%%%%%%%%%%%%%%%%%%%%%%%%%%%%


grad=zeros(dimJ);
for i = 2:dimJ(1)-1
    for j = 2:dimJ(2)-1
       grad(i,j)=S(i+1,j+1)-S(i,j);%%x gradient
    end
end
% figure('name','Object Plot'),imagesc(grad);
grad1=zeros(dimJ);
for i = 2:dimJ(1)-1
    for j = 2:dimJ(2)-1
       grad1(i,j)=grad(i+1,j+1)-grad(i,j);%%x gradient
    end
end
% figure('name','Object Plot'),imagesc(grad1);

%%try to fit a polynomial



























