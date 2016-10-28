% only for horizontal .... For vertical take transpose of image and pass !
% think about resizing the images before doing process
function[C_wi,JJ] = hough_cap(I,acc_threshold,cap_range);

% Pass non Binary image
% I = imread('circuit4.jpg');
k = size(I);
if length(k) == 3
    I = rgb2gray(I);
end
%imshow(I);
for im = 1:2
    if im == 1
        J = not(im2bw(I));
    else
        J = not(im2bw(I'));
    end
    %J = bwmorph(JJ,'thin'); % image converted to binary
    %figure,imshow(J);
    k = size(J);
    A = zeros(1,k(1)); % along rows ,, for horizontal lines
    for m = 1:k(1)
        for n = 1:k(2)   % horizontal => theta = 0 =>  m = p
            if J(m,n) == 1
                for p = 1 : k(1)
                    if abs(m - p) == 0
                        A(p) = A(p) + 1;
                    end
                end
            end
        end
    end
    %figure,plot(A);
    B = zeros(k(1),k(2));
    %ind = find((A >= max(A)/10) & (A <= max(A)/1.5))
    ind = find((A >= max(A)/acc_threshold)); % retrieving only that A's which are higher than a particular threshold
    C = zeros(size(A));
    C(ind) = 1;
    %figure,plot(C);
    %figure,imshow(C);
    count = 0;
    e = zeros(1,k(1));
    s = zeros(1,k(1));
    f = 1;
    o = 1;
    for h = (cap_range+1):(k(1)- cap_range) % applying a vertical window of 21 pixel length on A to check if parallel lines appear
        s = e;
        o = f;
        for g = -cap_range:cap_range
            if C(h+g) == 1
%               if g == -10
%                  count = count + 1;
%                  e(f) = h + g;
                if (h+g-1) ~= 0      %f = f+1;
                    if C(h+g-1) ~= 1 
                        count = count + 1;
                        e(f) = h+g;
                        f = f + 1;
                    end
                end
            end
        end
        if count <= 1
            e = s;
            f = o;
        end
    count = 0;
    end
% e will store required indices with zeros padded ahead, as length was
% getting more than required while saving
% here we have got the required 
    for o = 1:length(e)
        if e(o) == 0;
            h = o - 1;
            break;
        end
    end
    ind2 = e(1:h);
    ind2 = unique(ind2);
% w contains indices !!
% apply connected components to remove minor blobs
    v = size(ind2);
    if v(2) ~= 0
        for d = 1:v(2)
            for m = 1:k(1)
                for n=1:k(2)
                    if J(m,n) == 1
                        if abs(m - ind2(d)) == 0
                            B(m,n) = 1;
                        end
                    end
                end
            end
        end
    end

%figure,imshow(B) % B contains indices
% apply connected components to remove minor blob
% structuring element for DILATE
    if im == 1
        Nhood = [1 1 1];
        BB2= imdilate(B,Nhood);
%         figure,imshow(BB2);
        BB = bwmorph(BB2,'thin');
    else
        Nhood = [1;1;1];
        HH2 = imdilate(B.',Nhood);
%         figure,imshow(HH2);
        HH = bwmorph(HH2,'thin');
    end
end
% we have got images as BB and HH for horizontal and vertical
%figure,imshow(BB);

high_w = [0 0];
% our BB contains horizontal and HH contains vertical info
% till here it's correct
for im =1:2
    if im == 1 
        J = BB;
    else
        J = HH;
    end
    k = 1;
    kk = size(J);
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
    % R contains region now
    %figure,imshow(R);
    m = max(max(R)); % how many regions we got
    if m ~= 0
        weight_R = zeros(1,m);
        % got the no. regions till here 
        % now we need the size of regions
        for i = 1:m
            ind = find(R == i);
            weight_R(i) = length(ind);
        end
        %if m~= 0
        high_w(im) = max(weight_R);
    else
        high_w(im) = 0;
        weight_R = 0;
        
    end
        %else
        %    high_w(im) = 0;
        %end
    if im == 1
        RH = R;
        weight_RH = weight_R;
    else
        RV = R;
        weight_RV = weight_R;
    end
end

        
        
        

% we have RH and RV as separate storing regions 
high_w2 = max(high_w);
for im = 1:2
    if im == 1
        weight_R = weight_RH;
        R = RH;
    else
        weight_R = weight_RV;
        R = RV;
    end
    ind = find(weight_R < (high_w2/3));
    s = size(R);
    for i = 1:s(1)
        for j = 1:s(2)
            if R(i,j) ~= 0
                for h = 1:length(ind)
                    if R(i,j) == ind(h);
                        R(i,j) = 0;
                    end
                end
            end
        end
    end
    if im == 1
        IM_H = R;
    elseif im == 2
        IM_V = R;
    end
%     figure,imshow(R);
end

        
    

[a,b,c,d] = region_I(IM_H,1);
[aa,bb,cc,dd] = region_I(IM_V,2);

Y_start = horzcat(a,aa);
Y_end = horzcat(b,bb);
X_start = horzcat(c,cc);
X_end = horzcat(d,dd);

Y_start = Y_start(:,1:(end-1))
Y_end = Y_end(:,1:(end-1))
X_start = X_start(:,1:(end-1))
X_end = X_end(:,1:(end-1))

C_wi = [Y_end' X_end' Y_start' X_start'];
[m,n] = size(X_start);
figure, imshow(I), hold on;
% for i = 1 : n
%     for j = 1 : m
%     h(i) = imrect (gca, [X_start(j,i) Y_start(j,i) X_end(j,i) Y_end(j,i)]);
%     end
% end

% removing cap Part
J = not(im2bw(I));
f = length(Y_start);


for i = 1:f
    if X_start(i) ~= 0
        for x = (X_start(i)-2) : (X_end(i) + 2)
            for y = (Y_start(i)-2) : (Y_end(i)+2)
                J(x,y) = 0;
            end
        end
    end
end
JJ = J;





