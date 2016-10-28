function[y] = diff_x(Iu)
% I in Binary
I = im2bw(Iu);
k = size(I);
for m = 1:k(1)-1
    for n=1:k(2)-1
        J(m,n) = abs(I(m,n) - I(m,n+1));
        L(m,n) = abs(I(m,n)- I(m+1,n));
    end
end
H = zeros(k(1)-1,k(2)-1);
se = strel('line',5,20);
Jd = imdilate(J,se);
Ld = imdilate(L,se);
for m = 1:k(1)-1
    for n=1:k(2)-1
        if Jd(m,n) == 1
            if Ld(m,n) == 1
                H(m,n) = 1;
            end
        end
    end
end
se = strel('line',40,25);
Hd = imdilate(H,se);
figure,imshow(Jd)
figure,imshow(Ld)
figure,imshow(Hd)
