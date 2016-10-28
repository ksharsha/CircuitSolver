function draw_bound(I,YS,YE,XS,XE) % I is actual image grayscale
I = rgb2gray(I);
k = size(I)
R = zeros(k(1),k(2),3);
R(:,:,1) = I;
R(:,:,2) = I;
R(:,:,3) = I;

l = length(XS);
for i = 1:l
    if XS(i) ~=0
    for m = (XS(i)-5):(XE(i)+5)
        R(m,YS(i)-5,1) = 255;
        R(m,YE(i)+5,1) = 255;
        R(m,YS(i)-5,2) = 0;
        R(m,YE(i)+5,2) = 0;
        R(m,YS(i)-5,2) = 0;
        R(m,YE(i)+5,2) = 0;
    end
    for n = (YS(i)-5):(YE(i) + 5)
        R(XS(i)-5,n,1) = 255;
        R(XE(i)+5,n,1) = 255;
        R(XS(i)-5,n,1) = 0;
        R(XE(i)+5,n,1) = 0;
        R(XS(i)-5,n,1) = 0;
        R(XE(i)+5,n,1) = 0;
    end
    end
end
figure,imshow(R);


