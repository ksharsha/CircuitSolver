function draw_bound_node(I,YS,YE,XS,XE) % I is actual image grayscale
I = rgb2gray(I);
k = size(I)
R = zeros(k(1),k(2),3);
R(:,:,1) = I;
R(:,:,2) = I;
R(:,:,3) = I;

L = length(XS);
for l = 1:L
    if YS(l) == 0
        YS(l) = 1;
    end
    if YE(l) == 0
        YE(l) = 1;
    end
    if XS(l) == 0
        XS(l) = 1;
    end
    if XE(l) == 0
        XE(l) = 1;
    end
end
        
for i = 1:l 
    for m = XS(i):XE(i)
        R(m,YS(i),1) = 255;
        R(m,YE(i),1) = 255;
        R(m,YS(i),2) = 0;
        R(m,YE(i),2) = 0;
        R(m,YS(i),2) = 0;
        R(m,YE(i),2) = 0;
    end
    for n = YS(i):YE(i)
        R(XS(i),n,1) = 255;
        R(XE(i),n,1) = 255;
        R(XS(i),n,2) = 0;
        R(XE(i),n,2) = 0;
        R(XS(i),n,3) = 0;
        R(XE(i),n,3) = 0;
    end
   
end
figure,imshow(R);
