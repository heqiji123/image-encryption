function [a,b,c,d] = IWT(B)
im=double(B);
lshaarInt = liftwave('haar','int2int');
els = {'p',[-0.125 0.125],0};
lsnewInt = addlift(lshaarInt,els);
[a,b,c,d] = lwt2(im,lsnewInt);
% a=uint8(a);
% b=uint8(b);
% c=uint8(c);
% d=uint8(d);
end