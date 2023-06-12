function [ P ] = IIWT(a,b,c,d)
lshaarInt = liftwave('haar','int2int');
els = {'p',[-0.125 0.125],0};
lsnewInt = addlift(lshaarInt,els);
a=double(a);b=double(b);c=double(c);d=double(d);
    p= (ilwt2(a,b,c,d,lsnewInt));
    P=double(p);
end
