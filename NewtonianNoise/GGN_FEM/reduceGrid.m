function r_red = reduceGrid(r,model,fac,type)

midi = (model.ng+1)/2;

if strcmp(type,'uniform')
    ci1 = [midi(1):-fac:1 (midi(1)+fac):fac:model.ng(1)];
    ci2 = [midi(2):-fac:1 (midi(2)+fac):fac:model.ng(2)];
    ci3 = [midi(3):-fac:1 (midi(3)+fac):fac:model.ng(3)];
    si1 = find(ci1<=fac);
    si2 = find(ci2<=fac);
    si3 = find(ci3<=fac);
    ci1(1:si1) = ci1(si1:-1:1);
    ci2(1:si2) = ci2(si2:-1:1);
    ci3(1:si3) = ci3(si3:-1:1);
elseif strcmp(type,'cubic')
    i1 = 1:fac:(midi(1)-1);
    i2 = 1:fac:(midi(2)-1);
    i3 = 1:fac:(midi(3)-1);
    i1_c = i1.^3;
    i2_c = i2.^3;
    i3_c = i3.^3;
    i1_c = unique(round(i1_c/i1_c(end)*i1(end)));
    i2_c = unique(round(i2_c/i2_c(end)*i2(end)));
    i3_c = unique(round(i3_c/i3_c(end)*i3(end)));
    i1_c = i1_c(i1_c>0);
    i2_c = i2_c(i2_c>0);
    i3_c = i3_c(i3_c>0);
    
    ci1 = [(midi(1)-i1_c) midi(1) (midi(1)+i1_c)];
    ci2 = [(midi(2)-i2_c) midi(2) (midi(2)+i2_c)];
    ci3 = [(midi(3)-i3_c) midi(3) (midi(3)+i3_c)];
    
    si1 = length(i1_c);
    si2 = length(i2_c);
    si3 = length(i3_c);
    
    ci1(1:si1) = ci1(si1:-1:1);
    ci2(1:si2) = ci2(si2:-1:1);
    ci3(1:si3) = ci3(si3:-1:1);
elseif strcmp(type,'quadratic')
    i1 = 1:fac:(midi(1)-1);
    i2 = 1:fac:(midi(2)-1);
    i3 = 1:fac:(midi(3)-1);
    i1_c = i1.^2;
    i2_c = i2.^2;
    i3_c = i3.^2;
    i1_c = unique(round(i1_c/i1_c(end)*i1(end)));
    i2_c = unique(round(i2_c/i2_c(end)*i2(end)));
    i3_c = unique(round(i3_c/i3_c(end)*i3(end)));
    i1_c = i1_c(i1_c>0);
    i2_c = i2_c(i2_c>0);
    i3_c = i3_c(i3_c>0);
    
    ci1 = [(midi(1)-i1_c) midi(1) (midi(1)+i1_c)];
    ci2 = [(midi(2)-i2_c) midi(2) (midi(2)+i2_c)];
    ci3 = [(midi(3)-i3_c) midi(3) (midi(3)+i3_c)];
    
    si1 = length(i1_c);
    si2 = length(i2_c);
    si3 = length(i3_c);
    
    ci1(1:si1) = ci1(si1:-1:1);
    ci2(1:si2) = ci2(si2:-1:1);
    ci3(1:si3) = ci3(si3:-1:1);
end
    
r_red = r(ci1,ci2,ci3,:);