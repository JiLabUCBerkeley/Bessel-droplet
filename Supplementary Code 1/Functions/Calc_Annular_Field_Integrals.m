function [PSF,I0, I1, I2] = Calc_Annular_Field_Integrals(x, y, z,mask_obj,wavelength,r_angle,refind)
%CALC_ANNULAR_FIELD_INTEGRALS  Calculates the 3 Richards & Wolf integrals used to calc. EM field
%   CALC_ANNULAR_FIELD_INTEGRALS calculates the three integrals I0, I1,
%      and I2 used by Richards and Wolf in the determination of all six
%      components of the EM field in a vicinity of a focus, for a single
%      condition of annular illumination
%input arguments x, y, and z are coordinates near focal point; mask_obj is
%the eletrical field at the back focal plane of objective, innerD and
%outerD define inner diameter and outer diameter of annular ring. 
% 
xlen = length(x);
ylen = length(y);
zlen = length(z);
%
%create cylindrical coord. dimensionless arrays:
xmat = (x.') * ones(1, ylen);   %xlen by ylen matrix, each column filled with x vector
ymat = ones(xlen, 1) * y;       %xlen by ylen matrix, each row filled with y vector
rmat = sqrt((xmat .^ 2) + (ymat .^ 2)) + eps;  %projection of vector from focus to field pt onto xy plane
k=2*pi/wavelength;
uvec = k .* i .* refind .* z;   %uvec = ikz, k = mag. of wavevector in air
vmat = k .* refind .* rmat;  %vmat = k * rperp


%%
%
%define matrices for the definite integrals I0, I1, and I2:
I0 = ones(xlen, ylen, zlen);
I1 = ones(xlen, ylen, zlen);
I2 = ones(xlen, ylen, zlen);
phi=ones(xlen, ylen, zlen);
%
%define functions for integrands of I0, I1, and I2:
I0integrand = inline('sqrt(cos(th)) .* sin(th) .* (1 + cos(th)) .* besselj(0, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
I1integrand = inline('sqrt(cos(th)) .* (sin(th) .^ 2) .* besselj(1, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
I2integrand = inline('sqrt(cos(th)) .* sin(th) .* (1 - cos(th)) .* besselj(2, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
%
interr = 1e-3;  %fractional error permitted for integral evaluation using quad
nullv = [];     %null vector needed for trace argument in quad fcn
%
z_save=z;
for z = 1:zlen  %loop thru all z positions
%     disp(['z=',num2str(z_save(z))])
%    z
    for x = 1:xlen  %loop thru all x positions
%        x
        for y = 1:ylen  %loop thru all y positions
            a=I0integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;
            b=I1integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;
            c=I2integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;           
            I0_tmp=trapz(a);
            I1_tmp=trapz(b);
            I2_tmp=trapz(c);
            I0(x,y,z)=I0_tmp;
            I1(x,y,z)=I1_tmp;
            I2(x,y,z)=I2_tmp;    
            phi(x,y,z)=angle(x+1i*y);
        end
    end
end
% phi=0;
Ex=(I0+I2.*cos(2*phi));
Ey=(I2.*cos(2*phi));
Ez=(-2i*I1.*cos(phi));
PSF=abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
PSF=PSF.^2;
end