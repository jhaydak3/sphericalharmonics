function [x, y, z] = sHarmonic(l,m)
%%[x y z] = orbitalSTLmaker(l,m)
% l,m are quantum numbers
% [x y z] returns the xyz coordinates used for triangulation of the faces

%This program will generate a solid shell of a given orbital. 
% More accurately, it will create a shell of a given spherical harmonic.
% p orbitals can be generated with l = 1, d with l = 2, and so on. Consult
% any physical chemistry book for the right value of l and m necessary.

%nthickness is the normal shell thickness. The value is in mm, and 3mm has
%been very good for stability.
%cthickness is the thickness of the plane at the center used for
%stabilization. Without this, the final print is prone to falling apart.
%These values can be changed and experimented with.

nthickness = 1; 
innerThickness = 1;
clf %clear figure

%input angle for generating points
%These angles can be changed to make a partly hollow object
phiangle = 2.00*pi;
thangle = 1*pi;

global spoints innerThickness
spoints = 50; % This value determines the sampling size

%generate grid used for plotting and calculating points
[phi, theta]= meshgrid(0:(phiangle/spoints):phiangle,0:(thangle/spoints):thangle);

%spherical harmonic calculations
r = sCalc(l,m,phi,theta) ;
x = abs(r) .* sin(theta) .* cos(phi) ;
y = abs(r) .* sin(theta) .* sin (phi);
z = abs(r) .* cos (theta);
%Optional plotting. The code can be uncommented and breaked here to view a
%surface of the spherical harmonics.
hold on
 f = surf(x,y,z);
 set(f,'LineStyle','none')
 xlabel('x')
 ylabel('y')
 zlabel('z')
 %colormap(redgreencmap([2]))
 view(40,30)
axis square


[rows,c] = size(z);
th = zeros(rows,c);
climit = .0001 * range(x(:)); %defines area over which thickness editted
D = climit;

%%OPTIONAL SCALING%%
%Multiplying by 100 makes it more 
x = 100*x; y = 100*y; z = 100*z;
th(:) = nthickness;

%Call to function that makes STL file
[F,V] = surf2solid(x,y,z,'triangulation','f','thickness',th);


sh = surf2solid(x,y,z,'triangulation','f','thickness',th);
%second call to the surf2solid causes it to make a plot of what the solid
%will look like
surf2solid(x,y,z,'triangulation','f','thickness',th);

%write the file 
    s1 = sprintf('%d%d',l,m);
    s = ['sharmonic\Y' s1 'Harmonic.stl'];
    stlwrite(s,sh);
end

function z = sCalc(l,m,phi,theta)
%%This helper function calculates the value of the spherical harmonics. In

    N = 2 *l + 1;
    N = N/ (4*pi); 
    N = N *  factorial(l-abs(m));
    N = N / factorial(l + abs(m));
    N = N^(1/2);
    
    L = legendre(l,cos(theta));
    
    P = N * L(abs(m)+1,:,:);
    P = (-1)^m * P;
    P = reshape(P(:), size(theta));
    %As per one popular convention, only the real part of the function is
    %plotted.
    y2 = exp(1i*m.*phi);

    z = y2.*P;
    %change this bottom line to play with it!
    z = z;

end