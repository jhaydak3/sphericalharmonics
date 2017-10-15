function [x y z] = orbitalSTLmaker(l,m)
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
aux = 0; %% true/false for support harmomic
clf %clear figure

%input angle for generating points
%These angles can be changed to make a partly hollow object
phiangle = 2*pi;
thangle = 1*pi;

global spoints innerThickness
spoints = 80; % This value determines the sampling size

%generate grid used for plotting and calculating points
[phi, theta]= meshgrid(0:(phiangle/spoints):phiangle,0:(thangle/spoints):thangle);

%spherical harmonic calculations
r = sCalc(l,m,phi,theta);
x = abs(r) .* sin(theta) .* cos(phi) ;
y = abs(r) .* sin(theta) .* sin (phi);
z = abs(r) .* cos (theta);
% x = [x x(:,1)];
% x = [x; x(1,:)];
% y = [y y(:,1)];
% y = [y; y(1,:)];
% z = [z z(:,1)];
% z = [z; z(1,:)];
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

%Helper function that makes center stabalizing part(s)
if aux
    sh2 = makeAux(l,m);
end

sh = surf2solid(x,y,z,'triangulation','f','thickness',th);
%second call to the surf2solid causes it to make a plot of what the solid
%will look like
surf2solid(x,y,z,'triangulation','f','thickness',th);

%write the file
s1 = sprintf('%d%d',l,m);
s = ['orbitals\Y' s1 'orbital.stl'];
stlwrite(s,sh);
if aux
    for ndx = 1:length(sh2)
        s = ['Y' s1 'orbitalaux' num2str(ndx) '.stl'];
        stlwrite(s,sh2{ndx});
    end
end
end
function sh2 = makeAux(l,m)
global innerThickness
sh2 = {};
if m == 0
    supportAngles = 0;
else
   supportAngles = 0:pi/m:(pi/m)*(m-1)
end
for rotAngle = supportAngles
    phiangl = 2.00*pi;
    thangle = 1.0*pi;
    
    th = -innerThickness;
    npoints = 30; %no need for high level of detail in this part
    
    [phi, theta]= meshgrid(linspace(0,phiangl,npoints),linspace(0,thangle,npoints));
    r = 100 * sCalc(l,m,phi,theta);
    for ndx = 1:length(r(:))
        if abs(r(ndx)) > 10
            r(ndx) = 10;
        end
    end
    
    cx = abs(r) .* sin(theta) .* cos(phi) ;
    cy = 0 * abs(r) .* sin(theta) .* sin (phi);
    cz = abs(r) .* cos (theta);
    [r c] = size(cz);
    sth = th*ones(r,c);
    %The following code snippet just defines the range and thickness of the
    %support material. They can be adjusted and tweaked 
    for ndx1 = 1:r
        for ndx2 = 1:c
            d = (cx(ndx1,ndx2)^2 + cz(ndx1,ndx2)^2)^.5;
            if d < 10
                sth(ndx1,ndx2) = th;
            elseif d >=10
                sth(ndx1,ndx2) = 0;
            end
        end
        
    end
    for ndx = 1:numel(cx)
        Z  = [
            cos(rotAngle) -sin(rotAngle) 0;
            sin(rotAngle) cos(rotAngle) 0;
            0 0 1] * [cx(ndx);cy(ndx);cz(ndx)];
        cx(ndx) = Z(1);
        cy(ndx) = Z(2);
        cz(ndx) = Z(3);
    end
    thisSupport = surf2solid(cx,cy,cz,'triangulation','f','thickness',sth);
    sh2{end+1} = thisSupport;
    surf2solid(cx,cy,cz,'triangulation','f','thickness',sth);
end
end
function z = sCalc(l,m,phi,theta)
%%This helper function calculates the value of the spherical harmonics. In

    N = 2 *l + 1;
    N = N/ (4*pi); 
    N = N *  factorial(l-abs(m));
    N = N / factorial(l + abs(m));
    N = N^(1/2);
    
    L = legendre(l,cos(theta));
    if l > 0
        P = N * L(abs(m)+1,:,:);
    else
        P = N .* L(:,:);
    end
    P = (-1)^m * P;
    P = reshape(P(:), size(theta));
    %As per one popular convention, only the real part of the function is
    %plotted.
    if m == 0
            y2 = 1;
        elseif m > 0
            y2 = sqrt(2) * cos(m.*phi);
        else
            y2 = sqrt(2) * sin(abs(m).*phi);
    end
    z = y2.*P; 

end