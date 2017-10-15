function  AO(n,l,m)
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

nthickness = 3; 
innerThickness = 2;
clf %clear fig

a0 = 1;
%input angle for generating points
%These angles can be changed to make a partly hollow object
border = 125;


global spoints innerThickness
spoints = 400; % This value determines the sampling size
%generate grid used for plotting and calculating points
rdist = 100;
phiangle = 2.00*pi;
thangle = 1*pi;
%[phi th r] = ngrid(linspace(0,phiangle,spoints),linspace(0,thangle,spoints), ...
 %   linspace(0,rdist,spoints));
[x y z]= meshgrid(linspace(-border,border,spoints),linspace(-border,border,spoints) ... 
    , linspace(-border,border,spoints));
r = sqrt( x.^2 + y.^2 + z.^2);
th = acos( z./r);
phi = atan(y./x);
%spherical harmonic calculations

Psi = @(r,th,phi)  sCalc(l,m,phi,th) .* Rnl(n,l,r);




psi = Psi(r,th,phi);
colors = sign(psi);
p = Psi(r,th,phi);
prob = max(p(:))/4;
x = abs(r) .* sin(th) .* cos(phi) ;
y = abs(r) .* sin(th) .* sin (phi);
z = abs(r) .* cos (th);

%uncomment following line to toggle plot
[F,V] = isosurface(abs(p),prob,colors);

daspect([1 1 1]);



daspect([1 1 1]);
view(3);
camlight('left');
lighting phong;
axis vis3d;
rotate3d on;
 xlabel('x')
 ylabel('y')
 zlabel('z')

th = nthickness;

%Call to function that makes STL file
sh = surf2solid(F,V,'thickness',th);

surf2solid(F,V,'thickness',th)


%write the file 
    s1 = sprintf('%d%d%d',n,l,m);
    s = ['H' s1 'orbital.stl'];
    stlwrite(s,sh);
%     for ndx = 1:length(sh2)
%       s = ['Y' s1 'orbitalaux' num2str(ndx) '.stl'];
%       stlwrite(s,sh2{ndx});
%     end
end
function sh2 = makeAux(l,m)
global innerThickness
sh2 = {};
if m == 0
    supportAngles = 0;
else
   supportAngles = 0:pi/m:(pi/m)*(m-1);
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
    [rc] = size(cz);
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


global spoints
%%This helper function calculates the value of the spherical harmonics. In

   N = 2 *l + 1;
    N = N/ (4*pi); 
    N = N *  factorial(l-abs(m));
    N = N / factorial(l + abs(m));
    N = N^(1/2);
    
    L = legendre(l,cos(theta));
    if l == 0
        P = N .* L;
    else
        P = N * L(abs(m)+1,:,:,:);
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

function R = Rnl(n,l,r)
a0 = 1;
Z = 1;
rho = (2*Z*r)/(n*a0);

R = rho.^l .* Lnk(n-l-1,2*l+1,rho) .* exp(-rho./2);
N = ((2*Z)/(n*a0))^3 * ((factorial(n-l-1))/(2*n*(factorial(n+l)^3)));
N =  N.^.5;
R = R*N;
end
%assosciated lagurre
function sum = Lnk(n,k,x)
%L_n^k
sum = 0;
for indx = 0:n
    thisterm = (-1)^indx .*nchoosek(n+k,n-indx) .* x.^indx .* ...
    1/(factorial(indx));
    sum = sum + thisterm;
end
end