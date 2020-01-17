function AddEliptAtomicArray(xrad, yrad, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = (2*xrad - 1) * AtomSpacing;
W = (2*yrad - 1) * AtomSpacing;

%Generating all possible positions for particles evenly spaced in box
%between the langth and width bounds 
xp(1, :) = linspace(-L/2, L/2, 2*xrad);
yp(1, :) = linspace(-W/2, W/2, 2*yrad);

%Initialize the number of atoms 
%Incremeent when each possible position is withing the elipse bounds 
numAtoms = 0;

%Checking if each possible position is withing the goemoety by itterating through all the 
%possible positions and determining if its within the bounds of the shape 
for i = 1:2*xrad
    for j = 1:2*yrad
        if (X0 - xp(i))^2/xrad + (Y0 - yp(j))^2/yrad <= (AtomSpacing)^2 %yrad*xrad*pi
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms  + numAtoms) = yp(j);
        else
            i;
            j;
        end
    end
end


x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end