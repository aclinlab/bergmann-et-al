function [x,y,z] = montageXYToStackXYZ(x_in, y_in, targetXYSize, montageDims)
% montageXYToStackXYZ: turn xy coordinates in a montage into xyz coordinates of the corresponding stack
% usage:
% [x,y,z] = montageXYToStackXYZ(x_in, y_in, targetXYSize, montageDims)
% turn a series of x and y points in a montage into a series of x y z
% coordinates
% x_in and y_in might be the outputs of getpts(h)
% targetXYSize is the size of the x and y that you're trying to get to
% Important:
% x_in is horizontal and y_in is vertical (as is the convention with
% getpts)
% But: targetXYSize is rows, columns (as Matlab normally does)
% (I know, it's counterintuitive... perhaps we will come back to this to
% fix later)

mDim = targetXYSize(1);
nDim = targetXYSize(2);
x = mod(x_in, mDim);
y = mod(y_in, nDim);
z = montageDims(2)*floor(y_in/nDim) + (floor(x_in/mDim)+1);


end