
### needs attentio


# function [Ax,Ay,Az] = area_finder(xyz,conn);
# % This function estimates normalised 
# % net cross-section areas of the RCM 
# % relatively to the x,y,z axes. 
# % Only the path with weight 1 is considered 
# % for estimating the area to avoid ambiguity
# % in defining the edge of the area. 
import numpy as np

def area(xyz,conn)

  idx = find(conn(:,1) ~= 1);
  conn(idx,:) = [];

  X = xyz(conn(:,2),1);
  Y = xyz(conn(:,2),2);
  Z = xyz(conn(:,2),3);

  Ax = polyarea(Y,Z);
  Ay = polyarea(X,Z);
  Az = polyarea(X,Y);

  scale_factor = np.sqrt(Ax^2+Ay^2+Az^2);

  Ax = Ax/scale_factor;
  Ay = Ay/scale_factor;
  Az = Az/scale_factor;

  if 1
    if ~ispolycw(Z,Y)
        Ax = -1*Ax;

  if ~ispolycw(X,Z)
      Ay = -1*Ay;

  if ~ispolycw(Y,X)
      Az = -1*Az;











end
