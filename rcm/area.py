from sympy import Polygon
import numpy as np


def area_3d(polygon_xyz: np.array) -> float:

    Ax = area(polygon_xyz[:,1:3]);
    Ay = area(polygon_xyz[:,0:3:2]);
    Az = area(polygon_xyz[:,0:2]);

    return Ax,Ay,Az



def area(polygon_2d: np.array) -> float:
    return float(Polygon(*polygon_2d).area)






# X = xyz(conn(:,2),1);
# Y = xyz(conn(:,2),2);
# Z = xyz(conn(:,2),3);

# Ax = polyarea(Y,Z);
# Ay = polyarea(X,Z);
# Az = polyarea(X,Y);

# scale_factor = sqrt(Ax^2+Ay^2+Az^2);

# Ax = Ax/scale_factor;
# Ay = Ay/scale_factor;
# Az = Az/scale_factor;


# if 1
# if ~ispolycw(Z,Y)
#     Ax = -1*Ax;
# end

# if ~ispolycw(X,Z)
#     Ay = -1*Ay;
# end

# if ~ispolycw(Y,X)
#     Az = -1*Az;
# end
# end