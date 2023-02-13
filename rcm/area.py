from sympy import Polygon
import numpy as np
# from typing import Tuple


def area_3d(xyz: np.array) -> tuple[float]:

    Ax = _area(xyz[:,[1,2]])
    Ay = _area(xyz[:,[0,2]])
    Az = _area(xyz[:,[0,1]])
    return Ax,Ay,Az


def _area(polygon_2d: np.array) -> float:
    try: 
        return float(Polygon(*polygon_2d).area)
    except:
        return 10




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