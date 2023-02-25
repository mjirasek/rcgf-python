from sympy import Polygon
import numpy as np
import pandas as pd


def area_3d(xyz: pd.DataFrame) -> tuple[float]:

    _xyz = np.array(xyz.loc[::,['x','y','z']])
    _Ax = _area(_xyz[:,[1,2]])
    _Ay = _area(_xyz[:,[0,2]])
    _Az = _area(_xyz[:,[0,1]])
    scale_factor = np.sqrt(_Ax**2 + _Ay**2 + _Az**2)

    Ax = _Ax / scale_factor
    Ay = _Ay / scale_factor
    Az = _Az / scale_factor

    return Ax, Ay, Az


def _area(polygon_2d: np.array) -> float:
    try: 
        return float(Polygon(*polygon_2d).area)
        # this is very slow part, clearly the Polygon
        # is very slow function.
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