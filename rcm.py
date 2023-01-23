import numpy as np

def _Bz(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        t: float, 
        x: float, y: float, z:float, 
        J) -> np.float64:
    return (
        (J*1*(a2*b1 - b1*x + a1*(-b2 + y)))/(4.*np.pi*(a2**2*np.sqrt(a1**2 \
        + b1**2 + c1**2) + b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + \
        np.sqrt(a1**2 + b1**2 + c1**2)*c2**2 + 2*a1*a2*np.sqrt(a1**2 + b1**2 \
        + c1**2)*t + 2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2)*t + \
        2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*c2*t + a1**2*np.sqrt(a1**2 + \
        b1**2 + c1**2)*t**2 + b1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 + \
        c1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 - 2*a2*np.sqrt(a1**2 + \
        b1**2 + c1**2)*x - 2*a1*np.sqrt(a1**2 + b1**2 + c1**2)*t*x + \
        np.sqrt(a1**2 + b1**2 + c1**2)*x**2 - 2*b2*np.sqrt(a1**2 + b1**2 + \
        c1**2)*y - 2*b1*np.sqrt(a1**2 + b1**2 + c1**2)*t*y + np.sqrt(a1**2 + \
        b1**2 + c1**2)*y**2 - 2*np.sqrt(a1**2 + b1**2 + c1**2)*c2*z - \
        2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*t*z + np.sqrt(a1**2 + b1**2 + \
        c1**2)*z**2 - a1*a2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - b1*b2*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - c1*c2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - a1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - b1**2*t*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - c1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + a1*x*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) + b1*y*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + c1*z*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2)))
    )

def _Bx(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        t: float, 
        x: float, y: float, z:float, 
        J) -> np.float64:
    return (
        (J*1*(b2*c1 - c1*y + b1*(-c2 + z)))/(4.*np.pi*(a2**2*np.sqrt(a1**2 \
        + b1**2 + c1**2) + b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + \
        np.sqrt(a1**2 + b1**2 + c1**2)*c2**2 + 2*a1*a2*np.sqrt(a1**2 + b1**2 \
        + c1**2)*t + 2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2)*t + \
        2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*c2*t + a1**2*np.sqrt(a1**2 + \
        b1**2 + c1**2)*t**2 + b1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 + \
        c1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 - 2*a2*np.sqrt(a1**2 + \
        b1**2 + c1**2)*x - 2*a1*np.sqrt(a1**2 + b1**2 + c1**2)*t*x + \
        np.sqrt(a1**2 + b1**2 + c1**2)*x**2 - 2*b2*np.sqrt(a1**2 + b1**2 + \
        c1**2)*y - 2*b1*np.sqrt(a1**2 + b1**2 + c1**2)*t*y + np.sqrt(a1**2 + \
        b1**2 + c1**2)*y**2 - 2*np.sqrt(a1**2 + b1**2 + c1**2)*c2*z - \
        2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*t*z + np.sqrt(a1**2 + b1**2 + \
        c1**2)*z**2 - a1*a2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - b1*b2*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - c1*c2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - a1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - b1**2*t*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - c1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + a1*x*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) + b1*y*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + c1*z*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2)))
    )

def _By(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        t: float, 
        x: float, y: float, z:float, 
        J) -> np.float64:
    return (
        -0.25*(J*1*(a2*c1 - c1*x + a1*(-c2 + \
        z)))/(np.pi*(a2**2*np.sqrt(a1**2 + b1**2 + c1**2) + \
        b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + np.sqrt(a1**2 + b1**2 + \
        c1**2)*c2**2 + 2*a1*a2*np.sqrt(a1**2 + b1**2 + c1**2)*t + \
        2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2)*t + 2*c1*np.sqrt(a1**2 + b1**2 \
        + c1**2)*c2*t + a1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 + \
        b1**2*np.sqrt(a1**2 + b1**2 + c1**2)*t**2 + c1**2*np.sqrt(a1**2 + \
        b1**2 + c1**2)*t**2 - 2*a2*np.sqrt(a1**2 + b1**2 + c1**2)*x - \
        2*a1*np.sqrt(a1**2 + b1**2 + c1**2)*t*x + np.sqrt(a1**2 + b1**2 + \
        c1**2)*x**2 - 2*b2*np.sqrt(a1**2 + b1**2 + c1**2)*y - \
        2*b1*np.sqrt(a1**2 + b1**2 + c1**2)*t*y + np.sqrt(a1**2 + b1**2 + \
        c1**2)*y**2 - 2*np.sqrt(a1**2 + b1**2 + c1**2)*c2*z - \
        2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*t*z + np.sqrt(a1**2 + b1**2 + \
        c1**2)*z**2 - a1*a2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - b1*b2*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - c1*c2*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - a1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) - b1**2*t*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) - c1**2*t*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + a1*x*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2) + b1*y*np.sqrt(a2**2 + \
        b2**2 + c2**2 + (a1**2 + b1**2 + c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y \
        + y**2 + 2*t*(a1*(a2 - x) + b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + \
        z**2) + c1*z*np.sqrt(a2**2 + b2**2 + c2**2 + (a1**2 + b1**2 + \
        c1**2)*t**2 - 2*a2*x + x**2 - 2*b2*y + y**2 + 2*t*(a1*(a2 - x) + \
        b1*(b2 - y) + c1*(c2 - z)) - 2*c2*z + z**2)))
    )

def B(a1,a2,b1,b2,c1,c2,x,y,z,J) -> tuple[np.ndarray]:
    return (
        (
        _Bx(a1,a2,b1,b2,c1,c2,1,x,y,z,J) - 
        _Bx(a1,a2,b1,b2,c1,c2,0,x,y,z,J)
        ),
        (
        _By(a1,a2,b1,b2,c1,c2,1,x,y,z,J) - 
        _By(a1,a2,b1,b2,c1,c2,0,x,y,z,J)
        ),
        (
        _Bz(a1,a2,b1,b2,c1,c2,1,x,y,z,J) - 
        _Bz(a1,a2,b1,b2,c1,c2,0,x,y,z,J)
        )
    )




