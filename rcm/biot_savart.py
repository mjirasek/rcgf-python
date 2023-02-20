import numpy as np

def _Bz(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        x: float, y: float, z:float, 
        ) -> np.float64:
    return (
            (a2*b1 - b1*x + a1*(-b2 + y))*(-((c1*c2 + a1*(a2 - x) + b1*(b2 - y) \
            + np.sqrt((a1**2 + b1**2 + c1**2)*((a2 - x)**2 + (b2 - y)**2 + (c2 - \
            z)**2)) - c1*z)/(np.sqrt((a2 - x)**2 + (b2 - y)**2 + (c2 - \
            z)**2)*(b2**2*c1**2 + a2**2*(b1**2 + c1**2) - 2*b1*b2*c1*c2 + \
            b1**2*c2**2 + b1**2*x**2 + c1**2*x**2 - 2*b2*c1**2*y + 2*b1*c1*c2*y + \
            c1**2*y**2 + 2*a1*x*(b1*(b2 - y) + c1*(c2 - z)) - 2*a2*((b1**2 + \
            c1**2)*x + a1*b1*(b2 - y) + a1*c1*(c2 - z)) + a1**2*((b2 - y)**2 + \
            (c2 - z)**2) - 2*b1*(-(b2*c1) + b1*c2 + c1*y)*z + b1**2*z**2))) + \
            1/(a2**2*np.sqrt(a1**2 + b1**2 + c1**2) + b1**2*np.sqrt(a1**2 + b1**2 \
            + c1**2) + 2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2) + \
            b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + c1**2*np.sqrt(a1**2 + b1**2 + \
            c1**2) + 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*c2 + np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2**2 - 2*a2*np.sqrt(a1**2 + b1**2 + c1**2)*x + \
            np.sqrt(a1**2 + b1**2 + c1**2)*x**2 - 2*b1*np.sqrt(a1**2 + b1**2 + \
            c1**2)*y - 2*b2*np.sqrt(a1**2 + b1**2 + c1**2)*y + np.sqrt(a1**2 + \
            b1**2 + c1**2)*y**2 + a1**2*(np.sqrt(a1**2 + b1**2 + c1**2) - \
            np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) + \
            a1*(a2 - x)*(2*np.sqrt(a1**2 + b1**2 + c1**2) - np.sqrt((a1 + a2 - \
            x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) - b1**2*np.sqrt((a1 + \
            a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - b1*b2*np.sqrt((a1 \
            + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - \
            c1**2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) \
            - c1*c2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - \
            z)**2) + b1*y*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 \
            - z)**2) - 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*z - 2*np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2*z + c1*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 \
            + (c1 + c2 - z)**2)*z + np.sqrt(a1**2 + b1**2 + c1**2)*z**2))
        )

def _Bx(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        x: float, y: float, z:float, 
        ) -> np.float64:
    return (
            (b2*c1 - c1*y + b1*(-c2 + z))*((-(c1*c2) + a1*(-a2 + x) + b1*(-b2 + \
            y) - np.sqrt((a1**2 + b1**2 + c1**2)*((a2 - x)**2 + (b2 - y)**2 + (c2 \
            - z)**2)) + c1*z)/(np.sqrt((a2 - x)**2 + (b2 - y)**2 + (c2 - \
            z)**2)*(b2**2*c1**2 + a2**2*(b1**2 + c1**2) - 2*b1*b2*c1*c2 + \
            b1**2*c2**2 + b1**2*x**2 + c1**2*x**2 - 2*b2*c1**2*y + 2*b1*c1*c2*y + \
            c1**2*y**2 + 2*a1*x*(b1*(b2 - y) + c1*(c2 - z)) - 2*a2*((b1**2 + \
            c1**2)*x + a1*b1*(b2 - y) + a1*c1*(c2 - z)) + a1**2*((b2 - y)**2 + \
            (c2 - z)**2) - 2*b1*(-(b2*c1) + b1*c2 + c1*y)*z + b1**2*z**2)) + \
            1/(a2**2*np.sqrt(a1**2 + b1**2 + c1**2) + b1**2*np.sqrt(a1**2 + b1**2 \
            + c1**2) + 2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2) + \
            b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + c1**2*np.sqrt(a1**2 + b1**2 + \
            c1**2) + 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*c2 + np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2**2 - 2*a2*np.sqrt(a1**2 + b1**2 + c1**2)*x + \
            np.sqrt(a1**2 + b1**2 + c1**2)*x**2 - 2*b1*np.sqrt(a1**2 + b1**2 + \
            c1**2)*y - 2*b2*np.sqrt(a1**2 + b1**2 + c1**2)*y + np.sqrt(a1**2 + \
            b1**2 + c1**2)*y**2 + a1**2*(np.sqrt(a1**2 + b1**2 + c1**2) - \
            np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) + \
            a1*(a2 - x)*(2*np.sqrt(a1**2 + b1**2 + c1**2) - np.sqrt((a1 + a2 - \
            x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) - b1**2*np.sqrt((a1 + \
            a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - b1*b2*np.sqrt((a1 \
            + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - \
            c1**2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) \
            - c1*c2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - \
            z)**2) + b1*y*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 \
            - z)**2) - 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*z - 2*np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2*z + c1*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 \
            + (c1 + c2 - z)**2)*z + np.sqrt(a1**2 + b1**2 + c1**2)*z**2))
        )

def _By(a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        x: float, y: float, z:float, 
        ) -> np.float64:
    return (
            (-(a2*c1) + c1*x + a1*(c2 - z))*((-(c1*c2) + a1*(-a2 + x) + b1*(-b2 \
            + y) - np.sqrt((a1**2 + b1**2 + c1**2)*((a2 - x)**2 + (b2 - y)**2 + \
            (c2 - z)**2)) + c1*z)/(np.sqrt((a2 - x)**2 + (b2 - y)**2 + (c2 - \
            z)**2)*(b2**2*c1**2 + a2**2*(b1**2 + c1**2) - 2*b1*b2*c1*c2 + \
            b1**2*c2**2 + b1**2*x**2 + c1**2*x**2 - 2*b2*c1**2*y + 2*b1*c1*c2*y + \
            c1**2*y**2 + 2*a1*x*(b1*(b2 - y) + c1*(c2 - z)) - 2*a2*((b1**2 + \
            c1**2)*x + a1*b1*(b2 - y) + a1*c1*(c2 - z)) + a1**2*((b2 - y)**2 + \
            (c2 - z)**2) - 2*b1*(-(b2*c1) + b1*c2 + c1*y)*z + b1**2*z**2)) + \
            1/(a2**2*np.sqrt(a1**2 + b1**2 + c1**2) + b1**2*np.sqrt(a1**2 + b1**2 \
            + c1**2) + 2*b1*b2*np.sqrt(a1**2 + b1**2 + c1**2) + \
            b2**2*np.sqrt(a1**2 + b1**2 + c1**2) + c1**2*np.sqrt(a1**2 + b1**2 + \
            c1**2) + 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*c2 + np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2**2 - 2*a2*np.sqrt(a1**2 + b1**2 + c1**2)*x + \
            np.sqrt(a1**2 + b1**2 + c1**2)*x**2 - 2*b1*np.sqrt(a1**2 + b1**2 + \
            c1**2)*y - 2*b2*np.sqrt(a1**2 + b1**2 + c1**2)*y + np.sqrt(a1**2 + \
            b1**2 + c1**2)*y**2 + a1**2*(np.sqrt(a1**2 + b1**2 + c1**2) - \
            np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) + \
            a1*(a2 - x)*(2*np.sqrt(a1**2 + b1**2 + c1**2) - np.sqrt((a1 + a2 - \
            x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2)) - b1**2*np.sqrt((a1 + \
            a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - b1*b2*np.sqrt((a1 \
            + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) - \
            c1**2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - z)**2) \
            - c1*c2*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 - \
            z)**2) + b1*y*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 + (c1 + c2 \
            - z)**2) - 2*c1*np.sqrt(a1**2 + b1**2 + c1**2)*z - 2*np.sqrt(a1**2 + \
            b1**2 + c1**2)*c2*z + c1*np.sqrt((a1 + a2 - x)**2 + (b1 + b2 - y)**2 \
            + (c1 + c2 - z)**2)*z + np.sqrt(a1**2 + b1**2 + c1**2)*z**2))
        )

def B(
        a1: float, a2: float, 
        b1: float, b2: float, 
        c1: float, c2: float, 
        x: float, y: float, z: float, 
        J: float) -> tuple[np.ndarray]:
    return (
        J*_Bx(a1,a2,b1,b2,c1,c2,x,y,z),
        J*_By(a1,a2,b1,b2,c1,c2,x,y,z),
        J*_Bz(a1,a2,b1,b2,c1,c2,x,y,z)
    )


# print(B(1,1,1,2,5,8,0,0,0))

