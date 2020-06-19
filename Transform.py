from math import sin, cos, pi
import numpy as np


class Transform:
    def __init__(self, _t=None):
        if _t is None:
            self.t = np.identity(4)
        else:
            self.t = _t

    def translateX(self, offsetx):
        temp = np.matrix([
            [1, 0, 0, offsetx],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        self.t = temp * self.t
        return self

    def translateY(self, offsety):
        temp = np.matrix([
            [1, 0, 0, 0],
            [0, 1, 0, offsety],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        self.t = temp * self.t
        return self

    def translateZ(self, offsetz):
        temp = np.matrix([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, offsetz],
            [0, 0, 0, 1],
        ])
        self.t = temp * self.t
        return self

    def rotateX(self, rotx):
        rotx = rotx * pi / 180
        temp = np.matrix([
            [1, 0, 0, 0],
            [0, cos(rotx), -sin(rotx), 0],
            [0, sin(rotx), cos(rotx), 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t
        return self

    def rotateY(self, roty):
        roty = roty * pi / 180
        temp = np.matrix([
            [cos(roty), 0, sin(roty), 0],
            [0, 1, 0, 0],
            [-sin(roty), 0, cos(roty), 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t
        return self

    def rotateZ(self, rotz):
        rotz = rotz * pi / 180
        temp = np.matrix([
            [cos(rotz), -sin(rotz), 0, 0],
            [sin(rotz), cos(rotz), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t
        return self

    # DH = Rx(alfa) * Tx(a) * Rz(teta) * Tz(d)
    # if we have mDH = mDH(T1) * mDH(T2) * mDH(T3) we should use like this:
    # t = Transform()
    # t.mDH(T3).mDH(T2).mDH(T1)
    def mDH(self, alpha, a, d, theta):
        self.translateZ(d).rotateZ(theta).translateX(a).rotateX(alpha)
        return self

    def constantXYZ(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma)-sin(alpha) *
             cos(gamma), cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), 0],
            [sin(alpha)*cos(beta), sin(alpha)*sin(beta)*sin(gamma)+cos(alpha) *
             cos(gamma), sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), 0],
            [-sin(beta), cos(beta)*sin(gamma), cos(beta)*cos(gamma), 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t

        return self

    def constantXZY(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta), -cos(alpha)*sin(beta)*cos(gamma)+sin(alpha) *
             sin(gamma), cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma), 0],
            [sin(beta), cos(beta)*cos(gamma), -cos(beta)*sin(gamma), 0],
            [-sin(alpha) * cos(beta), sin(alpha)*sin(beta) * cos(gamma) + cos(alpha) *
             sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha) * cos(gamma), 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t

        return self

    def constantyxz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma), -sin(alpha)*cos(beta),
             sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma), 0],
            [cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma), cos(alpha)*cos(beta), -cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), 0],
            [-cos(beta)*sin(gamma), sin(beta),
            cos(beta)*cos(gamma), 0],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t

        return self

    def constantyzx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta)*cos(gamma),
            -sin(beta), cos(beta)*sin(gamma),
            0], [cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),
            cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma), 0], 
            [sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), sin(alpha)*cos(beta), sin(alpha)*sin(beta)*sin(gamma), 0 ]
            ,[0, 0, 0, 1]
        ])
        self.t = temp * self.t

        return self

    def constantzxy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
             [sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma), sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), cos(alpha)*cos(beta), 0], 
            [cos(beta)*sin(gamma), 
             cos(beta)*cos(alpha), -sin(beta), 0], 
            [cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),
             cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), cos(alpha)*cos(beta), 0],
            [0, 0, 0 ,1]
        ])
        self.t = temp * self.t

        return self
    
    def constantzyx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta =  beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta)*cos(gamma) ,
            -cos(beta)*sin(gamma) ,
            sin(beta) , 0],
            [sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma) ,-sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma) ,
            -sin(alpha)*sin(beta) , 0],[-cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma) ,
            cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma), 
            cos(alpha)*cos(beta), 0 ] ,
            [0 , 0, 0,1]
        ])
        self.t = temp * self.t

        return self

    def constantxyx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta) ,sin(beta)*sin(gamma) ,sin(beta)*cos(gamma) , 0],
            [sin(alpha)*sin(beta), 
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma), 
            -sin(alpha)*cos(beta)*cos(gamma)-cos(alpha)*sin(gamma) , 0], 
            [-cos(alpha)*sin(beta), 
            cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma) 
            ,cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) , 0],
            [0 , 0, 0 , 1]
        ])
        self.t = temp * self.t

        return self

    def constantxzx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta) ,-sin(beta)*cos(gamma), 
            sin(beta)*sin(gamma), 0],
            [cos(alpha)*sin(beta) ,cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma), 
            -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma) , 0],
            [sin(alpha)*sin(beta) ,sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma), 
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma), 0 ],
            [0, 0, 0, 1]
        ])
        self.t = temp * self.t

        return self

    def constantyxy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma) ,sin(alpha)*sin(beta) ,sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma), 0 ] ,
            [sin(beta)*sin(gamma) ,cos(beta), 
            -sin(beta)*cos(gamma), 0 ],[-cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma) ,
            cos(alpha)*sin(beta) ,
            cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma), 0] ,
            [ 0 , 0 , 0 , 1 ]
        ])
        self. t = temp * self.t

        return self

    def constantyzy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma), -cos(alpha)*sin(beta), 
            cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma), 0 ],
            [sin(beta)*cos(gamma) ,cos(beta) ,sin(beta)*sin(gamma) , 0 ],
            [-sin(alpha)*cos(beta)*cos(gamma)-cos(alpha)*sin(gamma) ,sin(alpha)*sin(beta) ,
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ] 
        ])
        self.t = temp * self.t

        return self

    def constantzxz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma), 
            -sin(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(beta), 0 ],
            [cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma) ,cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) ,-cos(alpha)*sin(beta), 0 ],
            [sin(beta)*sin(gamma) ,sin(beta)*cos(gamma), 
            cos(beta),0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def constantzyz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) ,
            -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma), 
            cos(alpha)*sin(beta), 0 ],
            [sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma), 
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma), 
            sin(alpha)*sin(beta) , 0 ],
            [-sin(beta)*cos(gamma) ,sin(beta)*sin(gamma) ,cos(beta),0],
            [ 0 , 0 , 0 , 1 ] 
        ])
        self.t = temp * self.t 

        return self


    def eulerxyz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta)*cos(gamma) ,-cos(beta)*sin(gamma) ,sin(beta), 0 ],
            [sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma),
            -sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),-sin(alpha)*cos(beta), 0 ],
            [-cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),
            cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma),
            cos(alpha)*cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t 

        return self

    def eulerxzy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta)*cos(gamma)  ,-sin(beta) ,cos(beta)*sin(gamma), 0 ],
            [cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),
            cos(alpha)*cos(beta),cos(alpha)* sin(beta)*sin(gamma)-sin(alpha)*cos(gamma), 0 ],
            [sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma),
            sin(alpha)*cos(beta),sin(alpha)*sin(beta)*sin(gamma), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def euleryxz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma) ,
            sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma) ,sin(alpha)*cos(beta), 0 ],
            [cos(beta)*sin(gamma),cos(beta)*cos(alpha),-sin(beta), 0 ],
            [cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),
            cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*sin(gamma),cos(alpha)*cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def euleryzx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta),-cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),
            cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma), 0 ],
            [sin(beta),cos(beta)*cos(gamma),-cos(beta)*sin(gamma), 0 ],
            [-sin(alpha)*cos(beta),sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma),
            -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), 0 ],
            [ 0 , 0 , 0 ,1 ]
        ])
        self.t = temp * self.t

        return self

    def eylerzxy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),-sin(alpha)*cos(beta),
            sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma), 0 ],
            [cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma) ,cos(alpha)*cos(beta),
            -cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), 0 ],
            [-cos(beta)*sin(gamma),sin(beta),cos(beta)*cos(gamma), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def eulerzyx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp=np.matrix([
            [cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma), 
            cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), 0 ], 
            [sin(alpha)*cos(beta),sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma), 
            sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), 0 ], 
            [-sin(beta), cos(beta)*sin(gamma), cos(beta)*cos(gamma), 0],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def eulerxyx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(beta),sin(beta)*sin(gamma),sin(beta)*cos(gamma), 0 ],
            [sin(alpha)*sin(beta),-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma),
            -sin(alpha)*cos(beta)*cos(gamma)-cos(alpha)*sin(gamma), 0 ],
            [-cos(alpha)*sin(beta),cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma),
            cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def eulerxzx(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180
        
        temp = np.matrix([
            [cos(beta), -sin(beta)*cos(gamma),sin(beta)*sin(gamma), 0 ],
            [cos(alpha)*sin(beta),cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma),
            -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma), 0 ],
            [sin(alpha)*sin(beta),sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma),
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def euleryxy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma),sin(alpha)*sin(beta),
            sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma), 0 ],
            [sin(beta)*sin(gamma),cos(beta),-sin(beta)*cos(gamma), 0 ],
            [-cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma),cos(alpha)*sin(beta),
            cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) , 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def euleryzy(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma),-cos(alpha)*sin(beta),
            cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma), 0 ],
            [sin(beta)*cos(gamma),cos(beta),sin(beta)*sin(gamma), 0 ],
            [-sin(alpha)*cos(beta)*cos(gamma)-cos(alpha)*sin(gamma),
            sin(alpha)*sin(beta) ,-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def eulerzxz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma),
            -sin(alpha)*cos(beta)*cos(gamma)-cos(alpha)*sin(gamma),sin(alpha)*sin(beta), 0 ],
            [cos(alpha)*cos(beta)*sin(gamma)+sin(alpha)*cos(gamma),
            cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma),-cos(alpha)*sin(beta), 0 ],
            [sin(beta)*sin(gamma),sin(beta)*cos(gamma),cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

    def eulerzyz(self, alpha, beta, gamma):
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180

        temp = np.matrix([
            [cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma),
            -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma),cos(alpha)*sin(beta), 0 ],
            [sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma),
            -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma),sin(alpha)*sin(beta), 0 ],
            [-sin(beta)*cos(gamma),sin(beta)*sin(gamma),cos(beta), 0 ],
            [ 0 , 0 , 0 , 1 ]
        ])
        self.t = temp * self.t

        return self

            


#r = Transform()

#r.rotateX(30)

...

#r.rotateY(10)
