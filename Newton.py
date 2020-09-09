import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

class fun_u:
    def __init__(object,x,y):
        object.x=x
        object.y=y

    def ordinary(abc):
        x=abc.x
        y=abc.y
        return x**4 - 6*x**2*y**2 + y**4 - 1

    def derivative_x(abc):
        x=abc.x
        y=abc.y
        return 4*x**3 - 12*x*y**2

    def derivative_y(abc):
        x=abc.x
        y=abc.y
        return 4*y**3 - 12*x**2*y

class fun_v:
    def __init__(object,x,y):
        object.x=x
        object.y=y
    
    def ordinary(abc):
        x=abc.x
        y=abc.y
        return 4*x**3*y - 4*x*y**3
    
    def derivative_x(abc):
        x=abc.x
        y=abc.y
        return 12*x**2*y - 4*y**3
    
    def derivative_y(abc):
        x=abc.x
        y=abc.y
        return 4*x**3 - 12*x*y**2

def color_decision(x_k):
    if np.dot((x_k-np.array([1,0])),(x_k-np.array([1,0])))<0.0001:
        return 1
    elif np.dot((x_k-np.array([0,1])),(x_k-np.array([0,1])))<0.001:
        return 2
    elif np.dot((x_k-np.array([-1,0])),(x_k-np.array([-1,0])))<0.0001:
        return 3
    elif np.dot((x_k-np.array([0,-1])),(x_k-np.array([0,-1])))<0.0001:
        return 4
    else:
        print("to many iteration at",x_k)
        return 10
#
#def fun_u(x,y):
#    return x**4 - 6*x**2*y**2 + y**4 - 1
#
#def fun_v(x,y):
#    return 4*x**3*y - 4*x*y**3
#
#def fun_u_x(x,y):
#    return 4*x**3 - 12*x*y**2
#
#def fun_u_y(x,y):
#    return 4*y**3 - 12*x**2*y
#
#def fun_v_x(x,y):
#    return 12*x**2*y - 4*y**3
#
#def fun_v_y(x,y):
#    return 4*x**3 - 12*x*y**2

def solve_by_newton(class_1,class_2,z_k): #make sure z_k is np.array
    z_k=np.array(z_k)
    global s
    # print(s)
    #print("z_k",z_k)
    u=class_1(*z_k)
    v=class_2(*z_k)

    uo=u.ordinary()
    ux=u.derivative_x()
    uy=u.derivative_y()
    vo=v.ordinary()
    vx=v.derivative_x()
    vy=v.derivative_y()

    f_f_prime_ratio=np.array([uo*ux+vo*vx,uo*uy+vo*vy])
    f_f_prime_ratio=f_f_prime_ratio/(ux**2+vx**2)
    #print("f/f_prime",f_f_prime_ratio)
    if np.dot(f_f_prime_ratio,f_f_prime_ratio)<0.00000001:
        return z_k-f_f_prime_ratio
    elif s>800:
        s = 0
        return np.array([9,10])
    else:
        s=s+1
        return solve_by_newton(class_1,class_2,z_k-f_f_prime_ratio)

def solve_by_newton_2(class_1,class_2,z_k,root):
    global s
    global log_x
    global ratio
    z_k=np.array(z_k)
    root=np.array(root)
    log_x.append(np.log10(np.sqrt(np.dot(z_k-root,z_k-root))))
    
    print("z_k",z_k)
    u=class_1(*z_k)
    v=class_2(*z_k)
    
    uo=u.ordinary()
    ux=u.derivative_x()
    uy=u.derivative_y()
    vo=v.ordinary()
    vx=v.derivative_x()
    vy=v.derivative_y()
    
    f_f_prime_ratio=np.array([uo*ux+vo*vx,uo*uy+vo*vy])
    f_f_prime_ratio=f_f_prime_ratio/(ux**2+vx**2)
    print("f/f_prime",f_f_prime_ratio)
    z_k_1=z_k-f_f_prime_ratio
    ratio.append(np.sqrt(np.dot(z_k_1-root,z_k_1-root))/np.dot(z_k-root,z_k-root))
    #if np.dot(f_f_prime_ratio,f_f_prime_ratio)<0.0000000000001:
    print("go into if")
    print("s",s)
    if s>1000:
        return np.array([9,10])
    else:
        s=s+1
        return solve_by_newton_2(class_1,class_2,z_k_1,root)


def solve_by_newton_limited_step(class_1,class_2,z_k,delta):
    global s
    z_k=np.array(z_k)
    #print("z_k",z_k)
    u=class_1(*z_k)
    v=class_2(*z_k)
    
    uo=u.ordinary()
    ux=u.derivative_x()
    uy=u.derivative_y()
    vo=v.ordinary()
    vx=v.derivative_x()
    vy=v.derivative_y()
    
    f_f_prime_ratio=np.array([uo*ux+vo*vx,uo*uy+vo*vy])
    f_f_prime_ratio=f_f_prime_ratio/(ux**2+vx**2)
    #print("f/f_prime",f_f_prime_ratio)
    
    if np.dot(f_f_prime_ratio,f_f_prime_ratio)<0.00000001:
        return z_k-f_f_prime_ratio
    elif s>200:
        return np.array([9,10])
    else:
        z_k_1=z_k-f_f_prime_ratio
        s=s+1
        if np.sqrt(np.dot(f_f_prime_ratio,f_f_prime_ratio))>delta:
            z_k_1=z_k-delta*f_f_prime_ratio/np.sqrt(np.dot(f_f_prime_ratio,f_f_prime_ratio))
            return solve_by_newton_limited_step(class_1,class_2,z_k_1,delta)
        else:
            return solve_by_newton_limited_step(class_1,class_2,z_k_1,delta)
    


def newton_basin():
    global s
    x=np.linspace(-2.0,2.0,num=41)
    y=np.linspace(-2.0,2.0,num=41)
    data_x=[]
    data_y=[]
    value=[]
    for i in range(len(x)):
        for j in range(len(y)):
            print("z_k",x[i],y[j])
            data_x.append(x[i])
            data_y.append(y[j])
            s=0
            result=solve_by_newton(fun_u,fun_v,[x[i],y[j]])
            print("recusion time to find solution:",s)
            #result=solve_by_newton_limited_step(fun_u,fun_v,[x[i],y[j]],0.5)
            value.append(color_decision(result))

    data_x=np.array(data_x)
    data_y=np.array(data_y)
    value=np.array(value)
    grid=value.reshape((len(x),len(y)))
    plt.imshow(grid, extent=(data_x.min(),data_x.max(),data_y.max(),data_y.min()),interpolation='nearest',cmap=plt.get_cmap('Set1'))
    plt.show()

def recursion_time():
    global s
    x=np.linspace(-1.2,-0.8,num=400)
    y=np.linspace(-1.2,-0.8,num=400)
    data_x=[]
    data_y=[]
    value=[]
    for i in range(len(x)):
        for j in range(len(y)):
            print("z_k",x[i],y[j])
            data_x.append(x[i])
            data_y.append(y[j])
            s=0
            result=solve_by_newton(fun_u,fun_v,[x[i],y[j]])
            value.append(s)

    data_x=np.array(data_x)
    data_y=np.array(data_y)
    value=np.array(value)
    grid=value.reshape((len(x),len(y)))
    plt.imshow(grid, extent=(data_x.min(),data_x.max(),data_y.max(),data_y.min()),interpolation='nearest',cmap=plt.get_cmap('Spectral'))
    plt.colorbar()
    plt.show()

def convergence_rate():
    global s
    global log_x
    global ratio
    ratio=[]
    log_x=[]
    s=0
    root=[1,0]
    result=solve_by_newton_2(fun_u,fun_v,[5,3],[1,0])
    print("result is",result)
    print("log_x",log_x)
    x=np.arange(len(log_x))
    plt.plot(x,log_x,'o')
    #plt.plot(x,ratio,'ro')
    plt.show()

s=0
# result=solve_by_newton(fun_u,fun_v,[2.0,2.0])
# print("result is",result)

newton_basin()
# recursion_time()
#convergence_rate()

















#
#def test(class_1,class_2,x_k):
#    print(class_1.ordinary(),class_1.derivative_x(),class_1.derivative_y())
#    print(class_2.ordinary(),class_2.derivative_x(),class_2.derivative_y())
#
#
#u=fun_u
#print(u.x)
#v=fun_v
#
#test(u,v,[2,5])
