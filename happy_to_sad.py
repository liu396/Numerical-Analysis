from graphics import *
import matplotlib.pyplot as plt
import time
import numpy as np
from math import acos,sqrt,atan,cos,sin,pi

def lagrange_interpolate(x,x_value,y_value):
    # This function can calculate the value of Pn(x)
    assert len(x_value)==len(y_value)
    n=len(x_value)
    Pn_x=0
    for k in range(0,n):
        basis=1
        nominator=1
        for j in range(0,n):
            if(j!=k):
                basis=basis*(x_value[k]-x_value[j])
                nominator=nominator*(x-x_value[j])
                l=nominator/basis
        Pn_x=Pn_x+l*y_value[k]
    
    return Pn_x

def divided_difference(x_value,y_value):
    if len(x_value)==1:
        return y_value[0]
    else:
        factor=1/(x_value[0]-x_value[-1])
        x_save=x_value.copy()
        y_save=y_value.copy()
        x_save_2=x_value.copy()
        y_save_2=y_value.copy()
        x_save_2.pop(0)
        y_save_2.pop(0)
        x_save.pop()
        y_save.pop()
        return factor*(divided_difference(x_save,y_save)-divided_difference(x_save_2,y_save_2))

def phi(x,x_value,k):
    #calculate the multiplication Term
    if k==-1:
        return 1
    else:
        return phi(x,x_value,k-1)*(x-x_value[k])

def newton_interpolate(x,x_value,y_value):
    Pn_x=0
    assert len(x_value)==len(y_value)
    n=len(x_value)
    #print("n",n)
    x_in=[]
    y_in=[]
    for k in range(n):
        x_in.append(x_value[k])
        x_new=x_in.copy()
        y_in.append(y_value[k])
        y_new=y_in.copy()
        #print("xin",x_in)
        #print("yin",y_in)
        Pn_x=Pn_x+divided_difference(x_new,y_new)*phi(x,x_value,k-1)
    
    return Pn_x

def s_coefficient(m,x_value,y_value):
    #To calculate alpha[]
    x_full_value=x_value.copy()
    x_full_value.append(8.3)#x_n+1
    x_full_value.append(8.6)#x_n+2
    x_full_value.append(8.9)#x_n+3
    x_full_value.append(-0.9)#x_-3
    x_full_value.append(-0.6)#x_-2
    x_full_value.append(-0.3)#x_-1
    A=[]
    print(x_full_value)
    
    #take consideration of the left boundary condition
    row=[]
    for k in range(-m+1,len(x_value)):
        x_aux=[]
        for i in range(k+m,k-2,-1):
            x_aux.append(x_full_value[i])
        #print(x_aux)
        row.append(divided_difference_t_x_prime(x_full_value[0],x_aux,m))
    A=row
#middle rows
    for line in range(0,len(x_value)):
        row=[]
        for k in range(-m+1,len(x_value)):
            x_aux=[]
            for i in range(k+m,k-2,-1):
                x_aux.append(x_full_value[i])
            row.append(divided_difference_t_x(x_full_value[line],x_aux,m))
        A=np.vstack([A,row])

#take consideration of the right boudary condition
    row=[]
    for k in range(-m+1,len(x_value)):
        x_aux=[]
        for i in range(k+m,k-2,-1):
            x_aux.append(x_full_value[i])
        row.append(divided_difference_t_x_prime(x_full_value[len(x_value)-1],x_aux,m))
    A=np.vstack([A,row])
    #print(A)
    B=y_value.copy()
    B.insert(0,0)#x'0
    B.insert(len(B),0)#x'n
    print('B',B)
    alpha=np.linalg.solve(A,B)
    alpha_move=[0]*len(alpha)
    
    for w in range(len(alpha)):
        alpha_move[w-m+1]=alpha[w]

    print(x_full_value)
    print('coefficients generated')
    return [alpha_move,x_full_value]

def s_coefficient_linear(m,x_value,y_value):
    #To calculate alpha[]
    x_full_value=x_value.copy()
    x_full_value.append(8.3)#x_n+1
    x_full_value.append(-0.3)#x_-1
    A=[]
    print(x_full_value)
    
    row=[]
    #middle rows
    for k in range(-m+1,len(x_value)):
        x_aux=[]
        for i in range(k+m,k-2,-1):
            x_aux.append(x_full_value[i])
        row.append(divided_difference_t_x(x_full_value[0],x_aux,m))
    A=row
    
    for line in range(1,len(x_value)):
        row=[]
        for k in range(-m+1,len(x_value)):
            x_aux=[]
            for i in range(k+m,k-2,-1):
                x_aux.append(x_full_value[i])
            row.append(divided_difference_t_x(x_full_value[line],x_aux,m))
        A=np.vstack([A,row])

    #print(A)
    B=y_value.copy()
    print('B',B)
    alpha=np.linalg.solve(A,B)
    
    print(x_full_value)
    print('coefficients generated')
    return [alpha,x_full_value]



def s(t,m,x_full_value,alpha):# need to add extra x
    result=0
    
    for k in range(-m+1,len(x_full_value)-6):
        x_aux=[]
        for i in range(k+m,k-2,-1):
            x_aux.append(x_full_value[i])
        phik=divided_difference_t_x(t,x_aux,m)
        result=result+alpha[k]*phik
    return result

def s_linear(t,m,x_full_value,alpha):# need to add extra x
    result=0
    
    for k in range(-m+1,len(x_full_value)-2):
        x_aux=[]
        for i in range(k+m,k-2,-1):
            x_aux.append(x_full_value[i])
        phik=divided_difference_t_x(t,x_aux,m)
        result=result+alpha[k]*phik
    return result

def divided_difference_t_x(t,x_aux,m):
    if len(x_aux)==1:
        if t>x_aux[0]:
            return (t-x_aux[0])**m
        else:
            return 0
    else:
        factor=1/(x_aux[0]-x_aux[-1])
        x_left=x_aux.copy()
        x_right=x_aux.copy()
        x_left.pop()
        x_right.pop(0)
        return factor*(divided_difference_t_x(t,x_left,m)-divided_difference_t_x(t,x_right,m))

def divided_difference_t_x_prime(t,x_aux,m):
    if len(x_aux)==1:
        if t>x_aux[0]:
            return m*(t-x_aux[0])**(m-1)
        else:
            return 0
    else:
        factor=1/(x_aux[0]-x_aux[-1])
        x_left=x_aux.copy()
        x_right=x_aux.copy()
        x_left.pop()
        x_right.pop(0)
        return factor*(divided_difference_t_x_prime(t,x_left,m)-divided_difference_t_x_prime(t,x_right,m))


def quadrant(x,y):
    if x>=0 and y>=0:
        return 1
    elif x<0 and y>=0:
        return 2
    elif x<0 and y<=0:
        return 3
    else:
        return 4

def angle(x,y,r):
    if quadrant(x,y)==1 or quadrant(x,y)==2:
        return acos(x/r)
    else:
        return 2*pi-acos(x/r)

def link_two(center,r,xa,ya,xb,yb,angle_a,angle_b,arc):
    dtheta=10/r
    if abs(angle_a-angle_b)<pi:
        if angle_b>angle_a:
            theta_start=angle_a
            theta_end=angle_b
        else:
            theta_start=angle_b
            theta_end=angle_a
    else:
        if angle_a<angle_b:
            theta_start=angle_b-2*pi
            theta_end=angle_a
        else:
            theta_start=angle_a-2*pi
            theta_end=angle_b

    i=theta_start
    while i<theta_end:
            line=Line(Point(r*cos(i)+center[0],r*sin(i)+center[1]),Point(r*cos(i+dtheta)+center[0],r*sin(i+dtheta)+center[1]))
            #print(line)
            arc.append(line)
            i=i+dtheta
#   print(arc)



def draw_arc(p1,p2,p3,win):
    #win.autoflush=False
    arc=[]
    xa=p1.getX()
    ya=p1.getY()
    xb=p2.getX()
    yb=p2.getY()
    xc=p3.getX()
    yc=p3.getY()
    a11=2*(xa-xb)
    a12=2*(ya-yb)
    a21=2*(xa-xc)
    a22=2*(ya-yc)
    b1=-(xb**2+yb**2-xa**2-ya**2)
    b2=-(xc**2+yc**2-xa**2-ya**2)

    a=np.array([[a11,a12],[a21,a22]])
    b=np.array([b1,b2])
    det=np.linalg.det(a)
    if abs(det)<0.1:
        arc.append(Line(p1,p2))
        arc.append(Line(p2,p3))
    else:
        center=np.linalg.solve(a,b)
        #pt=Point(center[0],center[1])
        #pt.draw(win)
        r=sqrt((xa-center[0])**2+(ya-center[1])**2)
        #print(r)
        theta_1=angle(xa-center[0],ya-center[1],r)
        theta_2=angle(xb-center[0],ya-center[1],r)
        theta_3=angle(xc-center[0],ya-center[1],r)
        
        link_two(center,r,xa,ya,xb,yb,theta_1,theta_2,arc)
        link_two(center,r,xb,yb,xc,yc,theta_2,theta_3,arc)
    
    for j in range(len(arc)):
        arc[j].draw(win)
    #update()
    #win.autoflush=True

    return arc

def erase(smth,win):
    win.autoflush=False
    for i in range(len(smth)):
        smth[i].undraw()

    del smth

def create_control_points_mouth_runge(number,start_1,end_1,start_2,end_2,start_3,end_3):
    control_points=[]
    amplitude_1=(end_1-start_1)
    amplitude_2=(end_2-start_2)
    amplitude_3=(end_3-start_3)
    delta_1=(end_1-start_1)/(number-1)
    delta_2=(end_2-start_2)/(number-1)
    delta_3=(end_3-start_3)/(number-1)
    for i in range(number):
        control_points.append([Point(-100,start_1+amplitude_1*(1/(1+25*(i/(number-1)-0.5)**2))),Point(-50,start_2+amplitude_2*(1/(1+25*(i/(number-1)-0.5)**2))),Point(0,start_3+amplitude_3*(1/(1+25*(i/(number-1)-0.5)**2))),Point(50,start_2+amplitude_2*(1/(1+25*(i/(number-1)-0.5)**2))),Point(100,start_1+amplitude_1*(1/(1+25*(i/(number-1)-0.5)**2)))])
    return control_points

def create_control_points_righteye_runge(number,start_1,end_1,start_2,end_2,start_3,end_3):
    control_points=[]
    amplitude_1=(end_1-start_1)
    amplitude_2=(end_2-start_2)
    amplitude_3=(end_3-start_3)
    delta_1=(end_1-start_1)/(number-1)
    delta_2=(end_2-start_2)/(number-1)
    delta_3=(end_3-start_3)/(number-1)
    for i in range(number):
        control_points.append([Point(70,start_1+amplitude_1*(1/(1+25*(i/number-0.5)**2))),Point(120,start_2+amplitude_2*(1/(1+25*(i/number-0.5)**2))),Point(170,start_3+amplitude_3*(1/(1+25*(i/number-0.5)**2)))])
    return control_points

def create_control_points_lefteye_runge(number,start_1,end_1,start_2,end_2,start_3,end_3):
    control_points=[]
    amplitude_1=(end_1-start_1)
    amplitude_2=(end_2-start_2)
    amplitude_3=(end_3-start_3)
    delta_1=(end_1-start_1)/(number-1)
    delta_2=(end_2-start_2)/(number-1)
    delta_3=(end_3-start_3)/(number-1)
    for i in range(number):
        control_points.append([Point(-170,start_1+amplitude_1*(1/(1+25*(i/number-0.5)**2))),Point(-120,start_2+amplitude_2*(1/(1+25*(i/number-0.5)**2))),Point(-70,start_3+amplitude_3*(1/(1+25*(i/number-0.5)**2)))])
    return control_points

def create_control_points_mouth_linear(number,start_1,end_1,start_2,end_2,start_3,end_3):
    control_points=[]
    delta_1=(end_1-start_1)/number
    delta_2=(end_2-start_2)/number
    delta_3=(end_3-start_3)/number
    for i in range(number+1):
        control_points.append([Point(-100,start_1+i*delta_1),Point(-50,start_2+i*delta_2),Point(0,start_3+i*delta_3),Point(50,start_2+i*delta_2),Point(100,start_1+i*delta_1)])
    return control_points

def create_time_series(start,stop,frames):
    t=[]
    d=(stop-start)/(frames-1)
    for i in range(frames):
        t.append(start+i*d)
    return t


def polynomial():
    n_points=int(input("Enter total number of control frames"))
    win = GraphWin('Happy to sad', 300, 300)
    print(win)
    win.setCoords(-300,-300,300,300)
    
    #draw face (invariant)====
    face=Circle(Point(0,0),250)
    face.setFill("yellow")
    #=========================
    #draw nose (invariant)====
    nose_point=[Point(0,40),Point(-20,-20),Point(20,-20)]
    nose=Polygon(nose_point)
    nose.setFill("red")
    face.draw(win)
    nose.draw(win)
    #=========================
    #draw left_eye (variant with time)
    left_eye_control_points=[]
    left_eye_control_points=create_control_points_lefteye_runge(n_points,110,130,130,110,110,130)
    #left_eye_points=[Point(-170,110),Point(-120,130),Point(-70,110)]
    #draw_arc(left_eye_points[0],left_eye_points[1],left_eye_points[2],win)
    
    #==================================
    
    #draw right_eye (variant with time)
    right_eye_control_points=[]
    right_eye_control_points=create_control_points_righteye_runge(n_points,110,130,130,110,110,130)
    #right_eye_points=[Point(70,110),Point(120,130),Point(170,110)]
    #draw_arc(right_eye_points[0],right_eye_points[1],right_eye_points[2],win)
    #==================================
    
    #draw mouth (variant with time)
    mouth_control_points=[]
#    mouth_control_points.append([Point(-100,-100),Point(-50,-140),Point(0,-150),Point(50,-140),Point(100,-100)])
#    mouth_control_points.append([Point(-100,-120),Point(-50,-135),Point(0,-140),Point(50,-135),Point(100,-120)])
#    mouth_control_points.append([Point(-100,-125),Point(-50,-125),Point(0,-125),Point(50,-125),Point(100,-125)])
#    mouth_control_points.append([Point(-100,-140),Point(-50,-115),Point(0,-110),Point(50,-115),Point(100,-140)])
#    mouth_control_points.append([Point(-100,-150),Point(-50,-110),Point(0,-100),Point(50,-110),Point(100,-150)])
#    test_points=[Point(-100,-110),Point(-50,-110),Point(0,-110),Point(50,-110),Point(100,-110)]

    mouth_control_points=create_control_points_mouth_runge(n_points,-100,-150,-140,-110,-150,-100)
    print(len(mouth_control_points))
    #print(mouth_control_points[20])
    mouth_points=[]
    left_eye_points=[]
    right_eye_points=[]
    t=create_time_series(0,8,len(mouth_control_points))
    #   print(t)
  
    y_value=[]
    y_left_eye_value=[]
    y_right_eye_value=[]
    
    for i in range(len(mouth_control_points[0])):
        y=[]
        for j in range(len(mouth_control_points)):
            #print(mouth_control_points[j][0].getY())
            y.append(mouth_control_points[j][i].getY())
        y_value.append(y)
        #print(y_value[i])
        
    for i in range(len(left_eye_control_points[0])):
        a=[]
        b=[]
        for j in range(len(left_eye_control_points)):
            a.append(left_eye_control_points[j][i].getY())
            b.append(right_eye_control_points[j][i].getY())
        y_left_eye_value.append(a)
        y_right_eye_value.append(b)

#    print(len(t))
#    print(len(y_value))
#    print(len(y_value[0]))

    total_frames=640
    start_time=t[0]
    end_time=t[-1]
    dt=(end_time-start_time)/(total_frames-1)
    for i in range(total_frames):
        print(i*dt)
        aux_mouth_points=[]
        aux_lefteye_points=[]
        aux_righteye_points=[]
        for j in range(len(mouth_control_points[0])):
            #aux_points.append(Point(mouth_control_points[0][j].getX(),newton_interpolate(i*dt,t,y_value[j])))
            aux_mouth_points.append(Point(mouth_control_points[0][j].getX(),lagrange_interpolate(i*dt,t,y_value[j])))
        mouth_points.append(aux_mouth_points)
        
        for j in range(len(left_eye_control_points[0])):
            aux_lefteye_points.append(Point(left_eye_control_points[0][j].getX(),lagrange_interpolate(i*dt,t,y_left_eye_value[j])))
            aux_righteye_points.append(Point(right_eye_control_points[0][j].getX(),lagrange_interpolate(i*dt,t,y_right_eye_value[j])))
        left_eye_points.append(aux_lefteye_points)
        right_eye_points.append(aux_righteye_points)

        #print(mouth_points[i])

    win.autoflush=False
    x=[]
    y=[]

    for i in range(total_frames):
        x.append(i)
        y.append(mouth_points[i][0].getY())
    print(x)
    print(y)


    line=plt.plot(x,y,'o')

    plt.show()
    time.sleep(1.0)
    win.getMouse()
    for i in range(total_frames):
        arc1=draw_arc(mouth_points[i][0],mouth_points[i][1],mouth_points[i][2],win)
        arc2=draw_arc(mouth_points[i][2],mouth_points[i][3],mouth_points[i][4],win)
        arc3=draw_arc(left_eye_points[i][0],left_eye_points[i][1],left_eye_points[i][2],win)
        arc4=draw_arc(right_eye_points[i][0],right_eye_points[i][1],right_eye_points[i][2],win)
        print(i)
        update()
        time.sleep(dt)
        erase(arc1,win)
        erase(arc2,win)
        erase(arc3,win)
        erase(arc4,win)

#    arc1=draw_arc(mouth_control_points[0][0],mouth_control_points[0][1],mouth_control_points[0][2],win)
#    arc2=draw_arc(mouth_control_points[0][2],mouth_control_points[0][3],mouth_control_points[0][4],win)
#    time.sleep(0.2)
#    erase(arc1,win)
#    erase(arc2,win)
#    arc1=draw_arc(mouth_control_points[1][0],mouth_control_points[1][1],mouth_control_points[1][2],win)
#    arc2=draw_arc(mouth_control_points[1][2],mouth_control_points[1][3],mouth_control_points[1][4],win)
#    time.sleep(0.2)
#    erase(arc1,win)
#    erase(arc2,win)
#    arc1=draw_arc(mouth_control_points[2][0],mouth_control_points[2][1],mouth_control_points[2][2],win)
#    arc2=draw_arc(mouth_control_points[2][2],mouth_control_points[2][3],mouth_control_points[2][4],win)
#    time.sleep(0.2)
#    erase(arc1,win)
#    erase(arc2,win)
#    arc1=draw_arc(mouth_control_points[3][0],mouth_control_points[3][1],mouth_control_points[3][2],win)
#    arc2=draw_arc(mouth_control_points[3][2],mouth_control_points[3][3],mouth_control_points[3][4],win)
#    time.sleep(0.2)
#    erase(arc1,win)
#    erase(arc2,win)
#    arc1=draw_arc(mouth_control_points[4][0],mouth_control_points[4][1],mouth_control_points[4][2],win)
#    arc2=draw_arc(mouth_control_points[4][2],mouth_control_points[4][3],mouth_control_points[4][4],win)

#    arc1=draw_arc(test_points[0],test_points[1],test_points[2],win)
#    arc2=draw_arc(test_points[2],test_points[3],test_points[4],win)
    #==================================
#    pt=[]
#    pt.append(Point(20,20))
#    pt.append(Point(40,100))
#    pt.append(Point(150,200))
#    pt[0].draw(win)
#    pt[1].draw(win)
#    pt[2].draw(win)

#   arc=draw_arc(pt[0],pt[1],pt[2],win)
#    time.sleep(5.0)
#    erase(arc,win)
    print("OK")
    
def cubic_spline():
    n_points=int(input("Enter total number of control frames"))
    m=3
    win = GraphWin('Happy to sad', 300, 300)
    print(win)
    win.setCoords(-300,-300,300,300)
    
    #draw face (invariant)====
    face=Circle(Point(0,0),250)
    face.setFill("yellow")
    #=========================
    #draw nose (invariant)====
    nose_point=[Point(0,40),Point(-20,-20),Point(20,-20)]
    nose=Polygon(nose_point)
    nose.setFill("red")
    face.draw(win)
    nose.draw(win)
    #=========================
    #draw left_eye (variant with time)
    left_eye_control_points=[]
    left_eye_control_points=create_control_points_lefteye_runge(n_points,110,130,130,110,110,130)
    #left_eye_points=[Point(-170,110),Point(-120,130),Point(-70,110)]
    #draw_arc(left_eye_points[0],left_eye_points[1],left_eye_points[2],win)
    
    #==================================
    
    #draw right_eye (variant with time)
    right_eye_control_points=[]
    right_eye_control_points=create_control_points_righteye_runge(n_points,110,130,130,110,110,130)
    #right_eye_points=[Point(70,110),Point(120,130),Point(170,110)]
    #draw_arc(right_eye_points[0],right_eye_points[1],right_eye_points[2],win)
    #==================================
    
    #draw mouth (variant with time)
    mouth_control_points=[]
    
    mouth_control_points=create_control_points_mouth_runge(n_points,-100,-150,-140,-110,-150,-100)
    print(len(mouth_control_points))
    print(mouth_control_points[0])
    print(mouth_control_points[len(mouth_control_points)-1])
    mouth_points=[]
    left_eye_points=[]
    right_eye_points=[]
    t=create_time_series(0,8,len(mouth_control_points))
    #   print(t)
    
    y_value=[]
    y_left_eye_value=[]
    y_right_eye_value=[]
    
    for i in range(len(mouth_control_points[0])):
        y=[]
        for j in range(len(mouth_control_points)):
            #print(mouth_control_points[j][0].getY())
            y.append(mouth_control_points[j][i].getY())
        y_value.append(y)
    #print(y_value[i])
    
    for i in range(len(left_eye_control_points[0])):
        a=[]
        b=[]
        for j in range(len(left_eye_control_points)):
            a.append(left_eye_control_points[j][i].getY())
            b.append(right_eye_control_points[j][i].getY())
        y_left_eye_value.append(a)
        y_right_eye_value.append(b)
    
    #    print(len(t))
    #    print(len(y_value))
    #    print(len(y_value[0]))
    
    total_frames=640
    start_time=t[0]
    end_time=t[-1]
    dt=(end_time-start_time)/(total_frames-1)
    info_mouth=[]
    info_left_eye=[]
    info_right_eye=[]
    for j in range(len(mouth_control_points[0])):
        info_mouth.append(s_coefficient(m,t,y_value[j]))
    
    for j in range(len(left_eye_control_points[0])):
        info_left_eye.append(s_coefficient(m,t,y_left_eye_value[j]))
        info_right_eye.append(s_coefficient(m,t,y_right_eye_value[j]))

    print(info_mouth[0][0],info_mouth[1][0],info_mouth[2][0],info_mouth[3][0],info_mouth[4][0])
    print(info_left_eye[0][0],info_left_eye[1][0],info_left_eye[2][0])
    print(info_right_eye[0][0],info_right_eye[1][0],info_right_eye[2][0])

    for i in range(total_frames):
        print(i*dt)
        aux_mouth_points=[]
        aux_lefteye_points=[]
        aux_righteye_points=[]

        for j in range(len(mouth_control_points[0])):
            # aux_mouth_points.append(Point(mouth_control_points[0][j].getX(),lagrange_interpolate(i*dt,t,y_value[j])))
            #print(info_mouth[j][0])
            #print(info_mouth[j][1])
            aux_mouth_points.append(Point(mouth_control_points[0][j].getX(),s(i*dt,m,info_mouth[j][1],info_mouth[j][0])))
        mouth_points.append(aux_mouth_points)
        
        for j in range(len(left_eye_control_points[0])):
            aux_lefteye_points.append(Point(left_eye_control_points[0][j].getX(),s(i*dt,m,info_left_eye[j][1],info_left_eye[j][0])))
            aux_righteye_points.append(Point(right_eye_control_points[0][j].getX(),s(i*dt,m,info_right_eye[j][1],info_right_eye[j][0])))
        left_eye_points.append(aux_lefteye_points)
        right_eye_points.append(aux_righteye_points)
    
    #print(mouth_points[i])
    
    win.autoflush=False
    x=[]
    y=[]
    for i in range(total_frames):
        x.append(i*dt)
        y.append(mouth_points[i][0].getY())
    #print(x)
    #print(y)
    x_control=[]
    y_control=[]
    for i in range(len(mouth_control_points)):
        x_control.append(i/(n_points-1)*8)
        y_control.append(mouth_control_points[i][0].getY())
    line_1=plt.plot(x_control,y_control,'+')
    line=plt.plot(x,y,'-')
    plt.show()
    time.sleep(1.0)
    win.getMouse()
    for i in range(total_frames):
        arc1=draw_arc(mouth_points[i][0],mouth_points[i][1],mouth_points[i][2],win)
        arc2=draw_arc(mouth_points[i][2],mouth_points[i][3],mouth_points[i][4],win)
        arc3=draw_arc(left_eye_points[i][0],left_eye_points[i][1],left_eye_points[i][2],win)
        arc4=draw_arc(right_eye_points[i][0],right_eye_points[i][1],right_eye_points[i][2],win)
        print(i)
        update()
#        if (i==190):
#            win.getMouse()
        time.sleep(dt)
        erase(arc1,win)
        erase(arc2,win)
        erase(arc3,win)
        erase(arc4,win)

    print("OK")

    
    #line_1.draw(win)
    win.autoflush=True
    message = Text(Point(0, -180), 'Click to quit.')
    message.draw(win)
    win.getMouse()
    win.close()

def linear_spline():
    n_points=int(input("Enter total number of control frames"))
    m=1
    win = GraphWin('Happy to sad', 300, 300)
    print(win)
    win.setCoords(-300,-300,300,300)
    
    #draw face (invariant)====
    face=Circle(Point(0,0),250)
    face.setFill("yellow")
    #=========================
    #draw nose (invariant)====
    nose_point=[Point(0,40),Point(-20,-20),Point(20,-20)]
    nose=Polygon(nose_point)
    nose.setFill("red")
    face.draw(win)
    nose.draw(win)
    #=========================
    #draw left_eye (variant with time)
    left_eye_control_points=[]
    left_eye_control_points=create_control_points_lefteye_runge(n_points,110,130,130,110,110,130)
    #left_eye_points=[Point(-170,110),Point(-120,130),Point(-70,110)]
    #draw_arc(left_eye_points[0],left_eye_points[1],left_eye_points[2],win)
    
    #==================================
    
    #draw right_eye (variant with time)
    right_eye_control_points=[]
    right_eye_control_points=create_control_points_righteye_runge(n_points,110,130,130,110,110,130)
    #right_eye_points=[Point(70,110),Point(120,130),Point(170,110)]
    #draw_arc(right_eye_points[0],right_eye_points[1],right_eye_points[2],win)
    #==================================
    
    #draw mouth (variant with time)
    mouth_control_points=[]
    
    mouth_control_points=create_control_points_mouth_runge(n_points,-100,-150,-140,-110,-150,-100)
    print(len(mouth_control_points))
    print(mouth_control_points[0])
    print(mouth_control_points[len(mouth_control_points)-1])
    mouth_points=[]
    left_eye_points=[]
    right_eye_points=[]
    t=create_time_series(0,8,len(mouth_control_points))
    #   print(t)
    
    y_value=[]
    y_left_eye_value=[]
    y_right_eye_value=[]
    
    for i in range(len(mouth_control_points[0])):
        y=[]
        for j in range(len(mouth_control_points)):
            #print(mouth_control_points[j][0].getY())
            y.append(mouth_control_points[j][i].getY())
        y_value.append(y)
    #print(y_value[i])
    
    for i in range(len(left_eye_control_points[0])):
        a=[]
        b=[]
        for j in range(len(left_eye_control_points)):
            a.append(left_eye_control_points[j][i].getY())
            b.append(right_eye_control_points[j][i].getY())
        y_left_eye_value.append(a)
        y_right_eye_value.append(b)
    
    #    print(len(t))
    #    print(len(y_value))
    #    print(len(y_value[0]))
    
    total_frames=640
    start_time=t[0]
    end_time=t[-1]
    dt=(end_time-start_time)/(total_frames-1)
    info_mouth=[]
    info_left_eye=[]
    info_right_eye=[]
    for j in range(len(mouth_control_points[0])):
        info_mouth.append(s_coefficient_linear(m,t,y_value[j]))
    
    for j in range(len(left_eye_control_points[0])):
        info_left_eye.append(s_coefficient_linear(m,t,y_left_eye_value[j]))
        info_right_eye.append(s_coefficient_linear(m,t,y_right_eye_value[j]))

    print(info_mouth[0][0],info_mouth[1][0],info_mouth[2][0],info_mouth[3][0],info_mouth[4][0])
    print(info_left_eye[0][0],info_left_eye[1][0],info_left_eye[2][0])
    print(info_right_eye[0][0],info_right_eye[1][0],info_right_eye[2][0])

    for i in range(total_frames):
        print(i*dt)
        aux_mouth_points=[]
        aux_lefteye_points=[]
        aux_righteye_points=[]
    
        for j in range(len(mouth_control_points[0])):
        # aux_mouth_points.append(Point(mouth_control_points[0][j].getX(),lagrange_interpolate(i*dt,t,y_value[j])))
        #print(info_mouth[j][0])
        #print(info_mouth[j][1])
            aux_mouth_points.append(Point(mouth_control_points[0][j].getX(),s_linear(i*dt,m,info_mouth[j][1],info_mouth[j][0])))
        mouth_points.append(aux_mouth_points)
        
        for j in range(len(left_eye_control_points[0])):
            aux_lefteye_points.append(Point(left_eye_control_points[0][j].getX(),s_linear(i*dt,m,info_left_eye[j][1],info_left_eye[j][0])))
            aux_righteye_points.append(Point(right_eye_control_points[0][j].getX(),s_linear(i*dt,m,info_right_eye[j][1],info_right_eye[j][0])))
        left_eye_points.append(aux_lefteye_points)
        right_eye_points.append(aux_righteye_points)

    #print(mouth_points[i])

    win.autoflush=False
    x=[]
    y=[]
    for i in range(total_frames):
        x.append(i*dt)
        y.append(mouth_points[i][0].getY())
    #print(x)
    #print(y)
    x_control=[]
    y_control=[]
    for i in range(len(mouth_control_points)):
        x_control.append(i/(n_points-1)*8)
        y_control.append(mouth_control_points[i][0].getY())
    line_1=plt.plot(x_control,y_control,'+')
    line=plt.plot(x,y,'-')
    plt.show()
    time.sleep(1.0)
        
    win.getMouse()

    for i in range(total_frames):
        arc1=draw_arc(mouth_points[i][0],mouth_points[i][1],mouth_points[i][2],win)
        arc2=draw_arc(mouth_points[i][2],mouth_points[i][3],mouth_points[i][4],win)
        arc3=draw_arc(left_eye_points[i][0],left_eye_points[i][1],left_eye_points[i][2],win)
        arc4=draw_arc(right_eye_points[i][0],right_eye_points[i][1],right_eye_points[i][2],win)
        print(i)
        update()
        time.sleep(dt)
        erase(arc1,win)
        erase(arc2,win)
        erase(arc3,win)
        erase(arc4,win)

    print("OK")


    #line_1.draw(win)
    win.autoflush=True
    message = Text(Point(0, -180), 'Click to quit.')
    message.draw(win)
    win.getMouse()
    win.close()

print("Choose a interpolation method")
print("1 for polynomial interpolation")
print("2 for cubic spline interpolation")
print("3 for linear spline interpolation")
option=int(input("Enter a number"))
if option==1:
    polynomial()
elif option==2:
    cubic_spline()
else:
    linear_spline()
