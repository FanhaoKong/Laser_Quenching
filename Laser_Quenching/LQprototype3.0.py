import numpy as np
import matplotlib.pyplot as plt

xlength=100;ylength=100;zlength=50;unit=0.5;
speedx=10.0;speedy=0.0;laserpower=30000;
lightradius=15;#radius:格数
density=8.0;heatcap=15.0;initialtemp=300;starttemp=301;conb=33;conk=-conb/2000;
h=10;
timestep=0.1;rounds=500;
startposition=np.array([xlength/5.0,ylength/2.0],dtype='float64')

def conduct(temp):
    #热导率与温度的关系
    return conb+conk*temp

def initialfield():
    #温度场初始化
    tempfield=np.full((xlength,ylength,zlength),starttemp,dtype='float64')
    return tempfield

def boundary():
    #边界条件，独立封装
    #boundary最后维度索引：
    #0：x0表面；1：x1表面；2：y0表面；3：y1表面；4：z0表面；5：z1表面（0：靠近负轴向；1：靠近正轴向）
    boundary=np.zeros((xlength,ylength,zlength,6),dtype='int32')
    for i in range(xlength):
        for j in range(ylength):
            for k in range(zlength):
                if i==0:
                    boundary[i,j,k,0]=1
                if i==xlength-1:
                    boundary[i,j,k,1]=1
                if j==0:
                    boundary[i,j,k,2]=1
                if j==ylength-1:
                    boundary[i,j,k,3]=1
                if k==0:
                    boundary[i,j,k,4]=1
                if k==zlength-1:
                    boundary[i,j,k,5]=1
    return boundary

def nextlight(position):
    #确定下一时刻激光位置，激光光照点作匀速直线运动
    return position+timestep*np.array([speedx,speedy])
#boundary最后维度索引：
#0：x0表面；1：x1表面；2：y0表面；3：y1表面；4：z0表面；5：z1表面（0：靠近负轴向；1：靠近正轴向）

def nextfield(currentfield,boundary,position):
    print(np.max(currentfield),end='   ')
    print(np.min(currentfield))
    for i in range(xlength):
        for j in range(ylength):
            for k in range(zlength):
                heatconduct=conduct(currentfield[i,j,k])
                heatflow=np.zeros(6,dtype='float64')
                sunit=unit**2;vunit=unit**3;
                if boundary[i,j,k,0]==1:
                    heatflow[0]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[0]=-heatconduct*(currentfield[i,j,k]-currentfield[i-1,j,k])/unit
                if boundary[i,j,k,1]==1:
                    heatflow[1]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[1]=-heatconduct*(currentfield[i+1,j,k]-currentfield[i,j,k])/unit
                if boundary[i,j,k,2]==1:
                    heatflow[2]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[2]=-heatconduct*(currentfield[i,j,k]-currentfield[i,j-1,k])/unit
                if boundary[i,j,k,3]==1:
                    heatflow[3]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[3]=-heatconduct*(currentfield[i,j+1,k]-currentfield[i,j,k])/unit
                if boundary[i,j,k,4]==1:
                    heatflow[4]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[4]=-heatconduct*(currentfield[i,j,k]-currentfield[i,j,k-1])/unit
                if boundary[i,j,k,5]==1:
                    if np.sqrt((i+0.5-position[0])**2+(j+0.5-position[1])**2)<lightradius:
                        heatflow[5]=laserpower
                    
                    elif np.sqrt((i+0.5-position[0])**2+(j+0.5-position[1])**2)<lightradius+1:
                        aaa=np.random.rand()
                        if aaa<=0.5:
                            heatflow[5]=laserpower
                        else:
                            heatflow[5]=h*(initialtemp-currentfield[i,j,k])
                    else:
                        heatflow[5]=h*(initialtemp-currentfield[i,j,k])
                else:
                    heatflow[5]=-heatconduct*(currentfield[i,j,k+1]-currentfield[i,j,k])/unit
                heat=np.sum(heatflow)*sunit
                currentfield[i,j,k]=currentfield[i,j,k]+timestep*(heat/(density*heatcap*vunit))
                if currentfield[i,j,k]<initialtemp:
                    currentfield[i,j,k]=initialtemp
    print(currentfield[xlength-1,int(ylength/2),int(zlength/2)],end='   ')
    print(currentfield[int(xlength/2),ylength-1,int(zlength/2)],end='   ')
    print(currentfield[int(xlength/2),int(ylength/2),0])
    position=nextlight(position)
    return currentfield,position

'''
def display(currentfield,step):
    z=np.zeros((xlength,ylength),dtype='float64')
    depth=np.zeros((xlength,zlength),dtype='float64')
    for i in range(xlength):
        for j in range(ylength):
            z[i,j]=currentfield[i,j,zlength-1]
    for i in range(xlength):
        for k in range(zlength):
            depth[i,k]=currentfield[i,int(ylength/2.0),k]
    depth=depth.T
    x,y=np.ogrid[0:100:100j,0:100:100j]
    x2,z2=np.ogrid[0:100:100j,0:50:50j]
    extent=[np.min(x),np.max(x),np.min(y),np.max(y)]
    extent2=[np.min(x2),np.max(x2),np.min(z2),np.max(z2)]
    plt.figure(figsize=(12,10))
    plt.subplot(221)
    cs=plt.contour(z,8,extent=extent)
    plt.clabel(cs)
    plt.subplot(222)
    plt.contourf(x.reshape(-1),y.reshape(-1),z,6,cmap='autumn')
    plt.colorbar()
    plt.subplot(223)
    cs2=plt.contour(depth,8,extent=extent2)
    plt.clabel(cs2)
    plt.subplot(224)
    plt.contourf(x2.reshape(-1),z2.reshape(-1),depth,6,cmap='autumn')
    plt.colorbar()
    plt.savefig("F:/LQresult/LQresult-{}.jpg".format(step),dpi=600)
    plt.show()
'''
def display(currentfield,step):
    z=np.zeros((xlength,ylength),dtype='float64')
    depth=np.zeros((xlength,zlength),dtype='float64')
    for i in range(xlength):
        for j in range(ylength):
            z[i,j]=currentfield[i,j,zlength-1]
    for i in range(xlength):
        for k in range(zlength):
            depth[i,k]=currentfield[i,int(ylength/2.0),k]
    depth=depth.T
    x,y=np.ogrid[0:100:100j,0:100:100j]
    x2,z2=np.ogrid[0:100:100j,0:50:50j]
    extent=[np.min(x),np.max(x),np.min(y),np.max(y)]
    extent2=[np.min(x2),np.max(x2),np.min(z2),np.max(z2)]
    plt.figure(figsize=(15,6))
    plt.subplot(121)
    plt.contourf(x.reshape(-1),y.reshape(-1),z,50,cmap='rainbow')
    plt.colorbar()
    plt.subplot(122)
    plt.contourf(x2.reshape(-1),z2.reshape(-1),depth,50,cmap='rainbow')
    plt.colorbar()
    plt.savefig("F:/LQresult/LQresult-{}.jpg".format(step),dpi=600)
    plt.show()

def main():
    currentfield=initialfield()
    bouncondition=boundary()
    position=startposition
    for m1 in range(rounds):
        print(m1,end='   ')
        currentfield,position=nextfield(currentfield,bouncondition,position)
        #if m1%10==0:
        #    display(currentfield,m1)
        if np.abs(position[0]-xlength)<lightradius*1.2 or \
            np.abs(position[1]-ylength)<lightradius*1.2 or \
            np.abs(position[0])<lightradius*1.2 or np.abs(position[1])<lightradius*1.2:
            break
        display(currentfield,m1)
    display(currentfield,m1)

main()

                    