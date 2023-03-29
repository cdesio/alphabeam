import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import random

def calculateEnergy(stepsInTrack,distanceTravelled):
    distances = [findDistanceBetween((stepsInTrack[i][0],stepsInTrack[i][1],stepsInTrack[i][2]), (stepsInTrack[0][0],stepsInTrack[0][1],stepsInTrack[0][2])) for i in range(len(stepsInTrack))]
    energies = [a[6] for a in stepsInTrack]

    for ind, value in enumerate(distances):
        if distanceTravelled > value:
            lower = ind
        if  distanceTravelled < value:
            upper = ind
            break

    energy = energies[lower]-(energies[lower]-energies[upper])*(distanceTravelled-distances[lower])/(distances[upper]-distances[lower])
    return energy

def checkCross(point, momentum, stepsInTrack, distanceTravelled):
    tXneg = (-.00017 - point[0]) / momentum[0]
    tXpos = (.00017 - point[0]) / momentum[0] 

    tYneg = (-.00017 - value[1]) / momentum[1] 
    tYpos = (.00017 - point[1]) / momentum[1] 

    tZneg = (-.00017 - point[2]) / momentum[2]
    tZpos = (.00017 - point[2]) / momentum[2] 

    posValues = [a for a in [tXneg,tXpos,tYpos,tZneg,tZpos] if a>0]
    if plotGraphs:
        ax3.quiver(point[0], point[1], point [2], momentum[0], momentum[1], momentum[2], length=.0001, normalize=True)

    if min(posValues) == tXneg:
        # print("-ve X")
        t = tXneg
        newPos = [t*a + o for a,o in zip(momentum,point)]
        newPos [0] = +.00017

    elif min(posValues) == tXpos:
        # print("+ve X")
        t = tXpos
        newPos = [t*a + o for a,o in zip(momentum,point)]
        newPos [0] = -.00017

    elif min(posValues) == tYneg:
        # print("-ve Y")
        t = tYpos
        newPos = [t*a + o for a,o in zip(momentum,point)]
        newPos [0] = -.00017

    elif min(posValues) == tYpos:
        # print("+ve Y")
        return False, False # leaves ring

    elif min(posValues) == tZneg:
        # print("-ve Z")
        t = tZneg
        newPos = [t*a + o for a,o in zip(momentum,point)]
        newPos [2] = +.00017

    elif min(posValues) == tZpos:
        # print("+ve Z")
        t = tZpos
        newPos = [t*a + o for a,o in zip(momentum,point)]
        newPos [2] = -.00017

    elif len(posValues)==0:
        print("no possible values")

    cross = [a*t for a in momentum]
    d = findDistanceBetween(cross, (0,0,0)) #distance to cross point

    # do particle steps exit box
    if type(stepsInTrack[0])==np.float64:# only 1 step in track, does not exit
        maxStepDistFromStart = 0 
    else:
        maxStepDistFromStart = np.max([findDistanceBetween((stepsInTrack[i][0],stepsInTrack[i][1],stepsInTrack[i][2]), (stepsInTrack[0][0],stepsInTrack[0][1],stepsInTrack[0][2])) for i in range(len(stepsInTrack))])


    if d+distanceTravelled<maxStepDistFromStart and abs(newPos[0])<=0.00017 and abs(newPos[1])<=0.00017 and abs(newPos[2])<=0.00017:
        # crosses into adjacent box
        # print("need to add", newPos, particleMap[stepsInTrack[0][8]])
        if plotGraphs:
            ax3.quiver(newPos[0], newPos[1], newPos [2], momentum[0], momentum[1], momentum[2], length=.0001, normalize=True, color = 'r')
    else:
        return False, False
        # check if crosses again

    if abs(newPos[0])>0.00017 or abs(newPos[1])>0.00017 or abs(newPos[2])>0.00017 :
        print()
    
    return newPos, distanceTravelled+d

def plot_cylinder_along_z(center_x,center_y,radius,height_z, ax):
    z = np.linspace(-height_z, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y

    ax.plot_surface(x_grid, y_grid, z_grid, alpha=0.5)


def get_cube(a,b,c):   
    phi = np.arange(1,10,2)*np.pi/4
    Phi, Theta = np.meshgrid(phi, phi)

    x = np.cos(Phi)*np.sin(Theta)
    y = np.sin(Phi)*np.sin(Theta)
    z = np.cos(Theta)/np.sqrt(2)
    return x,y,z

def plotAll(ax, data, mmTonm):
    x = data[:,0]*mmTonm
    y = data[:,1]*mmTonm
    z = data[:,2]*mmTonm

    particle = [int(data[i,8]) for i in range(data.shape[0])]
    trackID = [int(data[i,12]) for i in range(data.shape[0])]

    cmap = mpl.cm.nipy_spectral
    norm = mpl.colors.Normalize()
    norm.autoscale(np.linspace(1,10,10))


    numToplot = len(trackID)
    ax.scatter(x[:numToplot],y[:numToplot],z[:numToplot], '.', c = cmap(norm(particle[:numToplot])))

    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")

    ax.set_xlim(-2*mmTonm, 2*mmTonm)
    ax.set_ylim(-2*mmTonm, 2*mmTonm)
    ax.set_zlim(-3*mmTonm, 3*mmTonm)

    # numToplot = len(x)
    # trackIDToPlot = 19
    # ax = plt.figure().add_subplot(projection='3d')
    # # plot momentum direction for each particle
    # for i in range(1,11):
    #     # to plot different colours need to plot each particle seperately
    #     xt =   [x[a] for a in range(numToplot) if (particle[a]==i)]
    #     yt =   [y[a] for a in range(numToplot) if (particle[a]==i) ]
    #     zt =   [z[a] for a in range(numToplot) if (particle[a]==i) ]
    #     ut =   [u[a] for a in range(numToplot) if (particle[a]==i) ]
    #     vt =   [v[a] for a in range(numToplot) if (particle[a]==i) ]
    #     wt =   [w[a] for a in range(numToplot) if (particle[a]==i) ]
    #     ax.quiver(xt,yt,zt,ut,vt,wt, length=100, normalize=True, color = cmap(norm(i)))



    # ax.set_xlabel("x (nm)")
    # ax.set_ylabel("y (nm)")
    # ax.set_zlabel("z (nm)")

    legend_elements = [
                        Line2D([0], [0], marker='o', color = cmap(norm(1)), label='e-', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(2)), label='gamma', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(3)), label='alpha', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(4)), label='Rn220', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(5)), label='Po216', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(6)), label='Pb212', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(7)), label='Bi212', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(8)), label='Tl208', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(9)), label='Po212', markersize=6, linestyle="None"),
                        Line2D([0], [0], marker='o', color = cmap(norm(10)), label='Pb208', markersize=6, linestyle="None"),
    ]

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(-0.1, .85, 1.3, 0.0), 
        loc="lower left",
        mode="expand", 
        ncol=5,
        frameon=False
    )
    # plt.show() 

def findDistanceBetween(a,b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )**0.5

def updateMomentum(i):
    theta = np.arcsin(i[0]/(i[1]**2 + i[0]**2)**0.5)
    if i[0]>0 and i[1]>0:
        # positive-positive quadrant
        theta = theta
    elif i[0]>0 and i[1]<0:
        # positive-negative quadrant
        theta = np.pi - theta
    elif i[0]<0 and i[1]<0:
        # negative-negative quadrant
        theta = abs(theta)+np.pi
    elif i[0]<0 and i[1]>0:
        # negative-positive quadrant
        theta = theta
    normMomentum = [i[3],i[4],i[5]]/(i[3]**2 + i[4]**2 + i[5]**2)**0.5

    # newMomentum_ = [normMomentum[0]*np.cos(theta)-normMomentum[1]*np.sin(theta), normMomentum[0]*np.sin(theta)+normMomentum[1]*np.cos(theta), -1*normMomentum[2]]

    # ax.quiver(i[0]*mmTonm, i[1]*mmTonm, i[2]*mmTonm, i[3]*mmTonm, i[4]*mmTonm, i[5]*mmTonm, length=100000, normalize=True, color = 'b')

    # ax.quiver(i[0]*mmTonm*np.cos(theta)-i[1]*mmTonm*np.sin(theta), i[0]*mmTonm*np.sin(theta)+i[1]*mmTonm*np.cos(theta), 0, normMomentum[0]*np.cos(theta)-normMomentum[1]*np.sin(theta), normMomentum[0]*np.sin(theta)+normMomentum[1]*np.cos(theta), -1*normMomentum[2], length=100000, normalize=True, color = 'r')
    # print(newMomentum_[1])
    newMomentum_ = [normMomentum[0]*np.cos(theta)-normMomentum[1]*np.sin(theta), normMomentum[0]*np.sin(theta)+normMomentum[1]*np.cos(theta), normMomentum[2]]
    if plotGraphs:
        ax.quiver(i[0]*mmTonm, i[1]*mmTonm, i[2]*mmTonm, i[3]*mmTonm, i[4]*mmTonm, i[5]*mmTonm, length=100000, normalize=True, color = 'b')

        ax.quiver(i[0]*mmTonm*np.cos(theta)-i[1]*mmTonm*np.sin(theta), i[0]*mmTonm*np.sin(theta)+i[1]*mmTonm*np.cos(theta), 0, newMomentum_[0], newMomentum_[1], newMomentum_[2], length=100000, normalize=True, color = 'r')

    return newMomentum_
    
particleMap = {1: "e-",
               2: "gamma",
               3: "alpha",
               4: "Rn220",
               5: "Po216",
               6: "Pb212",
               7: "Bi212",
               8: "Tl208",
               9: "Po212",
               10: "Pb208"}

mmTonm = 1000*1000
try:
    with open("compareToBoxes.bin", "rb") as f:
        numpy_data = np.fromfile(f,np.float32)
except IOError:
    print('Error While Opening the file!')  

# boxCentre = 0.155 #mm
boxCentre = 0.235 #mm
boxCopyID = 8
boxMin = boxCentre - 0.00015
boxMax = boxCentre + 0.00015
data = numpy_data.reshape(int(len(numpy_data)/14),14)
plotGraphs = True

if plotGraphs:
    ax = plt.figure().add_subplot(projection='3d')
    plot_cylinder_along_z(0,0,0.15*mmTonm,3*mmTonm, ax)
    # plot_cylinder_along_z(0,0,0.155*mmTonm,0.1*mmTonm, ax)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")


eventIDs = list(set([int(data[i,7]) for i in range(data.shape[0])]))
tracksToConsider = []
for event in eventIDs:
    dataEvent = data[data[:,7]==event] 
    trackIDs = list(set([int(dataEvent[i,12]) for i in range(dataEvent.shape[0])]))

    for i in range(len(trackIDs)):
        d = dataEvent[dataEvent[:,12]==trackIDs[i]] #select data with given track ID
        d = d[d[:,9]==boxCopyID] #select data with given box copy number
        x = d[:,0]
        y = d[:,1]

        r = (x**2 + y**2)**0.5

        if len(r) ==0:
            continue


        a = []
        if np.isclose(boxMin,r[0].round(5)) or np.isclose(boxMax, r[0].round(5)): # track is entering volume not created in volume - need to do decay products differently

            for ind,value in enumerate(r):
                # if np.isclose(0.15485, value) or np.isclose(0.15515, value):
                if (value.round(5)>boxMin and value.round(5)<boxMax) or np.isclose(boxMin,value.round(5)) or np.isclose(boxMax, value.round(5)):
                    a.append(True)

                else:
                    a.append(False)

            # a = np.where((0.15485==r)|(r==0.15515), True, False)
            if sum(a)>0:
                # plotAll(ax, d[a], mmTonm)
                tracksToConsider.append(d[a])
                # print(r)

        else:
            if r[0]>boxMin and r[0] < boxMax and d[0][13]==1:
                print("track starts inside from radioactive decay - keep", particleMap[d[0][8]], "event = ", d[0][7], "creator process ", d[0][13], "track = ", d[0][12], "r=", r[0])
                tracksToConsider.append(d[0])



newMomentum = []
newAllOther = []
newPos = []
count = 0
for i in tracksToConsider:
    # print(len(i), "particle = ", particleMap[i[0][8]])
    # ax.scatter(i[0][0]*mmTonm, i[0][1]*mmTonm, i[0][2]*mmTonm, color ='b')
    if type(i[0])==np.float64:# only 1 step in track
        print()
        print("starts inside, particle = ", particleMap[i[8]], "event = ", i[7], "creator process ", i[13], "track = ", i[12])

        newMomentum_ = updateMomentum(i)
        r = (i[1]**2 + i[0]**2)**0.5

        newPos.append([random.uniform(-.00017,.00017),r-boxCentre,random.uniform(-.00017,.00017)])      

        newMomentum.append(newMomentum_)
        newAllOther.append(i[6:12])


    else:
        # print("theta = ", theta)
        # ax.scatter(i[0][0]*mmTonm*np.cos(theta)-i[0][1]*mmTonm*np.sin(theta), i[0][0]*mmTonm*np.sin(theta)+i[0][1]*mmTonm*np.cos(theta), i[0][2]*mmTonm, color ='r')


        r = (i[0][1]**2 + i[0][0]**2)**0.5

        newMomentum_ = updateMomentum(i[0])

        if np.isclose(r.round(5),boxMin):
            if newMomentum_[1]>0:
                # check is coming into box
                newPos.append([random.uniform(-.00017,.00017),-.00017, random.uniform(-.00017,.00017)])
            else:
                continue
        elif np.isclose(r.round(5),boxMax):
            if newMomentum_[1]<0:
                # check is coming into box
                newPos.append([random.uniform(-.00017,.00017),.00017,random.uniform(-.00017,.00017)])      
            else:
                continue
        
        else:
            print("????")

        newMomentum.append(newMomentum_)
        newAllOther.append(i[0][6:12])


    
    # print(f"/gps/direction {newMomentum[count][0]} {newMomentum[count][1]} {newMomentum[count][2]}")

    count +=1

if plotGraphs:
    cmap = mpl.cm.nipy_spectral
    norm = mpl.colors.Normalize()
    norm.autoscale(np.linspace(1,10,10))

# arrow colours are all arrows, then arrow 1 head 1, arrow 1 head2, arrow 2 head 1 ...
    colourMap = [x[2] for x in newAllOther]
    for j, _ in enumerate(newPos):
        colourMap.append(colourMap[j])
        colourMap.append(colourMap[j])

    colourMap = cmap(norm(colourMap))
    ax2 = plt.figure().add_subplot(projection='3d')
    ax2.quiver([x[0] for x in newPos],[x[1] for x in newPos],[x[2] for x in newPos],[x[0] for x in newMomentum],[x[1] for x in newMomentum],[x[2] for x in newMomentum], length=0.0001, normalize=True, color = colourMap)

    a = .00017*2
    b = .00017*2
    c = .00017*2
    x,y,z = get_cube(a,b,c)

    ax2.plot_surface(x*a, y*b, z*c, color='b', alpha=0.1)

    ax2.set_xlim(-.0002, .0002)
    ax2.set_ylim(-.0002, .0002)
    ax2.set_zlim(-.0002, .0002)
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("y (mm)")
    ax2.set_zlabel("z (mm)")

    ax3 = plt.figure().add_subplot(projection='3d')
    a = .00017*2
    b = .00017*2
    c = .00017*2
    x,y,z = get_cube(a,b,c)

    ax3.plot_surface(x*a, y*b, z*c, color='b', alpha=0.1)

    ax3.set_xlim(-.0002, .0002)
    ax3.set_ylim(-.0002, .0002)
    ax3.set_zlim(-.0002, .0002)
    ax3.set_xlabel("x (mm)")
    ax3.set_ylabel("y (mm)")
    ax3.set_zlabel("z (mm)")

toAdd = []
newEnergy = []
newOther = []
for ind, value in enumerate(newPos):
    # check where/if track exits cube
    point = value
    momentum = newMomentum[ind]
    add = True
    distanceTravelled = 0
    while add:
        add, distanceTravelled = checkCross(point, momentum, tracksToConsider[ind],distanceTravelled)
        if add:
            # workout energy by interpolating
            toAdd.append(add+ momentum + [calculateEnergy(tracksToConsider[ind],distanceTravelled)] + list(tracksToConsider[ind][0][7:12] ))

            point = add
            

newBin = np.hstack([np.array(newPos), np.array(newMomentum), np.array(newAllOther)])

newBin = np.vstack([newBin, np.array(toAdd)])
print("number of particles in box = ", len(np.array(newBin)))

newBin = newBin.flatten()

newBin.astype('float64').tofile(f"compareToBoxesTransformed_{int(boxCentre*1000)}um.bin")


# data = numpy_data.reshape(int(len(numpy_data)/12),12)

plt.show()