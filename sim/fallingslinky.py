from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
def main():
    saveFile = 'slinkySim.csv'
    length = 1.3
    totMass = .21
    nMasses = 4
    collisions = True 
    g = 9.8
    dt = .0002
    nIter = 2000
    plotAll = False
    #Calculated parameters
    m = totMass/nMasses
    k = (nMasses-1)*nMasses*m*g/(2*length)
    masses = np.array([m for ii in xrange(nMasses)])
    collapsed = [False for ii in xrange(nMasses)]
    y = np.zeros((nIter+1,nMasses))
    v = np.zeros_like(y)
    t = np.array([dt*ii for ii in xrange(nIter+1)])
    topInit = 0
    temp = topInit
    #Set intitial positions
    for ii in xrange(nMasses):
        y[0,ii] = -temp
        temp += (nMasses-ii-1)*m*g/k
    #Loop over timesteps and elements of the slinky
    for ii in xrange(nIter):
        for jj in xrange(nMasses-1,-1,-1):
            #If the mass has not yet hit, then update pos, vel w/ Newton's Law
            if not collapsed[jj]:
                if jj == 0:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(y[ii,jj+1]-y[ii,jj])*k)/masses[jj]*dt
                elif jj == nMasses-1:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(y[ii,jj-1]-y[ii,jj])*k)/masses[jj]*dt
                else:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(y[ii,jj-1]+y[ii,jj+1]-2*y[ii,jj])*k)/masses[jj]*dt
                y[ii+1,jj] = y[ii,jj]+v[ii+1,jj]*dt
            #If the mass has hit, then keep it stuck to lower mass
            else:
                v[ii+1,jj] = v[ii+1,jj+1]
                y[ii+1,jj] = y[ii+1,jj+1]
            #Check to see if it has hit lower mass, if so, do inelastic collision
            if collisions == True and not collapsed[jj] and jj != nMasses-1 and y[ii+1,jj] < y[ii+1,jj+1]:
                y[ii+1,jj] = y[ii+1,jj+1]
                v[ii+1,jj] = (masses[jj]*v[ii+1,jj]+masses[jj+1]*v[ii+1,jj+1])/(masses[jj]+masses[jj+1])
                v[ii+1,jj+1] = v[ii+1,jj]
                masses[jj+1] = masses[jj]+masses[jj+1]
                masses[jj] = 0
                collapsed[jj] = True


    if plotAll:
        yplots = ()
        vplots = ()
        for ii in xrange(nMasses):
            yplots += (t,y[:,ii])
            vplots += (t,v[:,ii])
    else:
        yplots = (t,y[:,0],t,y[:,nMasses-1])
        vplots = (t,v[:,0],t,v[:,nMasses-1])
    print 'Saving file: '+saveFile
    output = np.append(y,v,axis=1)
    with open(saveFile,'w+b') as f:
       np.savetxt(f,np.concatenate((np.array([t]).T,y,v),axis=1),delimiter=',') 

    axy = plt.subplot2grid((2,2),(0,0),colspan=2)
    axy.plot(*yplots)
    plt.ylabel('Height')
    axv = plt.subplot2grid((2,2),(1,0),colspan=2)
    plt.plot(*vplots)
    plt.ylabel('Velocity')
    plt.xlabel('Time')
    plt.show()

if __name__ == '__main__':
    main()
