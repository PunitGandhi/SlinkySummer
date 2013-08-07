from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
def main():
    saveFile = 'heavytop.csv'
    #M02082013/slinkym200l027_5r1ct.csv'
    unLength = .075
    length = unLength*20.
    totMass = .200
    nMasses = 100
    collisions = True 
    g = 9.8
    dt = .0002
    nIter = 2000
    plotAll = True 
    #Calculated parameters
    lPer = unLength/(nMasses-1)
    m = totMass/nMasses
    k = (nMasses-1)*nMasses*m*g/(2*(length-unLength))
    masses = np.array([m for ii in xrange(nMasses)])
    masses[0] = 100*masses[0]
    collapsed = [False for ii in xrange(nMasses)]
    y = np.zeros((nIter+1,nMasses))
    v = np.zeros_like(y)
    t = np.array([dt*ii for ii in xrange(nIter+1)])
    topInit = 0
    temp = topInit
    #Set intitial positions
    for ii in xrange(nMasses):
        y[0,ii] = -temp
        temp += lPer+(nMasses-ii-1)*m*g/k
    #Loop over timesteps and elements of the slinky
    for ii in xrange(nIter):
        for jj in xrange(nMasses-1,-1,-1):
            #If the mass has not yet hit, then update pos, vel w/ Newton's Law
            if not collapsed[jj]:
                if jj == 0:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(1. if np.abs(y[ii,jj+1]-y[ii,jj])>lPer else 0.)*(y[ii,jj+1]-y[ii,jj]-np.sign(y[ii,jj+1]-y[ii,jj])*lPer)*k)/masses[jj]*dt
                elif jj == nMasses-1:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(1. if np.abs(y[ii,jj-1]-y[ii,jj])>lPer else 0.)*(y[ii,jj-1]-y[ii,jj]-np.sign(y[ii,jj-1]-y[ii,jj])*lPer)*k)/masses[jj]*dt
                else:
                    v[ii+1,jj] = v[ii,jj]+(-masses[jj]*g+(1. if np.abs(y[ii,jj-1]-y[ii,jj])>lPer else 0.)*(y[ii,jj-1]-y[ii,jj]-np.sign(y[ii,jj-1]-y[ii,jj])*lPer)*k+
                                           (1. if np.abs(y[ii,jj+1]-y[ii,jj])>lPer else 0.)*(y[ii,jj+1]-y[ii,jj]-np.sign(y[ii,jj+1]-y[ii,jj])*lPer)*k)/masses[jj]*dt
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
    plt.title(saveFile[10:-4])
    axv = plt.subplot2grid((2,2),(1,0),colspan=2)
    plt.plot(*vplots)
    plt.ylabel('Velocity')
    plt.xlabel('Time')
    #plt.show()
    plt.savefig(saveFile[:-3]+'png')

if __name__ == '__main__':
    main()
