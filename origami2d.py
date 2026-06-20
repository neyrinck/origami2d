import numpy as np

def foldfinder(x,boxlength):
    """ finds out-of-order particles. Does it in two parts to straddle boundary """
    ng = x.shape[0]
    n2 = ng//2
    n14 = ng//4
    n34 = 3*n14
    narange = np.arange(ng)
    madder = np.zeros(ng,dtype='bool')
    madder[n14:n34] = (np.argsort(x) != narange)[n14:n34]
    xhalf = np.concatenate((x[n2:],x[:n2]+boxlength))
    halforder = (np.argsort(xhalf) != narange)
    madder[n34:] = halforder[n14:n2]
    madder[:n14] = halforder[n2:n34]

    return madder

def m2D_cart(x,y,boxlength):
    """
    Finds out-of-ortder particles along the two Cartesian axes.
    Assumes 2D arrays of x and y, such that x advances with the 1st index of the array,
    and y advances with the 2nd.
    Also assumes that particles are _not_ wrapped around a periodic boundary -- may need to unwrap them
    """
    
    N = x.shape[0] # assumes shape[0] == shape[1]
    if (x.shape[0] != x.shape[1]) | (x.shape != y.shape):
        print('warning! non-square array! we assume square')
        
    m = np.zeros((N,N),dtype='int8') # the cartesian directions

    # first along the x and y axes
    for i in range(N):
        m[:,i] += foldfinder(x[:,i],boxlength)
        m[i,:] += foldfinder(y[i,:],boxlength)
    return m    

def m2D_45deg(x,y,boxlength):
    #diagonals
    N = x.shape[0] # assumes shape[0] == shape[1]
    if (x.shape[0] != x.shape[1]) | (x.shape != y.shape):
        print('warning! non-square array! we assume square')

    narange = np.arange(N)
    mdiag = np.zeros((N,N),dtype='int8')
    xy4diag = np.zeros((N,2))
    matrot45clockwise = np.array([[1,-1],[1,1]])*0.5
    matrot45cclock = np.array([[1,1],[-1,1]])*0.5
    
    for i in np.arange(0,N):
        # first 45degrees up from x-axis
        xy4diag[:,0] = x[(narange+i)%N,narange]-(i/float(N)*boxlength)
        xy4diag[:,1] = y[(narange+i)%N,narange]
        xy4diag[N-i:,0] += boxlength
        xy4diag=np.dot(xy4diag,matrot45clockwise)
        mdiag[(narange+i)%N,narange] += foldfinder(xy4diag[:,0],boxlength)


        # now 45 degrees down from x-axis
        xy4diag[:,0] = x[(narange+i)%N,narange[::-1]]+(i/float(N)*boxlength)
        xy4diag[:,1] = y[(narange+i)%N,narange[::-1]]
        xy4diag[N-i:,0] += boxlength
        xy4diag=np.dot(xy4diag,matrot45cclock) # Could just compute based on x or y coordinate, but this may not be much slower

        mdiag[(narange+i)%N,narange[::-1]] += foldfinder(xy4diag[:,0],boxlength)
        
    return mdiag

def origami_order_tag(x,y,boxlength):
    maxpos = np.max((np.max(np.abs(x.flatten())),np.max(np.abs(y.flatten()))))
    if (maxpos > 2*boxlength):
        print('warning! the maximum position, ',maxpos,' > 2*boxlength; boxlength=',boxlength)
    if (maxpos < 0.5*boxlength):
        print('warning! the maximum position, ',maxpos,' < boxlength/2; boxlength=',boxlength)

    m_cart = m2D_cart(x,y, boxlength)
    m_45deg = m2D_45deg(x,y, boxlength)

    # take the maximum among Cartesian and 45-degree-diagonal directions
    return np.maximum(m_cart,m_45deg)
