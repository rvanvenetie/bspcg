import sys
import numpy as np
import matplotlib.pyplot as plt

def writemesh( fn, N, T, G):
  fp = open( fn, 'w')

  fp.write("%%sequential\n")
  fp.write("%d %d\n" % ( N.size/2, T.size/3))
  for i in range( N.size/2):
    fp.write( "%g %g %d\n" % (N[i][0], N[i][1], G[i]))
  for i in range( T.size/3):
    fp.write( "%d %d %d\n" % (T[i][0], T[i][1], T[i][2]))

def readmesh( fn):
    fp = open( fn)

    distributed = 0

    if fp.readline().find("distributed") > 0:
        distributed = 1

    if distributed:
        P = 0

    line = fp.readline().split(" ")
    nverts = int(line[0])
    ntris = int(line[1])
    if distributed:
        P = int(line[2])

    x = np.zeros(nverts)
    y = np.zeros(nverts)
    b = np.zeros(nverts)

    for i in range( nverts):
        [x[i], y[i], b[i]] = fp.readline().split(" ")

    t = np.zeros([ntris, 3])
    p = np.zeros(ntris)

    for i in range( ntris):
        if distributed:
            [t[i][0], t[i][1], t[i][2], p[i]] = fp.readline().split(" ")
        else:
            [t[i][0], t[i][1], t[i][2]] = fp.readline().split(" ")

    if distributed:
        return (1, nverts, ntris, x, y, b, t, P, p)
    else:
        return (0, nverts, ntris, x, y, b, t, 1, p)

if __name__ == "__main__":
    fn = sys.argv[1] if len( sys.argv) > 1 else "lshaped.m"

    mesh = readmesh( fn)
    print mesh
    distributed = mesh[0]


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect( 'equal')
    plt.tripcolor( mesh[3], mesh[4], mesh[6], 'go-', facecolors=mesh[8])
    plt.show()
