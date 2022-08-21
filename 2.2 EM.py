import numpy as np

def EY(phi,L,h,i,Y,instruction):

    if instruction == "top":
        return -(phi[-1][i][int(Y/h) +1] - phi[-1][i][int(Y/h)])/h

    elif instruction == "bottom":
        return (phi[-1][i][int(Y/h) -1] - phi[-1][i][int(Y/h)])/h

def EX(L,h,DX,DY,w,e,X,Y):
    phi = potential(L,h,DX,DY,w,e)

    if X > L:
        return -(phi[-1][int(X/h)+1][int(Y/h)] - phi[-1][int(X/h)][int(Y/h)])/h

    elif X <= L and Y != 1:
        return -(phi[-1][int(X/h)+1][int(Y/h)] - phi[-1][int(X/h)][int(Y/h)])/h

    elif Y == 1 and X <= L:
        return 0

def E(X,Y):
    return (EX(X,Y),EY(X,Y))

def potential(L,h,DX,DY,w,e):

    M = int(DX/h)+1
    N = int(DY/h)+1

    "initial set-up"
    phi0 = [[0.1 for j in range(N)] for i in range(M)]

    "rules to set boundary condition"
    for i in range(M):
        if i == M-1:
            phi0[i]=[0 for j in range(N)]

    for j in range(N):
        if j == 0:
            for i in range(M):
                phi0[i][j] = 0

        if j == N-1:
            for i in range(M):
                phi0[i][j] = 0

        if j == int(1/h):
            for i in range(int(L/h)):
                phi0[i][j] = 1/2

    phi = [phi0]

    "first-step iterative averaging"

    phi2 = [[0.1 for j in range(N)] for i in range(M)]

    for i in range(0,M):
        if i == M-1:
            phi2[i]=[0 for j in range(N)]

        else:    
            for j in range(0,N):

                    if j == 0:
                        phi2[i][j] = 0

                    elif j == N-1:
                        phi2[i][j] = 0

                    elif j == int(1/h):
                        phi2[i][j] = 1/2

                    else:
                        if i == 0 and j != int(1/h):
                            phi2[0][j] = (phi[0][1][j]+phi[0][1][j]+phi[0][0][j-1]+phi[0][0][j+1])/4

                        if i > 0 and i <= int(L/h) and j != int(1/h):
                            phi2[i][j] = (phi[0][i-1][j]+phi[0][i+1][j]+phi[0][i][j-1]+phi[0][i][j+1])/4

                        if i > int(L/h):
                            phi2[i][j] = (phi[0][i-1][j]+phi[0][i+1][j]+phi[0][i][j-1]+phi[0][i][j+1])/4

    phi.append(phi2)



    "w: SOR parameter"

    "Number of iterations"
    K = 5000

    "e: Acceptable Error"

    for k in range(2,K):
        phi3 = [[0.1 for j in range(N)] for i in range(M)]

        for i in range(0,M):
            if i == M-1:
                phi3[i]=[0 for j in range(N)]

            else: 
                for j in range(0,N):

                    if j == 0:
                        phi3[i][j] = 0

                    elif j == N-1:
                        phi3[i][j] = 0

                    elif i <= int(L/h) and j == int(1/h):
                        phi3[i][j] = 1/2

                    else:
                        if i == 0 and j != int(1/h):
                            phi3[i][j] = (1-w) * phi[k-1][i][j] + (w/4)*(phi[k-1][1][j] + phi[k-1][1][j] + phi3[i][j-1] + phi[k-1][i][j+1])

                        if i > 0 and i <= int(L/h) and j != int(1/h):
                            phi3[i][j] = (1-w) * phi[k-1][i][j] + (w/4)*(phi3[i-1][j] + phi[k-1][i+1][j] + phi3[i][j-1] + phi[k-1][i][j+1])

                        if i > int(L/h):
                            phi3[i][j] = (1-w) * phi[k-1][i][j] + (w/4)*(phi3[i-1][j] + phi[k-1][i+1][j] + phi3[i][j-1] + phi[k-1][i][j+1])

        phi.append(phi3)

        s = 0
        for i in range(M):
            for j in range(N):
                s = s + np.abs(phi[k][i][j] - phi[k-1][i][j])

        r = s / (N*M)

        if r < e:
            print("stop at k=",k)
            break



    return [phi,k]


def charge(phi,L,h,DX,DY,w,e):
    MR = round((DX-1/4)/h)
    ML = - MR
    NB = round(1/(4*h))
    NT = round((DY-1/4)/h)
    Q = 0

    def pot(x,y):
        if x<0 and y <0:
            return phi[-x][-y]

        elif x<0 and y>=0:
            return phi[-x][y]

        elif x>=0 and y<0:
            return phi[x][-y] 

        elif x>=0 and y>=0:
            return phi[x][y] 
    

    for i in range(ML,MR+1):
        Q = Q + pot(i,NB) - pot(i,NB-1) + pot(i,NT) - pot(i,NT+1)

    for j in range(NB,NT+1):
        Q = Q + pot(ML,j) - pot(ML-1,j) + pot(MR,j) - pot(MR+1,j)
        
    return Q


"vertical electric field on a horizontal cut"

def hcutEyfield(phi,L,h,DX,DY,w,e,Y,instruction):
    
    M = int(DX/h)+1
    N = int(DY/h)+1
    eylist = []
    for i in range(0,M):
        eylist.append((i*h,EY(phi,L,h,i,Y,instruction)))

    return eylist


# Q2
"""
phi = potential(1,1/2,2,2,1,0.00001)[0][-1]

print(phi[1][1])
print(phi[1][3])

np.savetxt("2.csv", potential(1,1/2,2,2,1,0.000001)[0][-1], delimiter = ',')
"""

"Producing potential, horizontal cuts, and Summing flux for charge"
"""

phi_names = ["3.csv","4_8.csv","4_12.csv"]
hcut_names = ["hcutEfield.csv","hcutEfield_8.csv","hcutEfield_12.csv"]
top_names = ["plateTopEfield.csv","plateTopEfield_8.csv","plateTopEfield_12.csv"]
bottom_names = ["plateBottomEfield.csv","plateBottomEfield_8.csv","plateBottomEfield_12.csv"]
phi_storage = []


charges = []
for i in range(3):
    h = [1/4,1/8,1/12][i]
    phi = potential(2,h,4,4,1.8,0.000001)[0]


    phi_storage.append(phi[-1])
    np.savetxt(phi_names[i], phi[-1], delimiter = ',')
    np.savetxt(hcut_names[i], hcutEyfield(phi,2,h,4,4,1,0.000001,0,"top"), delimiter = ',')
    np.savetxt(top_names[i], hcutEyfield(phi,2,h,4,4,1,0.000001,1,"top"), delimiter = ',')
    np.savetxt(bottom_names[i], hcutEyfield(phi,2,h,4,4,1,0.000001,1,"bottom"), delimiter = ',')


    Q = charge(phi[-1],2,h,4,4,1.8,0.000001)
    print("when h = ",h," the charge is", Q)
    charges.append((h,Q))


for h in [1/16,1/24,1/32,1/40,1/48]:
    phi = potential(2,h,4,4,1.8,0.000001)[0]
    Q = charge(phi[-1],2,h,4,4,1.8,0.000001)
    print("when h = ",h," the charge is", Q)
    charges.append((h,Q))

#   print(charges)
np.savetxt("charges.csv", charges, delimiter = ',')
"""

"Investigation of how termination of k depends on w"
"""
print("Investigation of how termination of k depends on w")

term = []
for h in [1/4, 1/8, 1/12]:
    print("for h = ", h)
    for w in np.arange(1,2,0.1):
        print("for w = ", w)
        term.append([h,w,potential(2,h,4,4,w,0.000001)[1]])

np.savetxt("term.csv", term, delimiter = ',')

"""

"Investigation of increasing DX DY"

"""
np.savetxt("81.csv", potential(2,1/12,8,8,1.8,0.000001)[0][-1], delimiter = ',')
np.savetxt("82.csv", potential(2,1/12,8,12,1.8,0.000001)[0][-1], delimiter = ',')
np.savetxt("83small.csv", potential(2,1/12,2.5,1.5,1.8,0.000001)[0][-1], delimiter = ',')
"""


"Equipoential Lines"

#SECTION: COMPLEX NUMBER OPERATIONS [X,Y] = X + i Y

def multiply(c1,c2):
    x1 = c1[0]
    y1 = c1[1]
    x2 = c2[0]
    y2 = c2[1]

    return np.array((x1*x2 - y1*y2,x1*y2 + x2*y1))

def complex_exponential(c):
    x = c[0]
    y = c[1]

    return np.array((np.exp(x) * np.cos(y),np.exp(x) * np.sin(y)))

def divide(c1,c2):
    x1 = c1[0]
    y1 = c1[1]
    x2 = c2[0]
    y2 = c2[1]
    return multiply(c1,np.array((x2,-y2))) / (x2**2 + y2**2)



def position(p,s,L):
# This returns (X-L)+iY
    i = np.array((0,1))
    W = np.array((-p,s))
    
    return (np.array((1,0)) + complex_exponential(-2*np.pi*multiply(i,W)))/np.pi - 2* multiply(i,W) + np.array((L,0))

def equipotential_line(p,L):
    return [position(p,i,L) for i in np.arange(-10,10,0.01)]

def field_line(s,L):
    return [position(i,s,L) for i in np.arange(-0.5,0.51,0.01)]

def Efield(p,s,L):
# This returns Ex - i Ey
    i = np.array((0,1))
    W = np.array((-p,s))
    return divide(np.array((0,1)), 2 * (complex_exponential(-2*np.pi*multiply(i,W)) + np.array((1,0))))

AnalyticalTop = []
for s in np.arange(0,5,0.01):
    AnalyticalTop.append((-position(1/2,s,0)[0], -Efield(1/2,s,0)[1]))

AnalyticalBottom = []
for s in np.arange(-50,0,0.01):
    AnalyticalBottom.append((-position(1/2,s,0)[0], -Efield(1/2,s,0)[1]))

np.savetxt("AnalyticalTop.csv",AnalyticalTop, delimiter = ',')
np.savetxt("AnalyticalBottom.csv",AnalyticalBottom, delimiter = ',')

# np.savetxt("equipotential48.csv", equipotential_line(0.48,0), delimiter = ',')

"""
for k in range(5):
    p = [-0.42,-0.44,-0.46,-0.48,-0.49][k]
    s = [-4,-2,-3,0.2,0.4][k]
    equipotential_name = ["equipotential42.csv","equipotential44.csv","equipotential46.csv","equipotential48.csv","equipotential49.csv"]
    field_name = ["field-40.csv","field-20.csv","field-30.csv","field02.csv","field04.csv"]
    np.savetxt(equipotential_name[k], equipotential_line(p,0), delimiter = ',')
    np.savetxt(field_name[k], field_line(s,0), delimiter = ',')

"""

"""
phi = potential(10,1/4,25,25,1.7,0.000001)[0]
hcutEyfield(phi,10,1/4,25,25,1.7,0.000001,1,"top")
"""

"Investigation of large values of L, DX, DY"

phi_names = ["11-05.csv","11-10.csv","11-20.csv","11-50.csv"]
top_names = ["Q11plateTopEfield_05.csv","Q11plateTopEfield_10.csv","Q11plateTopEfield_20.csv","Q11plateTopEfield_50.csv"]
bottom_names = ["Q11plateBottomEfield_05.csv","Q11plateBottomEfield_10.csv","Q11plateBottomEfield_20.csv","Q11plateBottomEfield_50.csv"]
phi_storage = []


for i in range(4):
    h = 1/2
    L = [5,10,20,50][i]
    DX = 100
    DY = 100
    phi = potential(L,h,DX,DY,1.8,0.0001)[0]


    phi_storage.append(phi[-1])
    np.savetxt(phi_names[i], phi[-1], delimiter = ',')
    np.savetxt(top_names[i], hcutEyfield(phi,L,h,DX,DY,1.8,0.0001,1,"top"), delimiter = ',')
    np.savetxt(bottom_names[i], hcutEyfield(phi,L,h,DX,DY,1.8,0.0001,1,"bottom"), delimiter = ',')





