import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt

# def create_L(X):
#     midL = np.zeros((2, 2))
#     radL = np.zeros((2, 2))
#     radL[0][0] = (0.8 * X.b[0] * 2 - 0.8 * X.a[0]*2)/2
#     midL[0][0] = 0.8 * X.a[0]*2 + radL[0][0]

#     radL[0][1] = (1.5 * X.b[1]*2 - 1.5 * X.a[1]*2)/2
#     midL[0][1] = 1.5 *X.a[1]*2 + radL[0][1]
    
#     radL[1][0] = (-1/X.b[0] + 1/X.a[0])/2
#     midL[1][0] = 1/X.b[0] + radL[1][0]

#     midL[1][1] = -1 
#     radL[1][1] = 0
#     return ip.Interval(midL, radL, midRadQ=True)

def create_L(X):
    print('=========')
    print(X)
    print(X.a)
    print(X.b)
    print('=========')
    midL = np.zeros((2, 2))
    radL = np.zeros((2, 2))
    radL[0][0] = (4 * X.b[0] - 4 * X.a[0])/2
    midL[0][0] = 4 * X.a[0] + radL[0][0]

    radL[0][1] = (2 * X.b[1] - 2 * X.a[1])/2
    midL[0][1] = 2 * X.a[1] + radL[0][1]
    
    radL[1][0] = (-1/(X.b[0]+2)**2 + 1/(X.a[0]+2)**2)/2
    midL[1][0] = -1/(X.a[0]+2)**2 + radL[1][0]

    midL[1][1] = -1 
    radL[1][1] = 0
    L = ip.Interval(midL, radL, midRadQ=True)
    print(L)
    return L

# def F_x(x):
#     F = np.zeros(2)
#     F[0] = 0.8 * x[0]**2 + 1.5 * x[1]**2 - 2
#     F[1] = np.log(x[0]) - x[1] + 1
#     return F

def F_x(x):
    F = np.zeros(2)
    F[0] = 2 * x[0]**2 + x[1]**2 - 2
    F[1] = 1/(x[0]+2) - x[1]
    return F

def inv_midA(A):
    midA = np.zeros((2, 2))
    for i in range (0,2):
        for j in range (0,2):
            midA[i][j] = A[i][j].mid
    return np.linalg.inv(midA)

def mid_X(X):
    midX = np.zeros(2)
    midX[0] = X[0].mid
    midX[1] = X[1].mid
    return midX



def matrix_norm(A):
    norm = 0
    sum = 0
    for i in range(0,2):
        sum = A[i][0].b + A[i][1].b
        if(sum > norm):
            norm = sum
    return norm

def b_norm(A):
    return A[1]

def dist(X, Y):
    d = np.zeros(2)
    for i in range(0,2):
        if(max(abs(X[i].b - Y[i].b), abs(X[i].a - Y[i].a)) > 0.0000001):
            return False
    return True

def Dist(X, Y):
    d = np.zeros(2)
    for i in range(0,2):
        d[i] = max(abs(X[i].b - Y[i].b), abs(X[i].a - Y[i].a))
    return d

def intersection(X, kr):
    midL = np.zeros(2)
    radL = np.zeros(2)
    a1 = kr.a[0] if (kr.a[0] >= X.a[0]) else X.a[0]
    b1 = kr.b[0] if (kr.b[0] <= X.b[0]) else X.b[0]
    a2 = kr.a[1] if (kr.a[1] >= X.a[1]) else X.a[1]
    b2 = kr.b[1] if (kr.b[1] <= X.b[1]) else X.b[1]

    radL[0] = (b1 - a1)/2
    midL[0] = a1 + radL[0]
    radL[1] = (b2 - a2)/2
    midL[1] = a2 + radL[1]
    print('L = ', radL[0], midL[0], radL[1], midL[1])
    return ip.Interval(midL, radL, midRadQ=True)


if __name__ == "__main__":
    midX = [1, 0]
    radX = [0.2, 0.2]
    I = [[1,0], [0,1]]
    X = ip.Interval(midX, radX, midRadQ=True)
    X = ip.Interval([[0.5, 1], [0, 1]])
    print(X)
    X_k = []
    while(True):
        X_k.append(X)
        A = create_L(X)
        x_av = mid_X(X)
        b = F_x(x_av)
        lambd = inv_midA(A)
        eta = matrix_norm(I - lambd @ A)
        teta = max(abs(lambd @ b))/(1 - eta)
        mid_kr = [0, 0]
        rad_kr = [teta, teta]
        print('rad_kr = ', rad_kr)
        X_kr = ip.Interval(mid_kr, rad_kr, midRadQ=True)
        kr = (lambd @ b + (I - lambd @ A) @ X_kr)
        new_X_kr = intersection(kr, X_kr)
        while(not dist(X_kr, new_X_kr)):
            X_kr = new_X_kr
            kr = (lambd @ b + (I - lambd @ A) @ X_kr)
            print('ints: ', (kr, X_kr))
            new_X_kr = intersection(kr, X_kr)
        X_kr = new_X_kr 
        N = x_av - X_kr
        print(N)
        newX = intersection(X, N)
        if(dist(X, newX)):
            print(Dist(X, newX))
            X = newX
            break
        X = newX
    print(newX)

    # kr_k.append(kravchik(x_av, F, Lambda, I, X, L))

    fig, ax = plt.subplots(figsize=(4, 4))


    for i in range(len(X_k)):
        # kr_one = abs(kr_k[i][0].b - kr_k[i][0].a) 
        # kr_two = abs(kr_k[i][1].b - kr_k[i][1].a) 
        one = abs(X_k[i][0].b - X_k[i][0].a)
        two = abs(X_k[i][1].b - X_k[i][1].a)
        # Rect = plt.Rectangle((kr_k[i][0].a, kr_k[i][1].a), kr_one , kr_two, edgecolor='red', facecolor='none', linewidth=0.7)
        
        iveRect = plt.Rectangle((X_k[i][0].a, X_k[i][1].a), one , two, edgecolor='black', facecolor='none', linewidth=1.5)
      
        plt.gca().add_patch(iveRect)
        # plt.gca().add_patch(Rect)

    x = np.arange(-1, 1, 0.0001)
    yp = np.sqrt((2 - 2*pow(x, 2)))
    yn = - np.sqrt((2 - 2*pow(x, 2)))
    # yp2 = np.sqrt((3 - 0.8*pow(x, 2))/1.5)
    # yn2 = - np.sqrt((3 - 0.8*pow(x, 2))/1.5)
    # x = np.arange(-10, 15)

    # y = 1 - x
    # plt.plot(x, y, color='g', linewidth=0.7, label = 'x + y = 1')
    plt.plot(x, yp, color='g', linewidth=1, label = r'$2x^2+y^2-2=0$')
    plt.plot(x, yn, color='g', linewidth=1)
    # plt.plot(x, yp2, color='g', linewidth=1)
    # plt.plot(x, yn2, color='g', linewidth=1)

    x = np.arange(0.001, 2, 0.001)
    y = 1/(x+2)
    # y2 = np.log(x)

    plt.plot(x, y, color='red', linewidth=1, label = r'$\frac{1}{x+2}-y=0$')
    # plt.plot(x, y2, color='red', linewidth=1)
    # y = 3 - x
    # plt.plot(x, y, color='olive', linewidth=0.7, label = 'x + y = 3')

    # y = x
    # plt.plot(x, y, color='blue', linewidth=0.7, label = 'x/y = 1')
    
    # y = x/6
    # plt.plot(x, y, color='darkblue', linewidth=0.7, label = 'x/y = 6')


    plt.grid()
    # plt.xlim(-5, 3)
    # plt.ylim(-2, 5.5)
    plt.xlim(-3, 3)
    plt.ylim(-2, 2)
    plt.legend()
    plt.show()