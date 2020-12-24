import heapq
import random
import matplotlib.pyplot as plt
import numpy as py
from scipy import stats
constant_multiplyer = 65539
alpha = 0.1

def Random_Number_Generation(Zi_1,n):
    randomNums = []
    randomNums.append(Zi_1 / (2 ** 31))

    for i in range(n-1):
        Zi = (constant_multiplyer * Zi_1) % (2 ** 31)
        Ui = Zi / (2 ** 31)
        randomNums.append(Ui)
        Zi_1 = Zi

    Zi = (65539 * Zi_1) % (2 ** 31)
    return randomNums, Zi

def criticalpoint_chiSquare(k,alpha):
    critical_chiSquare = stats.chi2.ppf(q=1-alpha, df=k-1)
    return critical_chiSquare

def Uniformity_Test(randomNums,n,k):
    subInterval = []
    for i in range(k+1):
        subInterval.append(i/k)

    fi = []
    for i in range(k+1):
        fi.append(0)

    i = 0
    while(i < len(randomNums)):
        for j in range(1,k+1):
            if (randomNums[i] >= subInterval[j-1]) and (randomNums[i] < subInterval[j]):
                fi[j] += 1
                break

        i += 1

    sum = 0.0
    for i in range(1,k+1):
        sum += (fi[i] - (n/k)) ** 2

    chiSquare = sum*k/n
    critical_chiSquare = criticalpoint_chiSquare(k,alpha)

    if chiSquare > critical_chiSquare:
        print(chiSquare,"/Rejected")
    else:
        print(chiSquare,"/Not Rejected")

def fj_generate(fj1_d__str,s,k,D):
    tempstr = s
    if len(s)<D:
        for i in range(1,k+1):
            tempstr += str(i)
            fj_generate(fj1_d__str,tempstr,k,D)
            tempstr = s
    else:
        fj1_d__str.append(s)


def Serial_Test(randomNums,N,K,D):
    Ul = []
    count = 0
    for i in range(int(N/D)):
        inner = []
        for j in range(D):
            inner.append(randomNums[count])
            count += 1
        Ul.append(inner)
    #print(randomNums)
    #print(Ul)

    subInterval = []
    for i in range(K + 1):
        subInterval.append(i / K)

    #print(subInterval)

    s = ""
    fj1_d__str = []
    fj_generate(fj1_d__str,s, K, D)
    #print(fj1_d__str)

    fj1_d = []
    for i in range(K**D):
        fj1_d.append(0)

    for l in range(len(Ul)):
        string = ""
        for i in range(len(Ul[l])):
            for j in range(1, K + 1):
                if (Ul[l][i] >= subInterval[j - 1]) and (Ul[l][i] < subInterval[j]):
                    string += str(j)
                    break
        ind = fj1_d__str.index(string)
        fj1_d[ind] += 1

    sum = 0.0
    for i in range(K**D):
        sum += (fj1_d[i] - (N / (K**D))) ** 2

    chiSquare = sum * (K**D) / N
    critical_chiSquare = criticalpoint_chiSquare((K**D), alpha)

    if chiSquare > critical_chiSquare:
        print(chiSquare,"/Rejected")
    else:
        print(chiSquare,"/Not Rejected")

def Runs_Test(randomNums,N):
    a = [[4529.4, 9044.9, 13568, 18091, 22615, 27892],
         [9044.9, 18097, 27139, 36187, 45234, 55789],
         [13568, 27139, 40721, 54281, 67852, 83685],
         [18091, 36187, 54281, 72414, 90470, 111580],
         [22615, 45234, 67852, 90470, 113262, 139476],
         [27892, 55789, 83685, 111580, 139476, 172860]]
    b = [1 / 6, 5 / 24, 11 / 120, 19 / 720, 29 / 5040, 1 / 840]

    r = [0,0,0,0,0,0]

    length = 1
    for i in range(1, len(randomNums)):
        if randomNums[i] >= randomNums[i - 1]:
            length += 1
        else:
            if(length > 6):
                r[5] += 1
            else:
                r[length-1] += 1
            length = 1

        if(i == len(randomNums)-1):
            if (length > 6):
                r[5] += 1
            else:
                r[length-1] += 1

    sum = 0.0
    for i in range(6):
        for j in range(6):
            sum += (a[i][j] * (r[i] - N*b[i]) * (r[j] - N*b[j]))

    Run = sum / N
    critical_chiSquare = criticalpoint_chiSquare(7, alpha)

    if Run >= critical_chiSquare:
        print(Run,"/Rejected")
    else:
        print(Run,"/Not Rejected")

def Correlation_Test(randomNums,N,J):
    h = int((N-1)/J - 1)
    sum = 0.0
    for k in range(h):
        sum += (randomNums[k*J] * randomNums[(k+1)*J])

    rho = (12*sum/(h+1)) - 3
    Var_rho = (13*h + 7) / ((h+1)**2)
    A = rho / (Var_rho ** 0.5)
    Z = stats.norm.ppf(q= 1-(alpha/2))

    if abs(A) > Z:
        print(abs(A),"/Rejected")
    else:
        print(abs(A),"/Not Rejected")



if __name__ == "__main__":
    seed = 1505084
    print("Uniformity Test")
    n = [20, 500, 4000, 10000]
    k = [10, 20]
    for i in range(len(n)):
        for j in range(len(k)):
            randomNums, seed = Random_Number_Generation(seed, n[i])
            Uniformity_Test(randomNums, n[i], k[j])

    print("\n")

    print("Serial Test")
    n = [20, 500, 4000, 10000]
    k = [4, 8]
    d = [2, 3]
    for i in range(len(n)):
        for l in range(len(d)):
            for j in range(len(k)):
                randomNums, seed = Random_Number_Generation(seed, n[i])
                Serial_Test(randomNums, n[i], k[j], d[l])

    print("\n")

    print("Runs Test")
    n = [20, 500, 4000, 10000]
    for i in range(len(n)):
        randomNums, seed = Random_Number_Generation(seed, n[i])
        Runs_Test(randomNums, n[i])

    print("\n")

    print("Correlation Test")
    n = [20, 500, 4000, 10000]
    J = [1, 3, 5]
    for i in range(len(n)):
        for j in range(len(J)):
            randomNums, seed = Random_Number_Generation(seed, n[i])
            Correlation_Test(randomNums, n[i], J[j])
