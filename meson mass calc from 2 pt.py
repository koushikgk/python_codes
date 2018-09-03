import numpy as np
import pandas as pd
import os
import fnmatch
import re
import math
import matplotlib.pyplot as plt
import decimal

if __name__ == "__main__":
    regexp = re.compile(r'meson([-.0-9-e]+)')
    regexp1 = re.compile(r'real=\[ *([-.0-9-e]+)')
    regexp2 = re.compile(r'([-.0-9-e]+),')
    regexp3 = re.compile(r'([-.0-9-e]+)]')
    corr = np.array([[]])
    gamma = [0]
    count = 0
    l = len(gamma)
    for filename in os.listdir('.'):
        if fnmatch.fnmatch(filename, 'pt_m0p1*'):
            with open(filename, 'r') as file:
                temp = np.array([decimal.Decimal(0) for i in range(32)])
                for line in file:
                    match = regexp.search(line)
                    if match:
                        for i in gamma:
                            if float(match.group(1))==i:
                                match1 = regexp1.search(line)
                                if match1:
                                    c = 0
                                    c = np.array(re.findall(regexp2, line))
                                    c = np.append(c, np.array(re.findall(regexp3, line)))
                                    c = np.array([decimal.Decimal(i) for i in c])
                                    temp = temp + c
                                    count = count+1
                                    print(line, c, count)
                                    if c.any()==float(0):
                                        print(c, file)

                if temp[1] == decimal.Decimal(0):
                    print(file)
                temp = temp/l
                print(temp)
                n_t = len(temp)
                corr = np.append(corr, temp)

    n_rows = int(len(corr)/n_t)                    
    corr = np.reshape(corr, (n_rows, n_t))

    corr1 = np.array([decimal.Decimal(0) for i in range(n_rows*int(n_t/2+1))])
    corr1 = np.reshape(corr1, (n_rows, int(n_t/2+1)))

    for i in range(n_rows):
        for j in range(int(n_t/2+1)):
            if j==0:
                corr1[i][j] = corr[i][j]
            else :
                corr1[i][j] = decimal.Decimal(0.5)*(corr[i][j] + corr[i][n_t-j])

    print("\nfolded correlator avg and std error\n")
    corr2 = np.array([decimal.Decimal(0) for i in range(int(n_t/2+1))])
    for i in range(int(n_t/2+1)):
        for j in range(n_rows):
            corr2[i] = corr2[i]+corr1[j][i]/n_rows

    r = np.array([])
    for j in range(int(n_t/2-1)):
        if(j==0):
            r = np.append(r, decimal.Decimal(0))
        else:
            r = np.append(r, decimal.Decimal(math.acosh(decimal.Decimal(0.5)*(corr2[j]+corr2[j+2])/corr2[j+1])))

    print("\njackknife estimate and std er\n")

    c = np.array([decimal.Decimal(0) for i in range(n_rows*int(n_t/2+1))])
    c = np.reshape(c, (n_rows, int(n_t/2+1)))

    for i in range(n_rows):
        for k in range(int(n_t/2+1)):
            temp = decimal.Decimal(0)
            for j in range(n_rows):
                if(j!=i):
                    temp = temp + corr1[j][k]
            temp = temp/(n_rows-1)
            c[i][k] = temp

    r1 = np.array([decimal.Decimal(0) for i in range(n_rows*int(n_t/2-1))])
    r1 = np.reshape(r1, (n_rows, int(n_t/2-1)))

    for i in range(n_rows):
        for j in range(int(n_t/2)-1):
            if j==0:
                r1[i][j] = decimal.Decimal(0)
            else:
                r1[i][j] = decimal.Decimal(math.acosh(decimal.Decimal(0.5)*(c[i][j]+c[i][j+2])/c[i][j+1]))

    r_er1 = np.array([])
    r_m = np.array([])
    m = np.array([])

    for i in range(int(n_t/2)-1):
        temp1 = decimal.Decimal(0)
        temp = decimal.Decimal(0)
        for j in range(n_rows):
            temp = temp + r1[j][i]
            temp1 = temp1 + (r1[j][i]-r[i])*(r1[j][i]-r[i])
        r_er1 = np.append(r_er1, decimal.Decimal(math.sqrt(temp1/n_rows)))
        r_m = np.append(r_m, decimal.Decimal(temp/n_rows))

    t = np.array([float(i) for i in range(int(n_t/2-1))])
    for i in range(int(n_t/2-1)):
        print(i, r[i], r_er1[i]*decimal.Decimal(math.sqrt(n_rows-1)))

    plt.xlabel('at')
    plt.ylabel('am')
    plt.ylim(0.7, 1.3)
    plt.xlim(0,15)
    plt.title("am_av_ll")
    plt.grid()
    plt.scatter(t, r)
    plt.errorbar(t,r, linestyle='None', yerr=r_er1*decimal.Decimal(math.sqrt(n_rows-1)))
    plt.show()
    plt.savefig('1.eps', format='eps', dpi=1000)
    plt.close() 

