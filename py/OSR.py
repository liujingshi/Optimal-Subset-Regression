
import numpy as np
from math import log

class OSR:

    def average(self, x, l, i):
        sumi = 0
        num = int(len(x) / l)
        for j in range(num):
            sumi += x[i - 1 + (j * l)]
        return sumi / num
    
    def MGF(self, x, n):
        F = []
        m = int(n / 3)
        for l in range(1, m + 1):
            H = []
            while len(H) < n:
                i = 1
                while i <= l and len(H) < n:
                    H.append(self.average(x, l, i))
                    i += 1
            F.append(H)
        return F

    def T(self, arr):
        return np.array(arr).T.tolist()

    def diff(self, x):
        return np.diff(np.array(x)).tolist()

    def sumAddSon(self, x, l, t):
        sumi = 0
        for i in range(t):
            sumi += x[l][i + 1]
        return sumi

    def sumAdd(self, x, x1):
        xa = x
        for i in range(len(xa)):
            xa[i][0] = x1
        for l in range(len(xa)):
            for i in range(1, len(xa[l])):
                xa[l][i] = x1 + self.sumAddSon(x, l, i)
        return xa

    def calcQx(self, x):
        xb = self.average(x, 1, 1)
        print(x[0])
        n = len(x)
        sumi = 0
        for t in range(n):
            sumi += (x[t] - xb)**2
        return sumi / n

    def calcQ(self, x, xt):
        n = len(xt)
        sumi = 0
        for t in range(n):
            sumi += (x[t] - xt[t])**2
        return sumi / n

    def calcS1(self, x, xk):
        Qx = self.calcQx(x)
        Qk = self.calcQ(x, xk)
        n = len(x)
        return n * (1 - (Qk / Qx))

    def calcUV(self, x):
        sumi = 0
        for i in x:
            sumi += abs(i)
        return sumi / len(x)

    def calcT(self, x, u, t):
        return 0 if x[t] > u else (1 if abs(x[t]) <= u else (2 if x[t] < -u else 0))

    def calcNij(self, x, f):
        x1 = self.diff(x)
        f1 = self.diff(f)
        U = self.calcUV(x1)
        V = self.calcUV(f1)
        Nij = [[0 for i in range(3)] for j in range(3)]
        for t in range(len(f1)):
            i = self.calcT(x1, U, t)
            j = self.calcT(f1, V, t)
            Nij[i][j] += 1
        return Nij

    def calcS2(self, x, f):
        Nij = self.calcNij(x, f)
        R1 = R2 = R3 = Ni = Nj = tem = 0
        n = len(x)
        for i in Nij:
            tem = 0
            for j in i:
                if j > 0:
                    tem += j * log(j)
            R1 += tem
        for i in Nij:
            Ni = 0
            for j in i:
                Ni += j
            if Ni > 0:
                R2 = Ni * log(Ni)
        for j in range(len(Nij)):
            Nj = 0
            for i in range(len(Nij)):
                Nj += Nij[i][j]
            if Nj > 0:
                R3 += Nj * log(Nj)
        return 2 * (R1 + (n - 1) * log(n - 1) - R2 - R3)

    def calcCSCson(self, x, f):
        S1 = self.calcS1(x, f)
        S2 = self.calcS2(x, f)
        return S1 + S2

    def calcCSC(self, x, F):
        f = self.comeBack(x, F)
        return self.calcCSCson(x, f)

    def calcCSCsOne(self, x, F):
        CSC = []
        for i in F:
            CSC.append(self.calcCSC(x, [i]))
        return CSC

    def calcCSCs(self, x, F):
        CSC = []
        for i in F:
            CSC.append(self.calcCSC(x, i))
        return CSC

    def matmul(self, arr1, arr2):
        return np.matmul(np.array(arr1), np.array(arr2)).tolist()

    def inv(self, arr):
        return np.linalg.inv(np.array(arr)).tolist()

    def Xb(self, F):
        result = []
        temp = [1 for i in range(len(F[0]))]
        result.append(temp)
        for i in F:
            result.append(i)
        return self.T(result)

    def comeBackP(self, y, Xb):
        p = [y]
        p = self.T(p)
        XbT = self.T(Xb)
        XbTMXb = self.matmul(XbT, Xb)
        XbTMXbI = self.inv(XbTMXb)
        XbTMXbIMXbT = self.matmul(XbTMXbI, XbT)
        return self.matmul(XbTMXbIMXbT, p)

    def comeBack(self, x, F):
        Xb = self.Xb(F)
        a = self.comeBackP(x, Xb)
        result = []
        Frow = len(F)
        Fcol = len(F[0])
        for i in range(Fcol):
            res = a[0][0]
            for j in range(Frow):
                res += a[j + 1][0] * F[j][i]
            result.append(res)
        return result

    def rouSelectSon(self, c, x):
        div = []
        for i in range(len(c)):
            if c[i] > x:
                div.append(i)
        return div

    def rouSelect(self, F, CSC, x):
        P = []
        for i in range(len(CSC)):
            temCSC = self.rouSelectSon(CSC[i], x)
            for j in range(len(temCSC)):
                P.append(F[i][temCSC[j]])
        return P

    def twoAddOne(self, two):
        n = len(two)
        t = 0
        two[n - 1] += 1
        for i in range(n)[::-1]:
            two[i] += t
            if two[i] > 1:
                t = 1
                two[i] = 0
            else:
                t = 0
                break
        return two

    def twoIsFull(self, two):
        for i in two:
            if i == 0:
                return False
        return True

    def group(self, P):
        son = []
        two = [0 for i in range(len(P))]
        while not self.twoIsFull(two):
            tem = []
            two = self.twoAddOne(two)
            for i in range(len(P)):
                if two[i] == 0:
                    continue
                tem.append(P[i])
            son.append(tem)
        return son

    def findMaxCSC(self, CSC):
        res = 0
        maxCSC = CSC[0]
        for i in range(1, len(CSC)):
            if CSC[i] > maxCSC:
                maxCSC = CSC[i]
                res = i
        return res

    def predictSon(self, A, son):
        result = []
        Frow = len(son)
        Fcol = len(son[0])
        for i in range(len(Fcol)):
            res = A[0][0]
            for j in range(len(Frow)):
                res += A[j + 1][0] * son[j][i]
            result.append(res)
        return result

    def predictEasy(self, x, q):
        n = len(x)
        F0 = self.MGF(x, n)
        x1 = self.diff(x)
        F1 = self.MGF(x1, n)
        x2 = self.diff(x1)
        F2 = self.MGF(x2, n)
        F3 = self.sumAdd(F1, x[0])
        F0p = self.MGF(x, n + q)
        F1p = self.MGF(x1, n + q)
        F2p = self.MGF(x2, n + q)
        F3p = self.sumAdd(F1p, x[0])
        CSC0 = self.calcCSCsOne(x, F0)
        CSC1 = self.calcCSCsOne(x, F1)
        CSC2 = self.calcCSCsOne(x, F2)
        CSC3 = self.calcCSCsOne(x, F3)
        F = [F0, F1, F2, F3]
        CSC = [CSC0, CSC1, CSC2, CSC3]
        print(CSC0[0])
        xx = 11.07
        xxs = [40.113, 38.885, 37.652, 36.415, 35.172, 33.924, 32.761, 31.41, 30.144, 28.869, 27.587, 26.296, 24.996, 23.685, 22.262, 21.026, 19.675, 18.307, 16.919, 15.507, 14.067, 12.592, 11.07, 9.488, 7.815, 5.911, 3.841]
        P = self.rouSelect(F, CSC, xx)
        for i in xxs:
            P = self.rouSelect(F, CSC, i)
            if len(P) > 6:
                break
        print(len(P))
        son = self.group(P)
        lastCSClist = self.calcCSCs(x, son)
        lastCSCi = self.findMaxCSC(lastCSClist)
        lastSon = lastCSClist[lastCSCi]
        lastXb = self.Xb(lastSon)
        lastA = self.comeBackP(x, lastXb)
        newSon = []
        for i in range(len(lastSon)):
            isFind = False
            for j in range(len(F0)):
                hasDif = False
                for k in range(len(F0[j])):
                    if abs(lastSon[i][k] - F0[j][k]) < 0.0001:
                        continue
                    hasDif = True
                if hasDif:
                    continue
                else:
                    isFind = True
                    newSon.append(F0p[j])
                    break
            if isFind:
                continue
            for j in range(len(F1)):
                hasDif = False
                for k in range(len(F1[j])):
                    if abs(lastSon[i][k] - F1[j][k]) < 0.0001:
                        continue
                    hasDif = True
                if hasDif:
                    continue
                else:
                    isFind = True
                    newSon.append(F1p[j])
                    break
            if isFind:
                continue
            for j in range(len(F2)):
                hasDif = False
                for k in range(len(F2[j])):
                    if abs(lastSon[i][k] - F2[j][k]) < 0.0001:
                        continue
                    hasDif = True
                if hasDif:
                    continue
                else:
                    isFind = True
                    newSon.append(F2p[j])
                    break
            if isFind:
                continue
            for j in range(len(F3)):
                hasDif = False
                for k in range(len(F3[j])):
                    if abs(lastSon[i][k] - F3[j][k]) < 0.0001:
                        continue
                    hasDif = True
                if hasDif:
                    continue
                else:
                    isFind = True
                    newSon.append(F3p[j])
                    break
            if isFind:
                continue
        return self.predictSon(lastA, newSon)

    def predict(self, x, q):
        nhz = self.predictEasy(x, q)
        e = []
        for i in range(len(x)):
            e.append(x[i] - nhz[i])
        ePre = self.predictEasy(e, q)
        result = []
        for i in range(len(x), len(nhz)):
            result.append(nhz[i] + ePre[i])
        return result
