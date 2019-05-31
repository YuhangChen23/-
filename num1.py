# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:27:48 2019

@author: 威化饼干
"""

import numpy as np
import matplotlib.pyplot as plt

def fun(x):#原函数
    if x == 0:
        return 0
    a = np.sqrt(x)*np.log(x)
   # a = x**1.5
    return a


def trapz(n):
    sum1 = 0
    for i in range(n-1):
        sum1 = sum1+fun((i+1)/n)
    output = (1/n)*0.5*(0+fun(1)+2*sum1)
    return output


def fuhesimpson(n):
    sum1 = 0
    for i in range(n):
        sum1 = sum1 + fun((i/n)+0.5*1/n)
    sum2 = 0 
    for i in range(n-1):
        sum2 = sum2 + fun((i+1)/n)
    output = (1/6)*(1/n)*(0+fun(1)+4*sum1+2*sum2)
    return output   


def romberg(k):#n为迭代次数
    h =1
    T =[]#T表
    T.append([0])
    for i in range(k):
        sum1 = 0
        lie = []
        h = 1/(2**(i+1))      
        for j in range(2**i):#求中值点的和
            sum1 += fun((2*(j+1)-1)*h)
        first = 0.5*T[i][0]+h*sum1
        lie.append(first)
        for m in range(i+1):#对于每一行
            zhi = lie[m]+(lie[m]-T[i][m])/(4**(m+1)-1)
            lie.append(zhi)   
        T.append(lie)
    return T



'''
自适应辛普森公式
'''
def Simpson(left,right): 
    mid = (left+right)/2.
    return (1/6)*(fun(left)+fun(right)+4*fun(mid))*(right-left)
            

def adaptive_simpson(l,r,S):
    m = l+(r-l)/2.
    left_mianji = Simpson(l,m)
    right_mianji = Simpson(m,r)
    if abs(left_mianji+right_mianji-S) <= 10**(-4):
        return left_mianji+right_mianji
    else:
        return (adaptive_simpson(l,m,left_mianji)+adaptive_simpson(m,r,right_mianji))
    
    
          


a=adaptive_simpson(0,1,Simpson(0,1)) 
b = Simpson(0,1)   
print(a)
print(b)        
    



 

h=[]
for i in range(20):
    h.append(i+1)
    
h_test = []
for i in range(200):
    h_test.append(i+1)

value_trapz = []
#value_simp = []
#value_romberg = []
#value_correct=[]
#for i in h:
#    #value_trapz.append(round(trapz(i),5))
#    value_simp.append(fuhesimpson(i))
#    value_romberg.append(romberg(i)[i][i])
#    value_correct.append(-4/9)
    
for i in h_test:
    value_trapz.append(round(trapz(i),3))
    
for i in value_trapz:
    print(i)

plt.figure(figsize=(9,6),dpi=100)
plt.plot(h_test,value_trapz,color='green',label='Composite Trapezoidal')
#plt.plot(h,value_simp,color='blue',label="Composite Simpson")
#plt.plot(h,value_romberg,color='black',label='Romberg')
#plt.plot(h,value_correct,color='red',label='Correct Value')
plt.legend()
plt.xlabel('n:iterations')
plt.ylabel('integral value')
        

