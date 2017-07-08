# -*- coding: utf-8 -*-
'''
Created on 2016-12-4
author: shaodonghang
'''
#WW Model snow albedo model
import  math
import  numpy            as     np
from    numpy            import array
from    random           import random
from    math             import sin, sqrt
import  time        
import  datetime         as     datetime
from    datetime         import timedelta
import  netCDF4          as     nc
from    netCDF4          import Dataset
from    decimal          import *
import  matplotlib.pylab as     plt
import  sympy
import  scipy.integrate  as integrate
import  scipy.special    as special
#利用MIE散射理论模拟单个粒子的光学特性
#根据经验公式可近似模拟单次散射反照率w,不对称参数g,消光截面Q_ext
w     = 1-exp(a0+a1*r**0.5+a2*r)
g     = b0+b1*r**0.5+b2*r
Q_ext = 2+exp(c0+c1*r**0.5+c2*r)
t1    = 3*SWE*Q_ext/(4*r*P_ice)
#SWE为雪水当量，P_ice为纯冰密度=917kg.m-3
#&-Eddington近似
g = g/(1+g)
w = (1-g**2)*w1/(1-w1*g**2)
t = (1-w1*g**2)*t1
#当入射光源为直射入射时
a = 1-w*g
b = g/a
s = (3*a*(1-w))**0.5
P = 2s/(3*a)
y = (1-A)/(1+A)
Q1= (y+P)*exp(s*t)
Q2= (y-P)*exp(-s*t)
Q = (1+P)*Q1-(1-P)*Q2
as_μ0   = (w*(1-b*s*μ0))/((1+P)*(1+s*μ0))
Q_as_μ0 = 2*(P*(1-y+w*b)+w*(1+b)*(y*s*μ0-P)/(1-s**2*μ0**2))*exp(-t/μ0)-w*b*(Q1-Q2)+w*(1+b)*(Q1/(1+s*μ0)-Q2/(1-s*μ0))

#当入射光源全为漫射入射时，以Q量表示的雪面反照率ad及t→∞时，ad_∞为
Q_ad = 2*P*((1-y+w*b)*(1-t)-y*w*(1+b)/(1-w))*exp(-t)-2*P*(w*(1+b)*(2/s**2+y*t/(1-w)+(1-y+w*b)*t**2))*Ei*(-t)+\
    2*w*(1+b)/s**2*(Q1*(Ei*(-(1+s)*t)+s-np.log(1+s))-Q2*(Ei*(-(1-s)*t)-s-np.log(1-s)))-w*b*(Q1-Q2)
#Wiscombe and Warren
#黑天空反照率R0
R0 = (2*y2*s*exp(-t/μ0)*(rs-w*R/(1-S**2*μ0**2))+w*(Q1*P1-Q2*P2))/(Q1*(y1+s)-Q2*(y1-s))
#t为比例光学厚度，w为单次散射反照率，g为不对称参数
t=t1*(1-w1*g1**2)
w=w1*(1-g1**2)/(1-w1*g1**2)
g=g1/(1+g1**2)
#&-Eddington近似参数
y1=(7-w*(4+3*g))/4
y2=(7-w*(4-3*g))/4
y3=(2-3*g*μ0)/4
y4=1-y3
#其他参数设置,rs是下层雪的朗伯反射率
Q1=exp(s*t)*(y2-rs*(y1+s))
Q2=exp(-s*t)*(y2-rs*(y1-s))
P1=(a2+s*y3)/(1+s*μ0)
P2=(a2-s*y3)/(1-s*μ0)
R =μ0*(rs*a1-a2)+rs*y4+y3
a1=y1*y4+y2*y3
a2=y2*y4+y2*y3
s =(y1**2+y2**2)**0.5
#m为散射辐射的比例，Rw为白天空反照率
integrate_Rw = 2*R0*np.cos(solar_zenith_angle)*np.sin(solar_zenith_angle)
Rw=integrate.quad(integrate_Rw, 0, pi)
#总反照率是直射辐射和散射辐射的加权
R_all=m*Rw+R0*(1-m)
#Choudhury和Chang创建了考虑粗糙度的修订模型
