import numpy as np
from scipy.optimize import curve_fit

#Converts data to proper format ([xs,ys])
def convert(points, index):
    new = []
    for i in range(len(points[0])):
        new.append([])
        for j in range(len(points[0][0])):
            new[i].append([])
    for pset in points:
        for i, point in enumerate(pset):
            for j, value in enumerate(point):
                new[i][j].append(point[j])
    return new[index]

#Sets the function to be performed regression on
def setFunction(f):
     global function
     function = f

def function1(x, *para):
    x = x / 1000
    y = function(x, *para) * 1000
    return y

def function2(x, *para):
    x = x * 1000
    y = function(x, *para) / 1000
    return y

#SAMPLES
#normaliser
function = lambda x, a, b, c, d : (c-(a/ (1 + c*np.exp(-b*(x-d)))))
#function = lambda x, a, b, c : ((a) / (1 + c*np.exp(-b*x)))
'''def function(x, a, b, c, d):
    #print("PRINTING:", a,b,c)
    return (c-(a/ (1 + c*np.exp(-b*(x-d)/100))))*100

def function1(x, a, b, c, d):
    #print("PRINTING:", a,b,c)
    return (c-(a/ (1 + c*np.exp(-b*(x-d)*100))))/100'''
#obs: removed sublist 1,len


#Custom regression. Der bruges scipy.
def calcConstants(data, initialGuess):
    xs = np.array([el for el in data[0]][::]).astype(dtype=np.float64)
    ys = np.array([el for el in data[1]][::]).astype(dtype=np.float64)
    return list(curve_fit(function, xs, ys, initialGuess, method='trf')[0])

def calcLinearSpace(constants, min, max):
    frekvens = 0.1
    x = [i/(1/frekvens) for i in range(int(min),int((int(max))*1/frekvens))]
    y = [function(i, *constants) for i in x]
    return [x,y]

def formatConstants(names, constants):
    result = ""
    for i in range(len(names)):
        result += names[i] + ": " + str(constants[i]) + "   "
    return result

def formatFormula(names, constants, f):
    for i in range(len(names)):
        f = f.replace(names[i],str(constants[i]))
    return f
