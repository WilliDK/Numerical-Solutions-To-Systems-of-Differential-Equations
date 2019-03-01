from NumericalSolutionsToDifferentialEquations import RungCutta, Euler
from NumericalSolutionsToDifferentialEquations import CustomRegression
import numpy as np
import math

# Reads data from file and converts to pairs of integers
def readData(filename):
    data = []
    with open(filename, mode='r') as file:
        data = list(map(lambda x:[int(x.split(",")[0]), int(x.split(",")[1])], file.readlines()))
    return data

# Essentially just adds up all y values sequentially.
def accumulateData(data):
    sum = 0
    for l in data:
        sum += l[1]
        l[1] = sum
    return data

# Data is read and accumulated
data = accumulateData(readData("data.csv"))
print(data)

# Function used for regression
def functionRegression(x, a, gamma, d):
    # First S and I is calculated with Runge-cutta
    RungCutta.setFs([
        str(-1 * a) + "*y*points[len(points)-1][1][1]",
        str(a) + "*points[len(points)-1][0][1]*y-y*" + str(gamma)
    ])
    x = list(x)
    RungCutta.iterative_points_multiple((
        (0,14865-2),
        (0,1)
    ), len(x))
    result = RungCutta.getPoints()
    # R is calculated since: R = N - S - I
    Rs = [1]+[14865 - a[0][1] - a[1][1] for a in result[1::]]
    # The new y values are calculcated (y is D, since we need to compare the result of this function with the dataset).
    y = [1] + [r*d for r in Rs[1::]]
    #print(a, gamma, d)
    return y

# Converts 2d arrays of pairs to 2d array like this: [xs,ys]
def customConvert(data):
    xs = []
    ys = []
    for point in data:
        xs.append(point[0])
        ys.append(point[1])
    return [xs, ys]

# This function controls the code.
def Regression(data):
    # The function of the customregression script is set to functionRegression
    CustomRegression.setFunction(functionRegression)
    # Precision must be the same as delta x in the dataset
    RungCutta.setPrecision(1)
    # Deviation before regression
    residuals = np.array(customConvert(data)[1]) - np.array(functionRegression(customConvert(data)[0], 0.000227, 3.2, 0.3))
    ss_res = np.sum(residuals ** 2)
    print("SUM " + str(ss_res))
    #Parameters are calculated
    constants = CustomRegression.calcConstants(customConvert(data), [0.000227, 3.2, 0.3])
    print("Found")
    # Deviation after regression
    residuals = np.array(customConvert(data)[1]) - np.array(functionRegression(customConvert(data)[0], constants[0], constants[1], constants[2]))
    ss_res = np.sum(residuals ** 2)
    print("SUM " + str(ss_res))
    return constants

print(Regression(accumulateData(readData('data.csv'))))

'''
#print(functionRegression([i for i in range(74)], 3.4239/14865, 3.2478, 0.3231))

def gradient_descent(data, parameters, learningRate):
    # I loop over every x and y value in the dataset
    for x,y in data:
        # A cost is calculated based on the current parameters
        c = cost(parameters, data, x)
        #All parameteres are updated based on this cost
        for i in range(len(parameters)):
            parameters[i] -= learningRate/(len(data)) * c
    return parameters

def cost(parameters, data, x):
    xs = [x for x, y in data]
    # New y values are calculated
    calculatedData = functionRegression(xs, parameters[0], parameters[1], parameters[2])
    sum = 0
    # The deviation is summed
    for i in range(len(data)):
        sum += ((calculatedData[data[i][0]] - data[i][1])*x)
    return sum / (2*len(data))

# This function is called.
def findParameters(intervals, precision):
    # Based on an initial interval a set of 10 values are chosen for each parameter.
    aVals = getValues(intervals[0])
    bVals = getValues(intervals[1])
    cVals = getValues(intervals[2])
    values = []
    # The higher the precision, the more precise the resulting paremeters are going to be.
    for i in range(precision):
        # Based on all combinations of parameter values, the combination with the least squared deviation is chosen recursively.
        optimalCombo = completeSearch([aVals, bVals, cVals], 0, [])
        values = [
            aVals[optimalCombo[0]],
            bVals[optimalCombo[1]],
            cVals[optimalCombo[2]]
        ]   
        # Based on these new intervals returned by completeSearch 10 new values are chosen within these intervals.
        aVals = getValues([aVals[optimalCombo[0]-1], aVals[optimalCombo[0]+1]])
        bVals = getValues([bVals[optimalCombo[0]-1], bVals[optimalCombo[0]+1]])
        cVals = getValues([cVals[optimalCombo[0]-1], cVals[optimalCombo[0]+1]])
    return values

# This function takes an interval and returns nVals evenly distributed values within this interval (in fact it is 12 but the 2 on the edges of the interval are just boundaries).
def getValues(interval):
    nVals = 10
    valueSpan = (interval[1]-interval[0])/nVals
    values = []
    for i in range(nVals+1):
        values.append(interval[0] + valueSpan * i)
    return values

# This functions works through all parameter values. Input sizes are fixed and therefore there will never be more than 10^3 combinations.
def completeSearch(values, index, combination):
    # Base case. All parameters are included in this combination.
    if(index == len(values)):
        return combination
    minDeviation = float('inf')
    minCombination = []
    # For all parameter values of a specific parameter the following is computed.
    for i in range(1,len(values[index])-2):
        # The combination based on the current value of the parameter is calculated.
        tempCombination = completeSearch(values, index + 1, combination + [i])
        if len(tempCombination) == 0: continue
        # Deviation of this set of parameters is calculated.
        tempDeviation = squaredDeviation(tempCombination)
        # If the deviation is less than the current minimum deviation the minDeviation is updated and the minCombination is set to be the current Combination.
        if tempDeviation < minDeviation: #it could be equal and then another branch should start
            minDeviation = tempDeviation
            minCombination = tempCombination
    return minCombination

# This function returns a summed squared deviation based on a set of parameters. The function used is functionRegression.
def squaredDeviation(calcParameters):
    global data
    reelData = [y for x,y in data]
    calculatedData = functionRegression([x for x,y in data], calcParameters[0], calcParameters[1], calcParameters[2])
    sum = 0.0
    for i in range(len(reelData)-1):
        sum += ((reelData[i]-calculatedData[i])*(reelData[i]-calculatedData[i]))
    if(math.isnan(sum)):
        return float('inf')
    return sum

print(gradient_descent(accumulateData(readData('data.csv')), [0.0001, 3, 0.5], 1))
#print(Regression(accumulateData(readData('data.csv'))))
'''
'''print(findParameters([
    (1,3),
    (1,3),
    (0,1)
], 10))'''