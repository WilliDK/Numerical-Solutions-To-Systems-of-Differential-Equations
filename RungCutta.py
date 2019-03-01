from numpy import inf
from numpy import nan
import time

precision = 0.1
nConstants = 10

def getPrecision():
  return precision

def setPrecision(new):
  global precision
  precision = new

# This function can calculate Runge Cutta on one differential equation recursively (bottom up tail recursion)
def recursive(preCoord, n):
  # Base case: if n is 0
  if(n == 0):
    return (preCoord)
  # Alpha constants are calculated.
  a_1 = f(preCoord[0],preCoord[1])
  a_2 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_1*precision)
  a_3 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_2*precision)
  a_4 = f(preCoord[0] + precision, preCoord[1]+ a_3*precision)
  # The coordinate is calculated
  newCoord = [
    preCoord[0] + precision,
    preCoord[1] + 1/6 * (a_1 + 2*a_2 + 2*a_3 + a_4) * precision
  ]
  # tail recursive call
  return recursive(newCoord, n - 1)

# This function can calculate Runge Cutta on one differential equation iteratively.
def iterative_point(preCoord, n):
  # n times the following is done
  for i in range(n):
    # Alpha constants are calculated.
    a_1 = f(preCoord[0],preCoord[1])
    a_2 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_1*precision)
    a_3 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_2*precision)
    a_4 = f(preCoord[0] + precision, preCoord[1]+ a_3*precision)
    # The coordinate is calculated
    newCoord = [
      preCoord[0] + precision,
      preCoord[1] + 1/6 * (a_1 + 2*a_2 + 2*a_3 + a_4) * precision
    ]
    preCoord = newCoord
  # After n times the newcoord variable will contain the last calculated point (the n'th point).
  return newCoord

points = []
def iterative_points(preCoord, n):
  for i in range(n):
    a_1 = f(preCoord[0],preCoord[1])
    a_2 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_1*precision)
    a_3 = f(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_2*precision)
    a_4 = f(preCoord[0] + precision, preCoord[1]+ a_3*precision)
    newCoord = [
      preCoord[0] + precision,
      preCoord[1] + 1/6 * (a_1 + 2*a_2 + 2*a_3 + a_4) * precision
    ]
    points.append(newCoord)
    preCoord = newCoord

# This function calculates Runge Cutta on a system of differential equations. All points are stored in the global points variable.
def iterative_points_multiple(preCoords, n):
  global points
  points = []
  # preCoords are the first points known (ex: (x_0, y_0)
  points.append(preCoords)
  # This following is done n-1 one times (the initial preCoords are n = 1.
  for j in range(n-1):
    newCoords = []
    # For all the points in the list precoords.
    for i in range(len(preCoords)):
      preCoord = preCoords[i]
      # The constants, alpha_1 .. alpha_4 are calculated.
      a_1 = calcFs(preCoord[0],preCoord[1],functions[i])
      a_2 = calcFs(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_1*precision,functions[i])
      a_3 = calcFs(preCoord[0] + 0.5*precision, preCoord[1]+0.5*a_2*precision,functions[i])
      a_4 = calcFs(preCoord[0] + precision, preCoord[1]+ a_3*precision,functions[i])
      # The new point is calculated and appended to newCoords
      newCoords.append([
        preCoord[0] + precision,
        preCoord[1] + 1/6 * (a_1 + 2*a_2 + 2*a_3 + a_4) * precision
      ])
    # Now the new points calculates can be appended to the global points list.
    points.append(newCoords)
    # The previous coords are not the points we just calculated.
    preCoords = newCoords

F = "{1}*{0}**2 - {1}"

functions = []

# ability to set and format the function used in runge cutta.
def setF(s):
  i = 0
  while s.find("y", i) != -1:
    i = s.find("y")
    s = s[0:i] + "{1}" + s[i+1:len(s)]
  i = 0
  while s.find("x", i) != -1:
    i = s.index("x")
    s = s[0:i] + "{0}" + s[i+1:len(s)]
  F = s

# This function creates a list of properly formattet functions to use in Runge Cutta.
def setFs(l):
  for i in range(len(l)):
    s = l[i]
    i = 0
    while s.find("y", i) != -1:
      i = s.find("y")
      s = s[0:i] + "{1}" + s[i+1:len(s)]
    i = 0
    while s.find("x", i) != -1:
      i = s.index("x")
      s = s[0:i] + "{0}" + s[i+1:len(s)]
    functions.append(s)

# This function can evaluate the F String as a function.
def f(x,y):
  #return y*x**2 - y
  return eval(F.format(x,y))

# This function can evaluate a String as a function.
def calcFs(x,y,f):
  return eval(f.format(x,y))

'''
y(x) = 2*exp((1/3)*x*(x^2-3))
y(0) = 2
y(1) = 1.026834238
y(2) = 3.895468082 
y(3) = 806.8575870
'''

def generateConstants(n,p,preK):
  # the n'th constant is calculated
  k = f(p[0] + (1-(n/nConstants)) * precision, p[1] + (1-(n/nConstants)) * precision * preK)
  # base case. If n is 0 (all constants are calculated).
  if(n == 0):
    return (getWeight(n)*k,getWeight(n))
  generation = generateConstants(n-1,p,k)
  return (getWeight(n)*k + generation[0], getWeight(n)+generation[1])

# The get weight function returns the weight of the n'th constants.
def getWeight(n):
  return abs((1/2 - (1-(n/nConstants)))*2)*-1+2

# Works just like the iterative solution in the next bilag.
def experimental(preCoord, n):
  for i in range(n):
    # Besides for the fact that the alpha constants are calculated recursively with the generateConstants function
    generation = generateConstants(nConstants, preCoord, 1)
    newCoord = [
      preCoord[0] + precision,
      preCoord[1] + 1/generation[1] * generation[0] * precision
    ]
    preCoord = newCoord
  return newCoord

def printPoints():
  for i in range(len(points)):
    print(points[i][0], ",", points[i][1])

def printMultiplePoints():
  for i in range(len(points)):
    for l in points[i]:
      print(l[0], ",", l[1])

def printSection(n):
  for i in range(len(points)):
    print(points[i][n][0], ",", points[i][n][1])

def getPoints():
  return points

'''rec4((0, 2), 100)
printPoints()'''

'''
start = time.time()*1000
recursive((0,2),500)
end = time.time()*1000
print(end-start)
start = time.time()*1000
iterative_point((0,2),500)
end = time.time()*1000
print(end-start)
'''