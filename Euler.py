precision = 0.1

def setPrecision(num):
  global precision
  precision = num

def iterative_point(preCoord, n):
  for i in range(n):
    newCoord = [
      preCoord[0] + precision,
      preCoord[1] + f(preCoord[0], preCoord[1])*precision
    ]
    preCoord = newCoord
  return newCoord

points = []

def iterative_points_multiple(preCoords, n):
  global points
  points = []
  points.append(preCoords)
  for i in range(n):
    newCoords = []
    for j in range(len(preCoords)):
      preCoord = preCoords[j]
      newCoords.append([
        preCoord[0] + precision,
        preCoord[1] + calcFs(preCoord[0], preCoord[1], functions[j]) * precision
      ])
    points.append(newCoords)
    preCoords = newCoords

functions = []

def f(x,y):
  return y*x**2 - y

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

def calcFs(x,y,f):
  return eval(f.format(x,y))

def printMultiplePoints():
  for i in range(len(points)):
    for l in points[i]:
      print(l[0], ",", l[1])

def getPoints():
  return points

def printSection(n):
  for i in range(len(points)):
    print(points[i][n][0], ",", points[i][n][1])

#print(iterative_point((0,2),1))