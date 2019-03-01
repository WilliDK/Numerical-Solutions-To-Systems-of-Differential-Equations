from NumericalSolutionsToDifferentialEquations import RungCutta
from NumericalSolutionsToDifferentialEquations import Visualize
from NumericalSolutionsToDifferentialEquations import Euler
from NumericalSolutionsToDifferentialEquations import CustomRegression
from NumericalSolutionsToDifferentialEquations import DataHandler
import numpy as np

#inputs
'''
a = 0.001
N = 1000
I_0 = 1
R_0 = 0
y = 0.33
''''''
N = 3550
I_0 = (0,50)
R_0 = (0,0)
y = 0.55
a = 1.065 / N
S_0 = (0,N - I_0[1] - R_0[1])
deathRate = 0.0155'''

N = 14865
I_0 = (0,1)
R_0 = (0,1)
y = 3.2041
a = 0.000227
S_0 = (0, N - I_0[1] - R_0[1])
dfr = 0.0000637145587767
dfr = 0.0155
d = 0.304454
'''
S_0 = (0,481000)
I_0 = (0,100)
R_0 = (0,18)
y = 0.2
N = S_0[1] + I_0[1] + R_0[1]
a = (R_0[1]*y)/N
deathRate = 0.0155'''


def runcuttaSIR(n):
    RungCutta.setFs([
        str(-1*a) + "*y*points[len(points)-1][1][1]",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y),
        #str(y)+"*points[len(points)-1][1][1]"
    ])
    RungCutta.setPrecision(0.1)
    RungCutta.iterative_points_multiple((
        S_0,
        I_0,
        #R_0
    ), n)
    Rs = [1] + [N - a[0][1] - a[1][1] for a in RungCutta.getPoints()[1::]]
    print(Rs)
    DataHandler.writeToFile(RungCutta.getPoints(), RungCutta.getPrecision())
    #RungCutta.printMultiplePoints()
    Visualize.plot([Visualize.convert(
        RungCutta.getPoints())[1]],
        ("dage", "antal personer"),
        ["smittede"],
        "SIR model for Kolera Epidemien på Vor Frelser Sogn, 1853, smittede"
    )
    Visualize.plot(Visualize.convert(
            RungCutta.getPoints()),
           ("dage", "antal personer"),
           ("smitbare", "smittede", "imune/døde"),
            "SIR model for Kolera Epidemien på Vor Frelser Sogn, 1853"
    )

def runcuttaSIRregression(n):
    RungCutta.setFs([
        str(-1*a) + "*y*points[len(points)-1][1][1]",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y),
        str(y)+"*points[len(points)-1][1][1]"
    ])
    RungCutta.setPrecision(0.1)
    RungCutta.iterative_points_multiple((
        S_0,
        I_0,
        R_0
    ), n)
    # S(t) regression
    functionString = "S(t) = c-(a/(1+c*e^(-b*(t-d)))"
    constantsString = ("a", "b", "c", "d")
    CustomRegression.setFunction(lambda x, a, b, c, d : (c-(a/ (1 + c*np.exp(-b*(x-d))))))
    Sconstants = CustomRegression.calcConstants(CustomRegression.convert(RungCutta.getPoints(), 0), [-3360, -0.001, 250, 90])
    linearSpaceS = [CustomRegression.calcLinearSpace(Sconstants, S_0[0], n * RungCutta.getPrecision())]
    print(CustomRegression.formatFormula(constantsString, Sconstants, functionString))
    xdata = np.array(CustomRegression.convert(RungCutta.getPoints(), 0)[0])
    ydata = CustomRegression.convert(RungCutta.getPoints(), 0)[1]
    residuals = np.array(ydata) - np.array(CustomRegression.function(xdata, *Sconstants))
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    print("R-squared: " + str(1 - (ss_res / ss_tot)))

    # R(t) regression
    functionString = "R(t) = a/(1+c*e^(-b*t))"
    constantsString = ("a", "b", "c")
    CustomRegression.setFunction(lambda x, a, b, c: (a / (1 + c * np.exp(-b * x))))
    Rconstants = CustomRegression.calcConstants(CustomRegression.convert(RungCutta.getPoints(), 2),[10, -0.001, 250])
    linearSpaceR = [CustomRegression.calcLinearSpace(Rconstants, S_0[0], n * RungCutta.getPrecision())]
    print(CustomRegression.formatFormula(constantsString, Rconstants, functionString))
    xdata = np.array(CustomRegression.convert(RungCutta.getPoints(), 2)[0])
    ydata = CustomRegression.convert(RungCutta.getPoints(), 2)[1]
    residuals = np.array(ydata) - np.array(CustomRegression.function(xdata, *Rconstants))
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    print("R-squared: " + str(1 - (ss_res / ss_tot)))

    # I(t)
    functionString = "I(t) = " + str(N) + " - a/(1+c*e^(-b*t)) - (g-(d/(1+g*e^(-f*(t-h))))"
    constantsString = ("a", "b", "c", "d", "f", "g", "h")
    CustomRegression.setFunction(lambda x, a, b, c, d, f, g, h: N - ((a) / (1 + c * np.exp(-b * x))) - ((g-(d/ (1 + g*np.exp(-f*(x-h)))))))
    linearSpaceI = [CustomRegression.calcLinearSpace(Rconstants + Sconstants, S_0[0], n * RungCutta.getPrecision())]
    print(CustomRegression.formatFormula(constantsString, Rconstants + Sconstants, functionString))
    xdata = np.array(CustomRegression.convert(RungCutta.getPoints(), 1)[0])
    ydata = CustomRegression.convert(RungCutta.getPoints(), 1)[1]
    residuals = np.array(ydata) - np.array(CustomRegression.function(xdata, *[*Rconstants,*Sconstants]))
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    print("R-squared: " + str(1 - (ss_res / ss_tot)))
    Visualize.plot(linearSpaceI,
           ("time", "people"),
           ("smitbare", "smittede", "imune/døde", "S(t) regression", "I(t) regression", "R(t) regression"),
                   "SIR og regression for Kolera Epidemien på Vor Frelser Sogn, 1853"
   )

def eulerSIR(n):
    Euler.setFs([
        str(-1*a) + "*y*points[len(points)-1][1][1]",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y),
        str(y)+"*points[len(points)-1][1][1]"
    ])
    Euler.iterative_points_multiple((
        S_0,
        I_0,
        R_0
    ), n)
    Euler.printSection(0)
    Euler.printSection(1)
    Euler.printSection(2)
    Visualize.plot(Visualize.convert(Euler.getPoints()), ("time", "people"), ["smitbare", "smittede", "imune/døde"])

def runcuttaSIR_Reallocation(n):
    RungCutta.setFs([
        str(dfr*N)+     "-"+str(a)+"*y*points[len(points)-1][1][1]"+"-"+str(dfr)+"*y",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y)+"-y*"+str(dfr),
        str(y)+"*points[len(points)-1][1][1]"+"-y*"+str(dfr)
    ])
    RungCutta.iterative_points_multiple((
        S_0,
        I_0,
        R_0
    ), n)
    #RungCutta.printSection(0)
    #RungCutta.printSection(1)
    #RungCutta.printSection(2)
    Visualize.plot([Visualize.convert(RungCutta.getPoints())[1]], ("dage", "antal personer"), ["smittede"],"med omfordeling")

def runcuttaSIR_Reallocation_Deaths(n):
    RungCutta.setFs([
        str(dfr)+"*("+str(N)+"-points[len(points)-1][3][1])"+     "-"+str(a)+"*y*points[len(points)-1][1][1]"+        "-"+str(dfr)+"*y",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y)+         "-y*"+str(dfr),
        str(y)+"*points[len(points)-1][1][1]"+          "-y*"+str(dfr),
        "("+str(y)+"*points[len(points)-1][1][1]"+          "-y*"+str(dfr)+")*"+str(d)
    ])
    RungCutta.iterative_points_multiple((
        S_0,
        I_0,
        R_0,
        (0,1)
    ), n)
    print(RungCutta.getPoints())
    #RungCutta.printSection(0)
    #RungCutta.printSection(1)
    #RungCutta.printSection(2)
    Visualize.plot([Visualize.convert(RungCutta.getPoints())[1]], ("dage", "antal personer"), ["smittede", "smittede", "imune/døde", "døde"],"med døde")

def runcuttaSIR_Deaths(n):
    RungCutta.setFs([
        "-"+str(a)+"*y*points[len(points)-1][1][1]",
        str(a)+"*points[len(points)-1][0][1]*y-y*"+str(y),
        str(y)+"*points[len(points)-1][1][1]",
        "("+str(y)+"*points[len(points)-1][1][1]"+")*"+str(d)
    ])
    RungCutta.iterative_points_multiple((
        S_0,
        I_0,
        R_0,
        (0,1)
    ), n)
    print(RungCutta.getPoints())
    #RungCutta.printSection(0)
    #RungCutta.printSection(1)
    #RungCutta.printSection(2)
    Visualize.plot(Visualize.convert(RungCutta.getPoints()), ("dage", "antal personer"), ["smitbare", "smittede", "imune/døde", "døde"],"med døde")

# kald en hvilken som helst konfiguration (funktion) her.
runcuttaSIR_Reallocation(10000)