import matplotlib.pyplot as plt

def convert(points):
    new = []
    for i in range(len(points[0])):
        new.append([])
        for j in range(len(points[0][0])):
            new[i].append([])
    for pset in points:
        for i, point in enumerate(pset):
            for j, value in enumerate(point):
                new[i][j].append(point[j])
    return new

def plot(points, axis_labels, function_labels, title):
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.title(title)
    for i, pointlist in enumerate(points):
        plt.scatter(pointlist[0], pointlist[1], s=10)
    plt.legend(function_labels)
    plt.show()
