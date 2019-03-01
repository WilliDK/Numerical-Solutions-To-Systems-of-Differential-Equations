def writeToFile(data, precision):
  with open("result.csv", "w") as file:
    for i in range(len(data)):
      file.write(str(precision * i) + ",")
      l = data[i]
      for j in range(len(l)):
        file.write(str(l[j][1])+("," if j != len(l)-1 else ""))
      file.write("\n")
