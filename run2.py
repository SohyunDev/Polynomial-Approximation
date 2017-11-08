import numpy as np
import numpy.linalg as lin
import math
import matplotlib.pyplot as plt

def readData(dataFileNum):
    if(int(dataFileNum)==1):
        f=open("data1.txt","r")
    elif(int(dataFileNum)==2):
        f=open("data2.txt","r")
    elif(int(dataFileNum)==3):
        f=open("data3.txt","r")
    else:
        print("Error Input\n")
        return -1
    coordinates = []
    lines = f.readlines()
    for line in lines:
        line = line.replace("\n", "")
        if(int(dataFileNum)==1):
            line = line.split(" ")
        elif(int(dataFileNum)==2|int(dataFileNum)):
            line = line.split('\t')
        line = list(line)
        coordinate = (float(line[0]), float(line[1]))
        coordinates.append(coordinate)
    f.close()
    return coordinates

def makeMatrixA(dataCoordinates, degreeNum):
    matrix = [[]for _ in range(len(dataCoordinates))]
    for coordinate in dataCoordinates:
        for squared in range(0,degreeNum):
            matrix[dataCoordinates.index(coordinate)].append(pow(coordinate[0],squared))
    return matrix

def makeMatrixb(dataCoordinates, degreeNum):
    matrix = [[]for _ in range(len(dataCoordinates))]
    for coordinate in dataCoordinates:
        matrix[dataCoordinates.index(coordinate)].append(coordinate[1])
    return matrix

def calculateQiRin(Qi, Rin):
    result = [[]for _ in range (len(Qi))]
    for rowindex in range(0,len(Qi)):
        result[rowindex].append(Qi[rowindex][0]*Rin)
    return result

def subColumns(column1, column2):
    result = [[] for _ in range(len(column1))]
    if(len(column1)==len(column2)):
        for index in range(0,len(column1)):
            result[index].append(column1[index][0]-column2[index][0])
    return result

def getColumn(matrix, index):
    column = [[]for _ in range(len(matrix))]
    for rowindex in range(0,len(matrix)):
        column[rowindex].append(matrix[rowindex][index-1])
    return column

def getRnn(RnnParameter):
    result = 0
    for rowindex in range(0,len(RnnParameter)):
        result += pow(RnnParameter[rowindex][0],2)
    result = math.sqrt(result)
    return result

def getGetRnnParameter(matrixA, matrixQ, matrixR, n):
    result = getColumn(matrixA, n)
    if(n == 1):
        return result
    else:
        for index in range(1,n):
            Qi = getColumn(matrixQ, index)
            Rin = matrixR[index-1][n-1]
            result = subColumns(result, calculateQiRin(Qi, Rin))
    return result

def getQn(parameter, Rnn):
    qn = [[]for _ in range(len(parameter))]
    for rowindex in range(0,len(parameter)):
        qn[rowindex].append(parameter[rowindex][0]/Rnn)
    return qn

def getRij(qi, aj):
    result = 0
    if(len(qi)==len(aj)):
        for index in range(0,len(qi)):
            result += qi[index][0]*aj[index][0]
    else:
        return -1
    return result

def addColumnInMatrix(matrix,column):
    for rowindex in range(0,len(matrix)):
        matrix[rowindex].append(column[rowindex][0])

def RniQnCycle(matrixA, matrixQ, matrixR, n):
    columnNum = len(matrixA[0])
    parameter = getGetRnnParameter(matrixA, matrixQ, matrixR, n+1)
    # add Rnn
    Rnn = getRnn(parameter)
    matrixR[n].append(Rnn)
    # add qn
    qn = getQn(parameter, Rnn)
    addColumnInMatrix(matrixQ, qn)
    # add Rni
    if(n==columnNum-1):
        return
    else:
        for i in range(n+2,columnNum+1):
            ai = getColumn(matrixA, i)
            Rni = getRij(qn,ai)
            matrixR[n].append(Rni)
    return

def printMatrix(matrix):
    for row in matrix:
        print(row)

def getXSolution(matrixQ, matrixR, matrixb):
    matrixRInv = lin.inv(matrixR)
    matrixQTrans = np.transpose(matrixQ)
    result = np.dot(matrixRInv,matrixQTrans)
    result = np.dot(result,matrixb)
    return result

def getDot(polynomial, x):
    y = 0.0
    for index in range(0,len(polynomial)):
        y += polynomial[index][0]*pow(x,index)
    return (x, y)

def getDots(polynomial, maximumX,minimumX):
    dots = []
    for index in range(int(minimumX), int(maximumX)):
        for repeat in range(0,100):
            dot = getDot(polynomial, index+((1/100)*repeat))
            dots.append(dot)
    return dots

def printResult(dots, xPoints, yPoints):
    x = []
    y = []
    for index in range(0,len(dots)):
        x.append(dots[index][0])
        y.append(dots[index][1])
    plt.plot(x,y,'b')
    for index in range(0,len(xPoints)):
        plt.scatter(xPoints,yPoints,s=1,color='r')
    plt.show()

dataFileNum = input("Choose Datafile Number(1~ ) : ")
degreeNum = input("Number of degree : ")
dataCoordinates = readData(dataFileNum)

if dataFileNum != -1:
    # 행렬 초기화
    matrixA = makeMatrixA(dataCoordinates,int(degreeNum))
    matrixb = makeMatrixb(dataCoordinates,int(degreeNum))
    columnNum = len(matrixA[0])
    rowNum = len(matrixA)
    matrixQ = [[]for _ in range(rowNum)]
    matrixR = [[]for _ in range(columnNum)]

    for index in range(1,len(matrixA[0])):
        for repeat in range(0,index):
            matrixR[index].append(0)

    #1번~matrixA의 열갯수만큼 반복
    for index in range(0,columnNum):
        RniQnCycle(matrixA,matrixQ,matrixR,index)
    print("\n")

    print("MatrixA")
    printMatrix(matrixA)
    print("\n")

    print("Matrixb")
    printMatrix(matrixb)
    print("\n")

    print("MatrixQ")
    printMatrix(matrixQ)
    print("\n")

    print("MatrixR")
    printMatrix(matrixR)

    #R^-1Q^Tb
    result = getXSolution(matrixQ, matrixR, matrixb)
    print("\n")
    print("Result")
    printMatrix(result)

    xColumn = getColumn(matrixA,2)
    x =[]
    y=[]
    for index in range(0,len(xColumn)):
        x.append(int(xColumn[index][0]))
        y.append(int(matrixb[index][0]))
    polynomialDots = getDots(result,max(x),min(x))
    printResult(polynomialDots,x,y)

