# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


#implementation of BAM
import pprint
import numpy
import numpy as np


def init_BAM(data):
	ab = []
	for ele in data:
		ab.append([generate_bipolar_form(ele[0]), generate_bipolar_form(ele[1])])

	x_length = len(ab[0][1])
	y_length = len(ab[0][0])

	_bam = [] #initialise empty bam array
	temp = []
	for ele in range(y_length):
		temp = [0] * x_length
		_bam.append(temp)

	return ab, x_length, y_length, _bam

def create_BAM(ab, _bam):
	for ele in ab:
		X = ele[0]
		Y = ele[1]
		for ix, xi in enumerate(X):
			for iy, yi in enumerate(Y):
				_bam[ix][iy] += xi * yi
	return _bam

def get_associations(A, _bam, x_length, y_length):
	A = multiply_vec(A, _bam, x_length, y_length)
	return threshold(A)

def multiply_vec(vec, _bam, x_length, y_length):
	result = [0] * x_length
	for x in range(x_length):
		for y in range(y_length):
			result[x] += vec[y] * _bam[y][x]
	return result
	

def testTargets(y, weight):
	# Multiply the target pattern with the weight matrix
    # (weight X y)
    x = np.dot(weight, y)
    x[x <= 0] = -1
    x[x > 0] = 1
    return np.array(x)

def transform_to_bipolar_vector(A):
    X = A.flatten()
    for i in range(X.size):
        if X.item(i) == 0:
            X.itemset(i, -1)
    X = X[:,np.newaxis]
    return X
def generate_bipolar_form(vec):

	result = []
	for ele in vec:
		if ele == 0:
			result.append(-1)
		else:
			result.append(1)
	return result

def threshold(vec):
	result = []
	for ele in vec:
		if ele < 0:
			result.append(0)
		else:
			result.append(1)
	return result


def read_matrix(file, n, m):
	matrix = np.zeros((n, m))

	for i in range(matrix.shape[0]):
		matrix[i] = [int(x) for x in file.readline().split()]

	return matrix


def write_matrix(matrix, file):
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			file.write(str(matrix.item((i, j))) + ' ')
		file.write('\n')
	return


print("\nTesting for input patterns: Set A")


def testInputs(x, weight):
	# Multiply the input pattern with the weight matrix
	# (weight.T X x)
	y = np.dot(weight.T, x)
	y[y < 0] = -1
	y[y >= 0] = 1
	return np.array(y)

def TestPattern():

	weight = weightMatrixCalculate()
	x1, x2, x3, x4 = xValue()
	print("\nOutput of input pattern 1")
	print(testInputs(x1, weight))
	print("\nOutput of input pattern 2")
	print(testInputs(x2, weight))
	print("\nOutput of input pattern 3")
	print(testInputs(x3, weight))
	print("\nOutput of input pattern 4")
	print(testInputs(x4, weight))



# Test for Target Patterns: Set B
print("\nTesting for target patterns: Set B")


def testTargets(y, weight):
	# Multiply the target pattern with the weight matrix
	# (weight X y)
	x = np.dot(weight, y)
	x[x <= 0] = -1
	x[x > 0] = 1
	return np.array(x)



def TargetPattern():

	weight=weightMatrixCalculate()
	y1,y2,y3,y4=yValue()
	print("\nOutput of target pattern 1")
	print(testTargets(y1, weight))
	print("\nOutput of target pattern 2")
	print(testTargets(y2, weight))
	print("\nOutput of target pattern 3")
	print(testTargets(y3, weight))
	print("\nOutput of target pattern 4")
	print(testTargets(y4, weight))

def xValue():

	x1 = np.array([1, 1, 1, 1, 1, 1]).reshape(6, 1)
	x2 = np.array([-1, -1, -1, -1, -1, -1]).reshape(6, 1)
	x3 = np.array([1, 1, -1, -1, 1, 1]).reshape(6, 1)
	x4 = np.array([-1, -1, 1, 1, -1, -1]).reshape(6, 1)

	return x1, x2, x3, x4

def yValue():

	y1 = np.array([1, 1, 1]).reshape(3, 1)
	y2 = np.array([-1, -1, -1]).reshape(3, 1)
	y3 = np.array([1, -1, 1]).reshape(3, 1)
	y4 = np.array([-1, 1, -1]).reshape(3, 1)

	return y1, y2, y3, y4


def weightMatrixCalculate():

	x1, x2, x3, x4=xValue()
	y1, y2, y3, y4 = yValue()

	inputSet = np.concatenate((x1, x2, x3, x4), axis=1)
	targetSet = np.concatenate((y1.T, y2.T, y3.T, y4.T), axis=0)
	print("\nWeight matrix:")
	weight = np.dot(inputSet, targetSet)
	print(weight)

	print("\n------------------------------")

	return weight

def main():


	print("sonrası")

	fileA = open("assets/5_1r", 'r')
	listOfLines = fileA.readlines()
	nA = [[int(num) for num in line.split()] for line in listOfLines]
	fileA.close()
	print("A",nA)

	
	fileB = open("assets/5_1", 'r')
	listOfLinesB = fileB.readlines()
	nB = [[int(num) for num in line.split()] for line in listOfLinesB]
	fileB.close()
	print("B",nB)
	
	data_associationsFile = [nA, nB]
	abFile, x_lengthFile, y_lengthFile, bamFile = init_BAM(data_associationsFile)
	bam_matrixFile = create_BAM(abFile, bamFile)
	ppFile = pprint.PrettyPrinter(indent=4)
	print("Bam Matrix: ")
	ppFile.pprint(bam_matrixFile)


	print("[1, 0, 1, 0, 1, 0] :- ", get_associations(nA[0], bam_matrixFile, x_lengthFile, y_lengthFile))
	print("[1, 1, 1, 0, 0, 0] :- ", get_associations(nB[0], bam_matrixFile, x_lengthFile, y_lengthFile))




if __name__ == '__main__':
	main()
	TestPattern()
	TargetPattern()
