import sys
import fileinput
import os
import matplotlib.pyplot as plt
def main():
	x = int(sys.argv[1])#Initial KPOINTS value
	y = int(sys.argv[2])#Final KPOINTS value
	z = int(sys.argv[3])#Stepping value
	kPoints = []
	for i in range(x, y + 1, z):
		kPoints.append(i)
	eCutOff = [] 
	forces = [] 
	for i in range(0, len(kPoints)):
		dir = "kMesh" + str(kPoints[i]) +"/out-jwl"
		with open(dir, 'r') as file:
			filedata = file.readlines()
		temp = filedata[len(filedata) - 2]
#		print(temp)
		temp2 = temp.split()
#		print(temp2)
		forces.append(float(temp2[2]))
	print(forces)
	makePlot(kPoints, forces)

def makePlot(kPoints, forces):
    plt.plot(kPoints, forces)
    plt.xlabel('KPOINTS')
    plt.ylabel('Energy')
    plt.show()

main()
