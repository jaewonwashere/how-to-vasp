import os
import sys
import fileinput
import subprocess
import time


def runTasks():#Main Call
	#command line arguments
	x = int(sys.argv[1])#Initial KPOINTS value
	y = int(sys.argv[2])#Final KPOINTS value
	z = int(sys.argv[3])#Stepping value

	dirX = os.getcwd()

	kPoints = []
	for i in range(x, y + 1, z):
		kPoints.append(i)

	createKPOINTS()#Creates the KPOINTS AutoMesh file

	time.sleep(1)	

	makeDir(kPoints)#Creates Dir to run the jobs in

	time.sleep(1)	#For some reason this is necessary

	qsubJobs(kPoints)#Runs Jobs

	rmKPOINTS(dirX)#Removes file "createKPOINTS" made	


def makeDir(kPoints):  #Copies all the necessary job files into their new Dir
	dirpath = os.getcwd()
	for i in range(0, len(kPoints)):
		tempI = "kMesh" + str(kPoints[i])
		os.makedirs(tempI)#Makes new dir
		subprocess.Popen(['cp','INCAR','job','KPOINTS','POSCAR','POTCAR', dirpath + "/" +tempI])


def qsubJobs(kPoints): #Changes the KPOINTS values and qsubs the jobs in the dir
	dirpath = os.getcwd()
	for i in range(0, len(kPoints)):
		print("Job " + str(i + 1))
		os.chdir(dirpath + "/kMesh" + str(kPoints[i]))
		with open('KPOINTS', 'r') as file:
			filedata = file.read()
		filedata = filedata.replace('1',str(kPoints[i]))
		with open('KPOINTS', 'w') as file:
			file.write(filedata)
		
		time.sleep(1)#Probably not necessary	
	
		process = subprocess.Popen(['qsub','job'], stdout=subprocess.PIPE)
		process.wait()
		print("Job "+str(i + 1) + " Running")


def createKPOINTS():
	filedata = "Automatic mesh\n0              ! number of k-points = 0 ->automatic generation scheme\nAuto           ! fully automatic\n  1           ! length (l)"
	with open('KPOINTS', 'w') as file:
		file.write(filedata)

def rmKPOINTS(x):
	os.chdir(x)
	subprocess.Popen(['rm','KPOINTS'])
	

runTasks()

