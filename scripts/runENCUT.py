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

	ENCUTs = []
	for i in range(x, y + 1, z):
		ENCUTs.append(i)

#	createINCAR()#Creates the KPOINTS AutoMesh file

	time.sleep(1)	

	makeDir(ENCUTs)#Creates Dir to run the jobs in

	time.sleep(1)	#For some reason this is necessary

	qsubJobs(ENCUTs)#Runs Jobs

#	rmINCAR(dirX)#Removes file "createKPOINTS" made	


def makeDir(kPoints):  #Copies all the necessary job files into their new Dir
	dirpath = os.getcwd()
	for i in range(0, len(kPoints)):
		tempI = "ENCUT" + str(kPoints[i])
		os.makedirs(tempI)#Makes new dir
		subprocess.Popen(['cp','INCAR','job','KPOINTS','POSCAR','POTCAR', dirpath + "/" +tempI])


def qsubJobs(kPoints): #Changes the KPOINTS values and qsubs the jobs in the dir
	dirpath = os.getcwd()
	for i in range(0, len(kPoints)):
		print("Job " + str(i + 1))
		os.chdir(dirpath + "/ENCUT" + str(kPoints[i]))
		with open('INCAR', 'r') as file:
			filedata = file.read()
		filedata = filedata.replace('SENTINEL',str(kPoints[i]))
		with open('INCAR', 'w') as file:
			file.write(filedata)
		
		time.sleep(1)#Probably not necessary	
	
		process = subprocess.Popen(['qsub','job'], stdout=subprocess.PIPE)
		process.wait()
		print("Job "+str(i + 1) + " Running")


def createINCAR():
	filedata = "Automatic mesh\n0              ! number of k-points = 0 ->automatic generation scheme\nAuto           ! fully automatic\n  1           ! length (l)"
	with open('KPOINTS', 'w') as file:
		file.write(filedata)

def rmINCAR(x):
	os.chdir(x)
	subprocess.Popen(['rm','INCAR'])
	

runTasks()

