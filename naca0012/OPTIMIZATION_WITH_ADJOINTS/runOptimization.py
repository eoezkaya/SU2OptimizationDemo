import os

#print("copying files for the training data...\n")
os.system("cp ./TrainingData/CD.csv ./")
os.system("cp ./TrainingData/CL.csv ./")
os.system("cp ./TrainingData/Area.csv ./")
os.system("./preprocessTrainingData.py")


f = open("./SimulationData/iterationNumber.dat", "w")
f.write("1")
f.close()

RODEO_HOME = "/home/eoezkaya/RoDeO"
BIN_RODEO = RODEO_HOME + "/build/rodeo"
configFilename = "NACA0012DragMinimization.cfg"
COMMAND = BIN_RODEO + " " +  configFilename

os.system(COMMAND)

