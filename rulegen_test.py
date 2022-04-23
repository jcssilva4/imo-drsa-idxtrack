from py4j.java_gateway import JavaGateway
import os


cwd = os.getcwd()

# set the input and output test files
inputFile = cwd + "/rulegen/data/drsa_countryRisk_2014.isf"
#outputFile = cwd + "/rulegen/data/drsa_countryRisk2.rules"

'''
## debbuging code (you can use this code to debbug the methods of your java class) ##
# run the gateway server in eclipse before executing this python script (https://www.py4j.org/getting_started.html)
# get the gateway server
gateway = JavaGateway() 

allRules = gateway.entry_point.generateRules(inputFile, outputFile)
'''


#'''
## production code (you can execute this after you generate the rulegen.jar) ##

# set the gateway server
# .jar file containing the VC-Domlem class
classpath = cwd + "/rulegen/rulegen.jar"
# dependencies
plg1 = cwd + "/rulegen/plugins/ForemkaCore_3.1.0.jar"
plg2 = cwd + "/rulegen/plugins/jRS.jar"
plg3 = cwd + "/rulegen/plugins/py4j0.10.9.5.jar"
plg4 = cwd + "/rulegen/plugins/trove-2.0.1 rc1.jar"
plg5 = cwd + "/rulegen/plugins/xsltc.jar"
classpath = os.pathsep.join((classpath, plg1, plg2, plg3, plg4, plg5))
# gateway server
gateway = JavaGateway.launch_gateway(classpath=classpath)
ggclass = gateway.jvm.rulegen.VCdomlemEntryPoint()

# get the rules
allRules = ggclass.generateRules(inputFile)
#'''

print(str(len(allRules)) + " Rules were generated: ")
print(allRules)
