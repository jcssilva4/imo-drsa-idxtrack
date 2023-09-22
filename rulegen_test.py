from py4j.java_gateway import JavaGateway
from py4j.protocol import Py4JError

import os

def get_rules(gateway):
	cwd = os.getcwd()
	inputFile = cwd + "/rulegen/data/currentDT.isf"
	# induce rules
	ggclass = gateway.jvm.rulegen.VCdomlemEntryPoint()

	# get the rules
	allRules = []
	try:
		allRules = ggclass.generateRules(inputFile)
	except Py4JError as e:
		print("Error...the rule set was not induced....trying again...")
	#'''


	#print(str(len(allRules)) + " Rules were generated: ")
	#print(allRules)
	return allRules

def get_rules_old(gateway):
	cwd = os.getcwd()
	inputFile = cwd + "/rulegen/data/currentDT.isf"
	# induce rules
	ggclass = gateway.jvm.rulegen.VCdomlemEntryPoint()

	# get the rules
	allRules = ggclass.generateRules(inputFile)
	#'''


	#print(str(len(allRules)) + " Rules were generated: ")
	#print(allRules)
	return allRules


def launch_gatewayserver():
	cwd = os.getcwd()
	# set the input and output test files
	inputFile = cwd + "/rulegen/data/currentDT.isf"
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
	return gateway

def debbug_gatewayserver():
	cwd = os.getcwd()
	# set the input and output test files
	inputFile = cwd + "/rulegen/data/currentDT.isf"

	## debbuging code (you can use this code to debbug the methods of your java class) ##
	# run the gateway server in eclipse before executing this python script (https://www.py4j.org/getting_started.html)
	# get the gateway server
	gateway = JavaGateway() 
	allRules = []
	allRules = gateway.entry_point.generateRules(inputFile)
	'''
	while len(allRules) == 0:
		try:
			allRules = gateway.entry_point.generateRules(inputFile)
		except Py4JError as e:
			#if e.cause:
			#	print(e.cause)
			print("currentRuleVar: " + str(allRules))

	print("final rulevar: " + str(allRules))	
	'''
	return allRules
	
