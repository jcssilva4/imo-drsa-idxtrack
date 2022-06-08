package rulegen;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import pl.poznan.put.cs.idss.jrs.core.mem.MemoryContainer;

/*
public class StackEntryPoint {

    private Stack stack;

    public StackEntryPoint() {
      stack = new Stack();
      stack.push("Initial Item");
    }

    public Stack getStack() {
        return stack;
    }

    public static void main(String[] args) {
        GatewayServer gatewayServer = new GatewayServer(new StackEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}
*/

import pl.poznan.put.cs.idss.jrs.rules.DomlemAlgorithm;
import pl.poznan.put.cs.idss.jrs.rules.Rule;
import pl.poznan.put.cs.idss.jrs.rules.RulesContainer;
import pl.poznan.put.cs.idss.jrs.utilities.ISFLoader;
import pl.poznan.put.cs.idss.jrs.wrappers.RulesGeneratorWrapper;
import pl.poznan.put.cs.idss.jrs.wrappers.VCdomLEMWrapperOpt;
//import pl.poznan.put.cs.idss.jrs.wrappers.VCdomLEMWrapper;
import py4j.GatewayServer;

//import com.sun.org.apache.xalan.internal.xsltc.runtime.Hashtable;

public class VCdomlemEntryPoint {
    
    public VCdomlemEntryPoint() throws Exception {
		// generate rules
		System.out.println("py4j gateaway server initialized...");
    }

    public List<String> generateRules(String inputPath) throws IOException {
    	List<String> allRules = new LinkedList<String>();
    	/*
    	The sequence of functions used in this code is based on the DomlemAlgorithm.class. 
    	The DomlemAlgorithm.class is the current VC-domlem rule generation code used in jMAF, according to the 
    	pl.poznan.put.cs.idss.jrs_1.0.0/plugin.xml file. 
    	*/
        File inputFile = new File(inputPath);
        //File outputFile = new File(outputPath);
        double confidenceLevel = 1.0D;
		int conditionsSelectionMethod = 1;
		boolean ruleType = false;
		int negativeExamplesTreatment = 0;
		System.out.println("reading inputFile");
        MemoryContainer learningMemoryContainer = ISFLoader.loadISFIntoMemoryContainer(inputFile.getAbsolutePath());
		if (learningMemoryContainer == null) {
			System.out.println("Error occured while reading file '" + inputFile.getAbsolutePath() + "'.");
			//return 2;
			return allRules;
		} else {
			System.out.println("generating rules...");
			Object wrapper;
			wrapper = new VCdomLEMWrapperOpt(learningMemoryContainer, confidenceLevel,
					conditionsSelectionMethod, negativeExamplesTreatment);
			((RulesGeneratorWrapper) wrapper).setInducePossibleRules(ruleType);
			RulesContainer rulesContainer = ((RulesGeneratorWrapper) wrapper).generateRules(learningMemoryContainer);
			System.out.println("saving rules...");
			ArrayList<Rule> certainAtLeast = rulesContainer.getRules(0,3);
			ArrayList<Rule> certainAtMost = rulesContainer.getRules(0,4);
			//String[] allRules = new String[certainAtLeast.size() + certainAtMost.size()];
			//int ruleCounter = 0;
			Iterator rulesIt;
			if (certainAtLeast.size() > 0){
			rulesIt = (Iterator) certainAtLeast.iterator();
				while (rulesIt.hasNext()) {
					Rule rule = (Rule) rulesIt.next();
			    	System.out.println(rule.toShortString());
			    	allRules.add(rule.toShortString());
				}
			}
			if (certainAtMost.size() > 0){
				rulesIt = (Iterator) certainAtMost.iterator();
				while (rulesIt.hasNext()) {
					Rule rule = (Rule) rulesIt.next();
			    	System.out.println(rule.toShortString());
			    	allRules.add(rule.toShortString());
				}
			}
			/*
			System.out.println("writting rules...");
			rulesContainer.getFileInfo().putParameterValue("DataFileDirectory",inputFile.getParentFile().getAbsolutePath());
			rulesContainer.getFileInfo().putParameterValue("DataFileName", inputFile.getName());
			rulesContainer.writeRules(outputFile.getAbsolutePath(), true, true);
			*/
			return allRules;
		}
    }
    

    public static void main(String[] args) throws Exception {
        GatewayServer gatewayServer = new GatewayServer(new VCdomlemEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }
   

}

