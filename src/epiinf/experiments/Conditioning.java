/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package epiinf.experiments;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import epiinf.SimulatedTrajectory;
import epiinf.TransmissionTreeSimulator;
import epiinf.models.SIRModel;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Conditioning {
    
    private static double getTreeLength (Tree tree) {
        double length = 0.0;
        for (Node node : tree.getNodesAsArray())
            length += node.getLength();
        
        return length;
    }
    
    public static void main (String [] args) throws Exception  {
        
        SIRModel model = new SIRModel();
        model.initByName(
                "S0", new IntegerParameter("999"),
                "infectionRate", new RealParameter("0.001"),
                "recoveryRate", new RealParameter("0.2"));
        
        SimulatedTrajectory trajSim = new SimulatedTrajectory();
        
        Tree tree = new Tree();
        RealParameter treeOrigin = new RealParameter();
        
        TransmissionTreeSimulator treeSim = new TransmissionTreeSimulator();

        double[] fracArray = {1.0, 0.5, 0.2, 0.1};
        int nLeaves = 20;
        int nTrees = 5000;

        PrintStream outf = new PrintStream("conditioning.txt");
        outf.println("f discard ratio");
        
        for (int fidx=0; fidx<20; fidx++) {
            double frac = 1.0 - fidx*(1.0-0.1)/(20-1);
            int nDiscard = 0;
            double ratioMean = 0.0;
            for (int i=0; i<nTrees; i++) {
                trajSim.initByName("model", model);
                treeSim.initByName(
                        "tree", tree,
                        "treeOrigin", treeOrigin,
                        "epidemicTrajectory", trajSim,
                        "nLeaves", nLeaves,
                        "model", model,
                        "sampleBeforeTime", frac);
                try {
                    treeSim.initStateNodes();
                    ratioMean += Conditioning.getTreeLength(tree)/nLeaves;
                } catch (Exception ex) {
                    nDiscard += 1;
                }
            }
            
            double discardFrac = nDiscard/(double)nTrees;
            ratioMean /= nTrees-nDiscard;
            
            outf.println(frac + " " + discardFrac + " " + ratioMean);
        }
        
        outf.close();
    }
}
