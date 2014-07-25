/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package epiinf.distribs;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.TrajectorySimulator;
import epiinf.TransmissionTreeSimulator;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import epiinf.models.EpidemicModel;
import epiinf.models.SIRModel;
import epiinf.models.SISModel;
import java.io.PrintStream;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
        + "parameters.")
public class SMCTreeDensity extends Distribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<TreeEventList> treeEventListInput = new Input<>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    EpidemicModel model;
    TreeEventList eventList;
    int nParticles;
    
    // DEBUG
    PrintStream debugOut;

    public SMCTreeDensity() { }

    @Override
    public void initAndValidate() throws Exception {
        model = modelInput.get();
        eventList = treeEventListInput.get();
        nParticles = nParticlesInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        
        // DEBUG
        debugOut = new PrintStream("SMCdebug.json");
        debugOut.println("{");
        
        logP = 0.0;
        
        List<Double> particleWeights = Lists.newArrayList();
        List<EpidemicState> particleStates = Lists.newArrayList();
        List<EpidemicState> particleStatesNew = Lists.newArrayList();
        
        // Initialize particles
        for (int p=0; p<nParticles; p++)
            particleStates.add(model.getInitialState());
        
        double t = 0.0;
        int k = 1;
        int interval = 0;
        for (TreeEvent treeEvent : eventList.getEventList()) {
            
            // DEBUG
            if (interval>0)
                debugOut.println(",");
            debugOut.println("\"interval" + interval + "\": {");
            
            // Update particles
            particleWeights.clear();
            double sumOfWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                
                // DEBUG
                if (p>0)
                    debugOut.println(", ");
                debugOut.println("\"p" + p + "\": {");
                
                double newWeight = updateParticle(particleStates.get(p), t, k, treeEvent);
                
                particleWeights.add(newWeight);                
                sumOfWeights += newWeight;
                
                // DEBUG
                debugOut.print("\n}");
            }
            
            // DEBUG
            debugOut.print("}");
            
            // Update marginal likelihood estimate
            logP += Math.log(sumOfWeights/nParticles);
            
            if (!(sumOfWeights>0.0))
                return Double.NEGATIVE_INFINITY;

            // Sample particle with replacement
            particleStatesNew.clear();
            for (int p=0; p<nParticles; p++) {
                double u = Randomizer.nextDouble()*sumOfWeights;
                
                int pChoice;
                for (pChoice = 0; pChoice<nParticles; pChoice++) {
                    u -= particleWeights.get(pChoice);
                    if (u<0.0)
                        break;
                }
                
                if (pChoice == nParticles)
                    System.err.println("sumOfWeights: " + sumOfWeights);
                
                particleStatesNew.add(particleStates.get(pChoice).copy());
            }
            
            // Switch particleStates and particleStatesNew
            List<EpidemicState> temp = particleStates;
            particleStates = particleStatesNew;
            particleStatesNew = temp;
            
            // Update lineage counter
            if (treeEvent.type == TreeEvent.Type.COALESCENCE)
                k += 1;
            else
                k -= 1;
            
            // Update start interval time
            t = treeEvent.time;
            
            interval += 1;
        }
        
        // DEBUG
        debugOut.println("}");
        
        double sumOfWeights = 0.0;
        for (double weight : particleWeights)
            sumOfWeights += weight;
        
        logP = Math.log(sumOfWeights/nParticles);
        return logP;
    }

    /**
     * Updates weight and state of particle.
     * 
     * @param particleState
     * @param startTime
     * @param lineages
     * @param finalTreeEvent
     * @return conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
            double startTime, int lineages, TreeEvent finalTreeEvent) {
        double conditionalP = 1.0;
        
        double t = startTime;
        
        // DEBUG
        List<Double> tList = Lists.newArrayList();
        List<Double> nList = Lists.newArrayList();
        tList.add(t);
        nList.add(particleState.I);
        
        while (true) {
            model.calculatePropensities(particleState);

            if (model.getTotalPropensity() > 0.0)
                t += Randomizer.nextExponential(model.getTotalPropensity());
            else
                t = Double.POSITIVE_INFINITY;
            
            if (t>finalTreeEvent.time)
                break;
            
            double u = model.getTotalPropensity()*Randomizer.nextDouble();
            
            EpidemicEvent.Type eventType = null;
            for (EpidemicEvent.Type type : model.getPropensities().keySet()) {
                u -= model.getPropensities().get(type);
                
                if (u<0) {
                    eventType = type;
                    break;
                }
            }
            
            model.incrementState(particleState, eventType);
            
            // Early exit if invalid state:
            if (particleState.I<lineages)
                return 0.0;
            
            // Increment conditional prob
            if (eventType == EpidemicEvent.Type.INFECTION)
                conditionalP *= model.getProbNoCoalescence(particleState, lineages);
            
            if (eventType == EpidemicEvent.Type.SAMPLE)
                return 0.0;
            
            // DEBUG
            tList.add(t);
            nList.add(particleState.I);
            
        }
        
        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.incrementState(particleState, EpidemicEvent.Type.INFECTION);
            model.calculatePropensities(particleState);
            conditionalP *= model.getProbCoalescence(particleState, lineages+1)
                    *model.getPropensities().get(EpidemicEvent.Type.INFECTION);
        } else {
            model.incrementState(particleState, EpidemicEvent.Type.SAMPLE);
            model.calculatePropensities(particleState);
            if (model.getPropensities().containsKey(EpidemicEvent.Type.SAMPLE))
                conditionalP *= model.getPropensities().get(EpidemicEvent.Type.SAMPLE);
        }
        
        // DEBUG
        tList.add(finalTreeEvent.time);
        nList.add(particleState.I);
        debugOut.print("\"t\": [");
        for (int s=0; s<tList.size(); s++) {
            if (s>0)
                debugOut.print(",");
            debugOut.print(tList.get(s));
        }
        debugOut.print("], \"n\": [");
        for (int s=0; s<nList.size(); s++) {
            if (s>0)
                debugOut.print(",");
            debugOut.print(nList.get(s));
        }
        debugOut.print("]");
        
        if (!particleState.isValid())
            return 0.0;
        else
            return conditionalP;
    }
    
    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
    
    /**
     * Main method for testing/debugging.
     * 
     * @param args 
     * @throws java.lang.Exception 
     */
    public static void main (String [] args) throws Exception {
        
        
        SISModel model = new SISModel();
        model.initByName(
                "S0", new IntegerParameter("999"),
                "infectionRate", new RealParameter("0.001"),
                "recoveryRate", new RealParameter("0.2"));
        
        TrajectorySimulator trajSim = new TrajectorySimulator();
        trajSim.initByName(
                "model", model,
                "maxDuration", 50.0,
                "samplingTime", 12.0,
                "sampleSize", 100,
                "fileName", "truth.txt");
        
        Tree tree = new Tree();
        RealParameter treeOrigin = new RealParameter();
        
        TransmissionTreeSimulator treeSim = new TransmissionTreeSimulator();
        treeSim.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin,
                "epidemicTrajectory", trajSim,
                "model", model,
                "fileName", "truth.newick");
        treeSim.initStateNodes();
        
        TreeEventList treeEventList = new TreeEventList();
        treeEventList.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin); 
        
        SMCTreeDensity treeDensity = new SMCTreeDensity();
        treeDensity.initByName(
                "treeEventList", treeEventList,
                "model", model,
                "nParticles", 100);
        
//        for (int bidx=-5; bidx<=5; bidx++) {
//            double beta = 0.001*Math.pow(1.3, bidx);
        double beta = 0.001;
            
            model.initByName(
                "S0", new IntegerParameter("999"),
                "infectionRate", new RealParameter(String.valueOf(beta)),
                "recoveryRate", new RealParameter("0.2"));
            
            treeDensity.initByName(
                "treeEventList", treeEventList,
                "model", model,
                "nParticles", 100);
            
            System.out.println("beta: " + beta
                    + " logP: " + treeDensity.calculateLogP());
//        }
    }
}
