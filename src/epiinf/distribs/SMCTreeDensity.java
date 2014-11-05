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
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.SimulatedTrajectory;
import epiinf.SimulatedTransmissionTree;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import epiinf.models.EpidemicModel;
import epiinf.models.SISModel;
import epiinf.util.EpiInfUtilityMethods;
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
//    PrintStream debugOut;
    
    public SMCTreeDensity() { }
    
    @Override
    public void initAndValidate() throws Exception {
        model = modelInput.get();
        eventList = treeEventListInput.get();
        nParticles = nParticlesInput.get();
    }
    
    @Override
    public double calculateLogP() throws Exception {
        
        logP = 0.0;
        
        List<Double> particleWeights = Lists.newArrayList();
        List<EpidemicState> particleStates = Lists.newArrayList();
        List<EpidemicState> particleStatesNew = Lists.newArrayList();
        
        // Initialize particles
        for (int p=0; p<nParticles; p++)
            particleStates.add(model.getInitialState());
        
        double t = 0.0;
        int k = 1;
        for (TreeEvent treeEvent : eventList.getEventList()) {
            
            // Update particles
            particleWeights.clear();
            double sumOfWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                
                double newWeight = updateParticle(particleStates.get(p), t, k, treeEvent);
                
                particleWeights.add(newWeight);
                sumOfWeights += newWeight;
            }
            
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
                k -= treeEvent.multiplicity;
            
            // Update start interval time
            t = treeEvent.time;
        }
        
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
        
        while (true) {
            model.calculatePropensities(particleState);
            
            if (model.getTotalPropensity() > 0.0)
                t += Randomizer.nextExponential(model.getTotalPropensity());
            else
                t = Double.POSITIVE_INFINITY;
            
            if (t>finalTreeEvent.time)
                break;
            
            // Check for missed rho sampling event
            if (t > model.getNextRhoSamplingTime(t))
                return 0.0;
            
            double u = model.getTotalPropensity()*Randomizer.nextDouble();
            
            EpidemicEvent event = new EpidemicEvent();
            for (EpidemicEvent.Type type : model.getPropensities().keySet()) {
                u -= model.getPropensities().get(type);
                
                if (u<0) {
                    event.type = type;
                    break;
                }
            }

            switch(event.type) {
                case INFECTION:
                    conditionalP *= 1.0 - lineages*(lineages-1)/particleState.I/(particleState.I+1);
                    break;

                case RECOVERY:
                    // Prob that given recovery is not on sampled lineage, sampling
                    // did not occur.
                    if (model.psiSamplingProbInput.get() != null)
                        conditionalP *= 1.0 - model.psiSamplingProbInput.get().getValue();
                    break;
            }
            
            model.incrementState(particleState, event);
            
            // Early exit if invalid state:
            if (conditionalP==0)
                return 0.0;

            if (particleState.I<lineages)
               return 0.0;
            
        }
        
        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.calculatePropensities(particleState);
            model.incrementState(particleState, EpidemicEvent.Infection);
            conditionalP *= 2.0/particleState.I/(particleState.I-1)
                *model.getPropensities().get(EpidemicEvent.Type.INFECTION);
        } else {
            
            // If the model contains an explicit sampling process, evaluate the
            // probability of the sampling event on the tree
            if (!model.rhoSamplingProbInput.get().isEmpty()
                || model.psiSamplingProbInput.get() != null) {
                
                double sampleProb = 0.0;
                
                // If model contains a rho sampling event at this time, calculate the probability
                // of sampling the number of samples in finalTreeEvent given the current
                // state.
                for (int i=0; i<model.rhoSamplingProbInput.get().size(); i++) {
                    double rhoProb = model.rhoSamplingProbInput.get().get(i).getValue();
                    double rhoTime = model.rhoSamplingTimeInput.get().get(i).getValue();
                    
                    if (Math.abs(rhoTime - finalTreeEvent.time)<model.getTolerance()) {
                        int I = (int)Math.round(particleState.I);
                        int k = finalTreeEvent.multiplicity;
                        sampleProb += Binomial.choose(I, k)
                            *Math.pow(rhoProb, k)*Math.pow(1.0-rhoProb, I-k);
                    }
                }
                
                // If the model contains a non-zero psi sampling rate, calculate the
                // probability of a sampled recovery occuring at the time of finalTreeEvent
                
                if (finalTreeEvent.multiplicity==1
                    && model.psiSamplingProbInput.get() != null) {
                    model.calculatePropensities(particleState);
                    sampleProb += model.getPropensities().get(EpidemicEvent.Type.RECOVERY)
                        *model.psiSamplingProbInput.get().getValue();
                }
                
                conditionalP *= sampleProb*Math.exp(GammaFunction.lnGamma(1+finalTreeEvent.multiplicity));
                
            }
            
            model.incrementState(particleState,
                EpidemicEvent.MultipleSamples(finalTreeEvent.multiplicity));
            
        }
        
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

//        Randomizer.setSeed(5);
        
        EpidemicModel model;

        model = new SISModel();
        model.initByName(
            "S0", new IntegerParameter("99"),
            "infectionRate", new RealParameter("0.01"),
            "recoveryRate", new RealParameter("0.2"),
            "rhoSamplingProb", new RealParameter("0.3"),
            "rhoSamplingTime", new RealParameter("4.0"));

        EpidemicTrajectory traj;
        Tree tree;

        do {
            do {
                traj = new SimulatedTrajectory();
                traj.initByName(
                    "model", model,
                    "maxDuration", 5.0);
            } while (!traj.hasSample());

            tree = new SimulatedTransmissionTree();
            tree.initByName(
                "epidemicTrajectory", traj,
                "fileName", "tree.newick");
        } while (false); //tree.getLeafNodeCount() != 10);
        
        RealParameter treeOrigin = new RealParameter("4.0");

        try (PrintStream ps = new PrintStream("expotree.txt")) {
            EpiInfUtilityMethods.writeExpoTreeFile(tree, treeOrigin.getValue(), ps);
        }
        
        TreeEventList treeEventList = new TreeEventList();
        treeEventList.initByName(
            "tree", tree,
            "treeOrigin", treeOrigin);
        
        SMCTreeDensity treeDensity = new SMCTreeDensity();
        
        try (PrintStream ps = new PrintStream("logLik.txt")) {
            ps.println("beta logP");
            for (double beta=0.005; beta<0.015; beta += 0.0005) {
                
                model = new SISModel();
                model.initByName(
                    "S0", new IntegerParameter("99"),
                    "infectionRate", new RealParameter(String.valueOf(beta)),
                    "recoveryRate", new RealParameter("0.2"),
                    "rhoSamplingProb", new RealParameter("0.3"),
                    "rhoSamplingTime", new RealParameter("4.0"));
                
                treeDensity.initByName(
                    "treeEventList", treeEventList,
                    "model", model,
                    "nParticles", 10000);
                
                double logP = treeDensity.calculateLogP();
                
                System.out.println("beta: " + beta
                    + " logP: " + logP);
                ps.println(beta + " " + logP);
            }
        }
    }
}
