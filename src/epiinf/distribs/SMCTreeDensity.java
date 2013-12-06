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
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import epiinf.models.EpidemicModel;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
        + "parameters.")
public class SMCTreeDensity extends Distribution {

    public Input<EpidemicModel> modelInput = new Input<EpidemicModel>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<Integer> nParticlesInput = new Input<Integer>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    EpidemicModel model;
    TreeEventList eventList;
    int nParticles;

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
        List<Double> particleWeightsNew = Lists.newArrayList();
        List<EpidemicState> particleStates = Lists.newArrayList();
        List<EpidemicState> particleStatesNew = Lists.newArrayList();
        
        // Initialize particles
        for (int p=0; p<nParticles; p++) {
            particleWeights.add(1.0);
            particleStates.add(model.getInitialState());
        }
        
        double t = 0.0;
        int k = 1;
        for (TreeEvent treeEvent : eventList.getEventList()) {
            
            // Update particles
            double sumOfWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                double newWeight = particleWeights.get(p) + 
                        updateParticle(particleStates.get(p), t, k, treeEvent);
                
                particleWeights.set(p, newWeight);                
                sumOfWeights += newWeight;
            }

            // Sample particle with replacement
            particleWeightsNew.clear();
            particleStatesNew.clear();
            for (int p=0; p<nParticles; p++) {
                double u = Randomizer.nextDouble()*sumOfWeights;
                
                int pChoice;
                for (pChoice = 0; pChoice<nParticles; pChoice++) {
                    u -= particleWeights.get(pChoice);
                    if (u<0.0)
                        break;
                }
                
                particleWeightsNew.add(1.0);
                particleStatesNew.add(particleStates.get(pChoice));
            }
            
            // Update lineage counter
            k += 1;
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
     * @return 
     */
    private double updateParticle(EpidemicState particleState,
            double startTime, int lineages, TreeEvent finalTreeEvent) {
        double conditionalLogP = 0.0;
        
        double t = startTime;
        
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
                    model.incrementState(particleState, type);
                    eventType = type;
                    break;
                }
            }
            
            // Increment conditional prob
            if (eventType == model.getCoalescenceEventType())
                conditionalLogP += Math.log(model.getProbNoCoalescence(particleState, lineages));
            
            if (eventType == model.getLeafEventType())
                conditionalLogP += Math.log(model.getProbNoLeaf());
        }
        
        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.incrementState(particleState, model.getCoalescenceEventType());
            conditionalLogP += Math.log(model.getProbCoalescence(particleState, lineages+1));
        } else {
            model.incrementState(particleState, model.getLeafEventType());
            conditionalLogP += Math.log(model.getProbLeaf());
        }
        
        return conditionalLogP;
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
    
}
