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
                
                particleStatesNew.add(particleStates.get(pChoice));
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
        }
        
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
                conditionalP *= model.getProbNoCoalescence(particleState, lineages);
            
            if (eventType == model.getLeafEventType())
                conditionalP *= model.getProbNoLeaf();
        }
        
        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.incrementState(particleState, model.getCoalescenceEventType());
            model.calculatePropensities(particleState);
            conditionalP *= model.getProbCoalescence(particleState, lineages+1)
                    *model.getPropensities().get(model.getCoalescenceEventType());
        } else {
            model.incrementState(particleState, model.getLeafEventType());
            model.calculatePropensities(particleState);
            conditionalP *= model.getProbLeaf()
                    *model.getPropensities().get(model.getLeafEventType());
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
    
}
