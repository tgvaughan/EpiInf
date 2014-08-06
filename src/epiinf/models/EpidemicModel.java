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

package epiinf.models;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Class representing an epidemic model.  Contains all the bits and bobs
 * necessary to calculate path probabilities and generate trajectories under
 * a particular model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EpidemicModel extends CalculationNode {
    
    public Input<Double> psiSamplingProbInput = new Input<>("psiSamplingProb",
            "Probability with which recoveries are translated into samples",
            0.0);
    
    public Input<List<Double>> rhoSamplingProbInput = new Input<>("rhoSamplingProb",
            "Probability with which a lineage at the corresponding time"
                    + "is sampled.", new ArrayList<>());
    
    public Input<List<Double>> rhoSamplingTimeInput = new Input<>("rhoSamplingTime",
            "Times at which rho sampling takes place", new ArrayList<>());
    
    public Input<Double> toleranceInput = new Input<>("tolerance",
            "Maximum absolute time difference between events on tree and "
                    + "events in epidemic trajectory for events to be"
                    + "considered compatible.  Default 1e-10.", 1e-10);
    
    protected List<EpidemicEvent> eventList;
    protected List<EpidemicState> stateList;
    
    protected Map<EpidemicEvent.Type, Double> propensities;
    protected double totalPropensity;
    protected double tolerance;
    
    EpidemicModel() {
        eventList = Lists.newArrayList();
        stateList = Lists.newArrayList();
        propensities = Maps.newHashMap();
    }
    
    @Override
    public void initAndValidate() {
        tolerance = toleranceInput.get();
    };
    
    /**
     * Retrieve tolerance with respect to difference between epidemic
     * and tree event times.
     * 
     * @return tolerance
     */
    public double getTolerance() {
        return tolerance;
    }
    
    /**
     * @return initial state of epidemic
     */
    public abstract EpidemicState getInitialState();
    
    /**
     * Calculate propensities of all possible reactions contained in
     * this model.
     * @param state state used to calculate propensities
     */
    public abstract void calculatePropensities(EpidemicState state);
    
    /**
     * Obtain the most recently calculated reaction propensities.
     * 
     * @return propensity map
     */
    public Map<EpidemicEvent.Type, Double> getPropensities() {
        return propensities;
    }
    
    /**
     * Obtain most recently calculated total reaction propensity.
     * 
     * @return propensity
     */
    public double getTotalPropensity() {
        return totalPropensity;
    }
    
    /**
     * Increment state according to reaction of chosen type.
     * 
     * @param state state to update
     * @param event
     */
    public abstract void incrementState(EpidemicState state,
            EpidemicEvent event);
    

    /**
     * Obtain probability of coalescence occurring on tree given a compatible
     * event occurred in the epidemic trajectory and conditional on the number
     * of extant lineages.
     * 
     * @param state
     * @param lineages
     * @return probability of coalescence event on the tree
     */
    public abstract double getProbCoalescence(EpidemicState state, int lineages);
    
    /**
     * Obtain probability of coalescence NOT occurring on tree given a compatible
     * event occurred in the epidemic trajectory and conditional on the number
     * of extant lineages.
     * 
     * @param state
     * @param lineages
     * @return probability of NO coalescence event on the tree
     */
    public abstract double getProbNoCoalescence(EpidemicState state, int lineages);
    
    /**
     * Retrieve the rho sampling time immediately following t.
     * 
     * @param t
     * @return next rho sampling time (+infinity if there is none)
     */
    public double getNextRhoSamplingTime(double t) {
        for (int i=0; i<rhoSamplingTimeInput.get().size(); i++) {
            if (rhoSamplingTimeInput.get().get(i)>t)
                return rhoSamplingTimeInput.get().get(i);
        }
        
        return Double.POSITIVE_INFINITY;
    }
    
    /**
     * Generate a sequence of events between startTime and endTime conditional
     * on the startState.  The results are retrieved using subsequent calls
     * to getEventList() and getStateList().  Note that the latter only provides
     * the list of states _following_ startState: i.e. there are as many states
     * as events in this list.
     * 
     * @param startState Starting state of trajectory
     * @param startTime Starting time of trajectory
     * @param endTime End time of trajectory
     */
    public void generateTrajectory(EpidemicState startState,
            double startTime, double endTime) {
        
        eventList.clear();
        stateList.clear();
        
        stateList.add(startState);
        
        EpidemicState thisState = startState.copy();
        thisState.time = startTime;

        double nextRhoSamplingTime = Double.POSITIVE_INFINITY;
        int nextRhoSamplingIndex = -1;
        if (rhoSamplingTimeInput.get().isEmpty()) {
            nextRhoSamplingTime = rhoSamplingTimeInput.get().get(0);
            nextRhoSamplingIndex = 0;
        }


        while (true) {
            calculatePropensities(thisState);
            
            double dt;
            if (totalPropensity>0.0)
                dt = Randomizer.nextExponential(totalPropensity);
            else
                dt = Double.POSITIVE_INFINITY;
            
            thisState.time += dt;
                    
            if (thisState.time>=endTime)
                break;
            
            EpidemicEvent nextEvent = new EpidemicEvent();
            
            if (thisState.time > nextRhoSamplingTime) {
                // Simultaneous sampling from the extant population

                nextEvent.type = EpidemicEvent.Type.SAMPLE;

                // Got to be a better way of sampling from a binomial distribution
                nextEvent.multiplicity = 0;
                for (int i=0; i<thisState.I; i++) {
                    if (Randomizer.nextDouble()<rhoSamplingProbInput.get().get(nextRhoSamplingIndex))
                        nextEvent.multiplicity += 1;
                }
                
                nextRhoSamplingIndex += 1;
                
                if (nextRhoSamplingIndex<rhoSamplingTimeInput.get().size())
                    nextRhoSamplingTime = rhoSamplingTimeInput.get().get(nextRhoSamplingIndex);
                else
                    nextRhoSamplingTime = Double.POSITIVE_INFINITY;

            } else {

                nextEvent.time = thisState.time;
            
                double u = totalPropensity*Randomizer.nextDouble();
                
                for (EpidemicEvent.Type type : propensities.keySet()) {
                    u -= propensities.get(type);
                    
                    if (u<0) {
                        nextEvent.type = type;
                        break;
                    }
                }

                // Replace recovery events with sampling events with fixed probability
                if (nextEvent.type == EpidemicEvent.Type.RECOVERY
                        && Randomizer.nextDouble()<psiSamplingProbInput.get())
                    nextEvent.type = EpidemicEvent.Type.SAMPLE;
            }
            
            incrementState(thisState, nextEvent);                
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }
    }
    
    /**
     * @return last simulated event list
     */
    public List<EpidemicEvent> getEventList() {
        return eventList;
    }
    
    /**
     * @return last simulated state list
     */
    public List<EpidemicState> getStateList() {
        return stateList;
    }
}
