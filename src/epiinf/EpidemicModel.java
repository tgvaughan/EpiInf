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

package epiinf;

import beast.core.CalculationNode;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
    
    protected List<EpidemicEvent> eventList;
    protected List<EpidemicState> stateList;
    
    protected Map<EpidemicEvent.EventType, Double> propensities;
    protected double totalPropensity;
    
    EpidemicModel() {
        eventList = Lists.newArrayList();
        stateList = Lists.newArrayList();
        propensities = Maps.newHashMap();
    }
    
    @Override
    public void initAndValidate() { };
    
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
     * Increment state according to reaction of chosen type.
     * 
     * @param state state to update
     * @param type type of reaction to implement
     */
    public abstract void incrementState(EpidemicState state,
            EpidemicEvent.EventType type);

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
     * 
     * @return log of path probability
     */
    public double generateTrajectory(EpidemicState startState,
            double startTime, double endTime) {
        
        double logP = 0.0;
        
        eventList.clear();
        stateList.clear();
        
        EpidemicState thisState = startState.copy();
        double t = startTime;

        while (true) {
            calculatePropensities(thisState);
            
            double dt;
            if (totalPropensity>0.0)
                dt = Randomizer.nextExponential(totalPropensity);
            else
                dt = Double.POSITIVE_INFINITY;
            
            logP += -Math.min(dt, endTime-t)*totalPropensity;
            
            t += dt;
                    
            if (t>=endTime)
                break;
            
            EpidemicEvent nextEvent = new EpidemicEvent();
            nextEvent.time = t;
            
            double u = totalPropensity*Randomizer.nextDouble();
            
            for (EpidemicEvent.EventType type : propensities.keySet()) {
                u -= propensities.get(type);
                
                if (u<0) {
                    nextEvent.type = type;
                    incrementState(thisState, type);
                    logP += Math.log(propensities.get(type));
                    break;
                }
            }
            
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }
        
        return logP;
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
    
    /**
     * Probability density of a path under this model, where the path is
     * defined by the parameters.
     * 
     * @param startTime start time of the path
     * @param endTime end time of the path
     * @param startState state at the beginning of the path
     * @param eventList list of events along the path
     * @return log of probability density
     */
    public double getPathProbability(double startTime, double endTime,
            EpidemicState startState, List<EpidemicEvent> eventList) {
        double logP = 0.0;
        
        double t = startTime;
        EpidemicState thisState = startState.copy();
        
        for (EpidemicEvent event : eventList) {
            
            calculatePropensities(thisState);
            logP += -totalPropensity*(event.time-t)
                    + Math.log(propensities.get(event.type));
            
            incrementState(thisState, event.type);
            t = event.time;
        }
        
        calculatePropensities(thisState);
        if (totalPropensity>0.0)
            logP += -totalPropensity*(endTime-t);
        
        return logP;
    }
    
    /**
     * Probability of a given interval length under model.
     * 
     * @param state state of the system within interval
     * @param waitingTime length of interval
     * @return log of probability density
     */
    public double getIntervalProbability(EpidemicState state,
            double waitingTime) {

        calculatePropensities(state);
        return -totalPropensity*waitingTime;
        
    }
    
    /**
     * Joint probability of an interval followed by an event.
     * 
     * @param startState state at start of interval
     * @param startTime time at start of interval
     * @param epidemicEvent object describing time and type of event
     * @return log of probability density
     */
    public double getJointReactionProbability(EpidemicState startState,
            double startTime, EpidemicEvent epidemicEvent) {
        
        double logP = getIntervalProbability(startState, epidemicEvent.time-startTime);
        logP += Math.log(propensities.get(epidemicEvent.type));
        
        return logP;
    }

}
