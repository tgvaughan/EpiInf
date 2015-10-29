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
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.ModelEvent;

import java.util.ArrayList;
import java.util.HashMap;
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
    
    public Input<RealParameter> psiSamplingRateInput = new Input<>(
            "psiSamplingRate",
            "Rate at which (destructive) psi-sampling is performed.");

    public Input<RealParameter> rhoSamplingProbInput = new Input<>(
            "rhoSamplingProb",
            "Probability with which a lineage at the corresponding time"
                    + "is sampled.");
    
    public Input<RealParameter> rhoSamplingTimeInput = new Input<>(
            "rhoSamplingTime",
            "Times at which rho sampling takes place");
    
    public Input<Double> toleranceInput = new Input<>("tolerance",
            "Maximum absolute time difference between events on tree and "
                    + "events in epidemic trajectory for events to be"
                    + "considered compatible.  Default 1e-10.", 1e-10);
    
    protected List<EpidemicEvent> eventList = new ArrayList<>();
    protected List<EpidemicState> stateList = new ArrayList<>();
    protected List<ModelEvent> modelEventList = new ArrayList<>();
    protected List<Double> modelEventTimes = new ArrayList<>();

    protected boolean modelEventListDirty;
    protected double tolerance;

    public Map<EpidemicEvent.Type,Double> propensities = new HashMap<>();


    EpidemicModel() { }
    
    @Override
    public void initAndValidate() {
        tolerance = toleranceInput.get();

        modelEventListDirty = true;
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
    

    public final void calculatePropensities(EpidemicState state) {
        propensities.put(EpidemicEvent.Type.RECOVERY, calculateRecoveryPropensity(state));
        propensities.put(EpidemicEvent.Type.INFECTION, calculateInfectionPropensity(state));

        if (psiSamplingRateInput.get() != null)
            propensities.put(EpidemicEvent.Type.PSI_SAMPLE, state.I * psiSamplingRateInput.get().getValue());
        else
            propensities.put(EpidemicEvent.Type.PSI_SAMPLE, 0.0);
    }

    public abstract double calculateRecoveryPropensity(EpidemicState state);
    public abstract double calculateInfectionPropensity(EpidemicState state);

    
    /**
     * Increment state according to reaction of chosen type.
     * 
     * @param state state to update
     * @param event
     */
    public abstract void incrementState(EpidemicState state,
            EpidemicEvent event);


    /**
     * Assemble list of model events.
     *
     * @return the event list
     */
    public List<ModelEvent> getModelEventList() {

        if (modelEventListDirty) {
            modelEventList.clear();

            if (rhoSamplingProbInput.get() != null) {
                for (int i = 0; i < rhoSamplingProbInput.get().getDimension(); i++) {
                    ModelEvent event = new ModelEvent();
                    event.type = ModelEvent.Type.RHO_SAMPLING;
                    event.rho = rhoSamplingProbInput.get().getValue(i);
                    event.time = rhoSamplingTimeInput.get().getValue(i);
                    modelEventList.add(event);
                    modelEventTimes.add(event.time);
                }
            }

            modelEventListDirty = false;
        }

        return modelEventList;
    }

    /**
     * Generate a sequence of events between startTime and endTime conditional
     * on the startState.  The results are retrieved using subsequent calls
     * to getEpidemicEventList() and getEpidemicStateList().  Note that the latter only provides
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

        List<ModelEvent> theseModelEvents = new ArrayList<>();
        getModelEventList().stream()
                .filter(e -> e.time > startTime && e.time < endTime)
                .forEach(theseModelEvents::add);

        while (true) {
            calculatePropensities(thisState);

            double totalPropensity = propensities.get(EpidemicEvent.Type.INFECTION)
                    + propensities.get(EpidemicEvent.Type.RECOVERY)
                    + propensities.get(EpidemicEvent.Type.PSI_SAMPLE);

            double dt;
            if (totalPropensity>0.0)
                dt = Randomizer.nextExponential(totalPropensity);
            else
                dt = Double.POSITIVE_INFINITY;
            
            thisState.time += dt;

            EpidemicEvent nextEvent = new EpidemicEvent();

            if (!theseModelEvents.isEmpty() && thisState.time > theseModelEvents.get(0).time) {
                ModelEvent event = theseModelEvents.get(0);
                theseModelEvents.remove(0);
                if (event.type == ModelEvent.Type.RHO_SAMPLING) {
                    nextEvent.type = EpidemicEvent.Type.RHO_SAMPLE;

                    // Got to be a better way of sampling from a binomial distribution
                    nextEvent.multiplicity = 0;
                    for (int i = 0; i < thisState.I; i++) {
                        if (Randomizer.nextDouble() < event.rho)
                            nextEvent.multiplicity += 1;
                    }


                    nextEvent.time = event.time;
                    thisState.time = event.time;

                    incrementState(thisState, nextEvent);
                    eventList.add(nextEvent);
                    stateList.add(thisState.copy());
                }

                continue;
            }

            if (thisState.time>=endTime)
                break;


            nextEvent.time = thisState.time;

            double u = totalPropensity*Randomizer.nextDouble();

            for (EpidemicEvent.Type type : propensities.keySet()) {
                u -= propensities.get(type);

                if (u<0) {
                    nextEvent.type = type;
                    break;
                }
            }

            incrementState(thisState, nextEvent);
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }
    }
    
    /**
     * @return last simulated event list
     */
    public List<EpidemicEvent> getEpidemicEventList() {
        return eventList;
    }
    
    /**
     * @return last simulated state list
     */
    public List<EpidemicState> getEpidemicStateList() {
        return stateList;
    }

    @Override
    protected boolean requiresRecalculation() {
        modelEventListDirty = true;
        return true;
    }

    @Override
    protected void restore() {
        modelEventListDirty = true;
    }

    /**
     * Test whether two times are equal to within the model's tolerance.
     *
     * @param timeA first time
     * @param timeB second time
     * @return true if times are to be considered equal
     */
    public boolean timesEqual(double timeA, double timeB) {
        return Math.abs(timeA - timeB) < tolerance;
    }

}
