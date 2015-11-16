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
import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.ModelEvent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Class representing an epidemic model.  Contains all the bits and bobs
 * necessary to calculate path probabilities and generate trajectories under
 * a particular model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EpidemicModel extends CalculationNode {

    public Input<Function> infectionRateInput = new Input<>(
            "infectionRate", "Infection rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> infectionRateShiftTimesInput = new Input<>(
            "infectionRateShiftTimes", "Infection rate shift times.");

    public Input<Function> recoveryRateInput = new Input<>(
            "recoveryRate", "Recovery rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> recoveryRateShiftTimesInput = new Input<>(
            "recoveryRateShiftTimes", "Recovery rate shift times.");

    public Input<Function> psiSamplingRateInput = new Input<>(
            "psiSamplingRate",
            "Rate at which (destructive) psi-sampling is performed.");

    public Input<RealParameter> psiSamplingRateShiftTimesInput = new Input<>(
            "psiSamplingRateShiftTimes",
            "Times at which psi-sampling rate changes.");

    public Input<Function> removalProbInput = new Input<>(
            "removalProb",
            "Probability that sample individual is removed from population (Default 1).");

    public Input<RealParameter> removalProbShiftTimesInput = new Input<>(
            "removalProbShiftTimes",
            "Times at which removal probability changes.");

    public Input<Function> rhoSamplingProbInput = new Input<>(
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
    protected List<Double[]> rateCache = new ArrayList<>();

    protected boolean ratesDirty;
    protected double tolerance;

    public double[] propensities = new double[EpidemicEvent.nTypes];


    EpidemicModel() { }
    
    @Override
    public void initAndValidate() {
        tolerance = toleranceInput.get();

        if (psiSamplingRateShiftTimesInput.get() != null) {
            if (psiSamplingRateInput.get() == null
                    || psiSamplingRateInput.get().getDimension()
                    != psiSamplingRateShiftTimesInput.get().getDimension() + 1)
                throw new IllegalArgumentException(
                        "Psi sampling rate and rate shift time dimensions " +
                                "don't match.");
        }

        if (removalProbShiftTimesInput.get() != null) {
            if (removalProbInput.get() == null
                    || removalProbInput.get().getDimension()
                    != removalProbShiftTimesInput.get().getDimension() + 1)
                throw new IllegalArgumentException(
                        "Removal prob and and removal prob shift time dimensions " +
                                "don't match.");
        }

        ratesDirty = true;
    }
    
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
        update();
        propensities[EpidemicEvent.RECOVERY] = calculateRecoveryPropensity(state);
        propensities[EpidemicEvent.INFECTION] = calculateInfectionPropensity(state);
        propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] = calculatePsiSamplingRemovePropensity(state);
        propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE] = calculatePsiSamplingNoRemovePropensity(state);
    }

    protected double getCurrentRate(Function rateParam, RealParameter rateShiftTimeParam,
                          double time) {
        if (rateParam != null) {
            if (rateShiftTimeParam != null) {
                int idx = Arrays.binarySearch(
                        rateShiftTimeParam.getValues(), time);
                if (idx<0)
                    idx = -(idx+1);

                return rateParam.getArrayValue(idx);
            } else
                return rateParam.getArrayValue();
        } else
            return 0.0;
    }

    protected double calculatePsiSamplingRemovePropensity(EpidemicState state) {
        return state.I * rateCache.get(state.intervalIdx)[EpidemicEvent.PSI_SAMPLE_REMOVE];
    }
    protected double calculatePsiSamplingNoRemovePropensity(EpidemicState state) {
        return state.I * rateCache.get(state.intervalIdx)[EpidemicEvent.PSI_SAMPLE_NOREMOVE];
    }
    protected abstract double calculateRecoveryPropensity(EpidemicState state);
    protected abstract double calculateInfectionPropensity(EpidemicState state);

    /**
     * Increment state according to reaction of chosen type.
     *
     * @param state state to update
     * @param event epidemic event
     */
    public abstract void incrementState(EpidemicState state,
            EpidemicEvent event);


    /**
     * Update model event list and reaction rate caches.
     */
    public void update() {
        if (!ratesDirty)
            return;

        updateModelEventList();

        if (rateCache.size() == 0) {
            for (int i=0; i< modelEventList.size()+1; i++)
                rateCache.add(new Double[EpidemicEvent.nTypes]);
        }

        for (int i=0; i< modelEventList.size()+1; i++) {
            double t = i>0 ? modelEventList.get(i-1).time : 0.0;

            rateCache.get(i)[EpidemicEvent.INFECTION] = getCurrentRate(
                    infectionRateInput.get(), infectionRateShiftTimesInput.get(), t);
            rateCache.get(i)[EpidemicEvent.RECOVERY] = getCurrentRate(
                    recoveryRateInput.get(), recoveryRateShiftTimesInput.get(), t);

            double psiSamplingRate, removalProb;

            if (psiSamplingRateInput.get() != null)
                psiSamplingRate = getCurrentRate(
                        psiSamplingRateInput.get(), psiSamplingRateShiftTimesInput.get(), t);
            else
                psiSamplingRate = 0.0;

            if (removalProbInput.get() != null)
                removalProb = getCurrentRate(
                        removalProbInput.get(), removalProbShiftTimesInput.get(), t);
            else
                removalProb = 1.0;

            rateCache.get(i)[EpidemicEvent.PSI_SAMPLE_REMOVE] = psiSamplingRate*removalProb;
            rateCache.get(i)[EpidemicEvent.PSI_SAMPLE_NOREMOVE] = psiSamplingRate*(1.0 - removalProb);
        }

        ratesDirty = false;
    }


    /**
     * Assemble list of model events.
     *
     * @return the event list
     */
    public void updateModelEventList() {

        modelEventList.clear();

        if (rhoSamplingProbInput.get() != null) {
            for (int i = 0; i < rhoSamplingProbInput.get().getDimension(); i++) {
                ModelEvent event = new ModelEvent();
                event.type = ModelEvent.Type.RHO_SAMPLING;
                event.rho = rhoSamplingProbInput.get().getArrayValue(i);
                event.time = rhoSamplingTimeInput.get().getValue(i);
                modelEventList.add(event);
            }
        }

        addRateShiftEvents(psiSamplingRateShiftTimesInput.get());
        addRateShiftEvents(removalProbShiftTimesInput.get());
        addRateShiftEvents(infectionRateShiftTimesInput.get());
        addRateShiftEvents(recoveryRateShiftTimesInput.get());

        Collections.sort(modelEventList, (e1, e2) -> {
            if (e1.time < e2.time)
                return -1;

            if (e1.time > e2.time)
                return 1;

            return 0;
        });

    }

    public void addRateShiftEvents(RealParameter rateShiftTimesParam) {

        if (rateShiftTimesParam != null) {
            for (int i=0; i<rateShiftTimesParam.getDimension(); i++) {
                ModelEvent event = new ModelEvent();
                event.type = ModelEvent.Type.RATE_CHANGE;
                event.time = rateShiftTimesParam.getValue(i);
                modelEventList.add(event);
            }
        }
    }

    public double getNextModelEventTime(EpidemicState state) {
        update();

        if (state.intervalIdx<modelEventList.size())
            return modelEventList.get(state.intervalIdx).time;
        else
            return Double.POSITIVE_INFINITY;
    }

    public ModelEvent getNextModelEvent(EpidemicState state) {
        update();

        if (state.intervalIdx < modelEventList.size())
            return modelEventList.get(state.intervalIdx);
        else
            return null;
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
        ratesDirty = true;
        return true;
    }

    @Override
    protected void restore() {
        ratesDirty = true;
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
