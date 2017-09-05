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
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.ModelEvent;

import java.util.*;

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

    public Input<Function> infectionRateShiftTimesInput = new Input<>(
            "infectionRateShiftTimes", "Infection rate shift times.");

    public Input<Boolean> infectionRateShiftTimesBackwardInput = new Input<>(
            "infectionRateShiftTimesBackward",
            "If true, infection rate shift times are assumed to be expressed" +
                    " as times before most recent sample.",
            false);


    public Input<Function> recoveryRateInput = new Input<>(
            "recoveryRate", "Recovery rate.", Input.Validate.REQUIRED);

    public Input<Function> recoveryRateShiftTimesInput = new Input<>(
            "recoveryRateShiftTimes", "Recovery rate shift times.");

    public Input<Boolean> recoveryRateShiftTimesBackwardInput = new Input<>(
            "recoveryRateShiftTimesBackward",
            "If true, recovery rate shift times are assumed to be expressed" +
                    " as times before most recent sample.",
            false);


    public Input<Function> psiSamplingVariableInput = new Input<>(
            "psiSamplingVariable",
            "Represents either the rate of psi-sampling or the " +
                    "proportion of individuals that generate at least " +
                    "one psi-sample while infectious, depending on " +
                    "the value of useSamplingProportion.");

    public Input<Function> psiSamplingVariableShiftTimesInput = new Input<>(
            "psiSamplingVariableShiftTimes",
            "Times at which psi-sampling variable changes.");

    public Input<Boolean> psiSamplingVariableShiftTimesBackwardInput = new Input<>(
            "psiSamplingVariableShiftTimesBackward",
            "If true, psi sampling variable shift times are assumed to be expressed" +
                    " as times before most recent sample.",
            false);

    public Input<Boolean> usePsiSamplingProportionInput = new Input<>(
            "usePsiSamplingProportion",
            "If true, psi sampling variable represents proportion psi/(mu+psi)" +
                    "of individuals that are sampled while infectious, " +
                    "otherwise it represents the per-indidivual rate of " +
                    "psi sampling",
            false);


    public Input<Function> removalProbInput = new Input<>(
            "removalProb",
            "Probability that sampled individual is removed from population.",
            Input.Validate.REQUIRED);

    public Input<Function> removalProbShiftTimesInput = new Input<>(
            "removalProbShiftTimes",
            "Times at which removal probability changes.");

    public Input<Boolean> removalProbShiftTimesBackwardInput = new Input<>(
            "removalProbShiftTimesBackward",
            "If true, removal probability change times are assumed to be expressed" +
                    " as times before most recent sample.",
            false);


    public Input<Function> rhoSamplingProbInput = new Input<>(
            "rhoSamplingProb",
            "Probability with which a lineage at the corresponding time"
                    + "is sampled.");
    
    public Input<Function> rhoSamplingTimeInput = new Input<>(
            "rhoSamplingTime",
            "Times at which rho sampling takes place");

    public Input<Boolean> rhoSamplingTimesBackwardInput = new Input<>(
            "rhoSamplingTimesBackward",
            "If true, rho sampling times are assumed to be expressed" +
                    " as times before most recent sample.",
            false);


    public Input<Function> originInput = new Input<>(
            "origin",
            "Time before end of sampling process that epidemic began.",
            Input.Validate.REQUIRED);

    public Input<Double> toleranceInput = new Input<>("tolerance",
            "Maximum absolute time difference between events on tree and "
                    + "events in epidemic trajectory for events to be"
                    + "considered compatible.  Default 1e-10.", 1e-10);
    
    protected List<EpidemicEvent> eventList = new ArrayList<>();
    protected List<EpidemicState> stateList = new ArrayList<>();
    protected List<ModelEvent> modelEventList = new ArrayList<>();
    protected List<Double[]> rateCache = new ArrayList<>();
    protected double[] removalProbCache = null;

    protected boolean ratesDirty;
    protected double tolerance;

    public double[] propensities = new double[EpidemicEvent.nTypes];
    public double currentRemovalProb;

    EpidemicModel() { }
    
    @Override
    public void initAndValidate() {
        tolerance = toleranceInput.get();

        if (psiSamplingVariableShiftTimesInput.get() != null) {
            if (psiSamplingVariableInput.get() == null
                    || psiSamplingVariableInput.get().getDimension()
                    != psiSamplingVariableShiftTimesInput.get().getDimension() + 1)
                throw new IllegalArgumentException(
                        "Psi sampling variable and variable shift time dimensions " +
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
     * @return initial state of epidemic
     */
    public abstract EpidemicState getInitialState();

    public final void calculatePropensities(EpidemicState state) {
        update();
        propensities[EpidemicEvent.RECOVERY] = calculateRecoveryPropensity(state);
        propensities[EpidemicEvent.INFECTION] = calculateInfectionPropensity(state);
        propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] = calculatePsiSamplingRemovePropensity(state);
        propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE] = calculatePsiSamplingNoRemovePropensity(state);
        currentRemovalProb = calculateCurrentRemovalProb(state);
    }

    /**
     * @return age of epidemic start relative to last tree event
     */
    public double getOrigin() {
        if (originInput.get() == null)
            throw new UnsupportedOperationException(
                    "Origin requested from EpidemicModel not initialized " +
                    "with an origin parameter.");

        return originInput.get().getArrayValue();
    }

    /**
     * Convert time parameter element to forwards direction if originally
     * expressed in backwards direction, also reversing the element order.
     * This ensures that if elements are in order of increasing backwards
     * time the return values of this function call will be in order of
     * increasing forwards time.
     *
     * @param timeParam time parameter
     * @param i index into parameter
     * @param isReversed if true, original times are backwards
     * @return forward time
     */
    protected double getForwardTime(Function timeParam, int i, boolean isReversed) {
        if (isReversed) {
            return originInput.get().getArrayValue()
                    - timeParam.getArrayValue(timeParam.getDimension()-1-i);
        } else {
            return timeParam.getArrayValue(i);
        }
    }

    /**
     * Obtain ith element of rate parameter where the index i counts the elements
     * forward in time.  If the rate parameter is expressed in the backward direction,
     * this means that the order of elements in the parameter will be reversed,
     * so that earlier time intervals come before later intervals.
     *
     * @param rateParam rate parameter
     * @param i index into parameter
     * @param isReversed if true, order is reversed
     * @return ith element, with order reversal in the case that isReversed is true
     */
    protected double getRateInInterval(Function rateParam, int i, boolean isReversed) {
        if (isReversed)
            return rateParam.getArrayValue(rateParam.getDimension()-i-1);
        else
            return rateParam.getArrayValue();
    }

    /**
     * Retrieve rate at given time.
     *
     * @param rateParam parameter defining rate
     * @param rateShiftTimeParam parameter defining rate change times
     * @param paramTimesBackwards if true, rate shift times are backwards
     * @param time time at which to evaluate rate.
     *
     * @return rate effective at given time
     */
    protected double getRateAtTime(Function rateParam, Function rateShiftTimeParam,
                                   boolean paramTimesBackwards, double time) {
        if (rateParam != null) {
            if (rateShiftTimeParam != null) {
                return rateParam.getArrayValue(
                        binarySearch(rateShiftTimeParam, paramTimesBackwards, time));
            } else
                return rateParam.getArrayValue();
        } else
            return 0.0;
    }

    /**
     * Use binary search to identify interval index corresponding to given time.
     *
     * @param rateShiftTimeParam function identifying rate shift times
     * @param paramTimesBackwards if true, rate shift times are backwards
     * @param time time to place
     * @return index into corresponding interval
     */
    protected int binarySearch(Function rateShiftTimeParam,
                                      boolean paramTimesBackwards, double time) {

        int N = rateShiftTimeParam.getDimension()+1;
        int imin=0, imax=N-1;

        while (imax>imin) {
            int imid = imin + (imax-imin)/2;

            if (imid==N-1 || getForwardTime(rateShiftTimeParam, imid, paramTimesBackwards)>time) {
                if (imid==0 || getForwardTime(rateShiftTimeParam, imid-1, paramTimesBackwards)<=time)
                    if (paramTimesBackwards)
                        return N-1 - imin;
                    else
                        return imid;
                else
                    imax = imid;
            } else {
                imin = imid+1;
            }
        }

        if (paramTimesBackwards)
            return N-1 - imin;
        else
            return imin;
    }

    protected double calculatePsiSamplingRemovePropensity(EpidemicState state) {
        return state.I * rateCache.get(state.modelIntervalIdx)[EpidemicEvent.PSI_SAMPLE_REMOVE];
    }
    protected double calculatePsiSamplingNoRemovePropensity(EpidemicState state) {
        return state.I * rateCache.get(state.modelIntervalIdx)[EpidemicEvent.PSI_SAMPLE_NOREMOVE];
    }
    protected abstract double calculateRecoveryPropensity(EpidemicState state);
    protected abstract double calculateInfectionPropensity(EpidemicState state);

    protected double calculateCurrentRemovalProb(EpidemicState state) {
        return removalProbCache[state.modelIntervalIdx];
    }

    /**
     * Increment state according to reaction of chosen type.
     *
     * @param state state to update
     * @param event epidemic event
     */
    public abstract void incrementState(EpidemicState state,
            EpidemicEvent event);

    /**
     * Uses algorithm outlined in Cao et al. (JCP, 2006) to select the next
     * tau leaping step size.
     *
     * @param epsilon relative change in propensity to allow
     * @param state epidemic state
     * @param infectProp infection propensity to use
     * @param recovProp recovery propensity to use
     * @return selected tau
     */
    public abstract double getTau(double epsilon, EpidemicState state,
                                  double infectProp, double recovProp);

    /**
     * Transform various inference parameterizations into a uniform
     * simulation parameterization.
     *
     * @param transformedRates Array where simulation parameters are recorded.
     * @param currentRates Map where inference parameters are recorded.
     */
    protected void transformRates(Double[] transformedRates, Map<ModelEvent.RateVariableType, Double> currentRates) {

        transformedRates[EpidemicEvent.INFECTION] = currentRates.get(ModelEvent.RateVariableType.INFECTION_RATE);
        transformedRates[EpidemicEvent.RECOVERY] = currentRates.get(ModelEvent.RateVariableType.RECOVERY_RATE);

        double psiSamplingVariable, psiSamplingRate, removalProb;
        psiSamplingVariable = currentRates.get(ModelEvent.RateVariableType.PSI_SAMPLING_VARIABLE);
        if (usePsiSamplingProportionInput.get()) {
            if (psiSamplingVariable > 0.0) {
                psiSamplingRate = currentRates.get(ModelEvent.RateVariableType.RECOVERY_RATE)
                        / (1.0 / psiSamplingVariable - 1.0);
            } else {
                psiSamplingRate = 0.0;
            }
        } else {
            psiSamplingRate = psiSamplingVariable;
        }

        removalProb = currentRates.get(ModelEvent.RateVariableType.REMOVAL_PROB);

        transformedRates[EpidemicEvent.PSI_SAMPLE_REMOVE] = psiSamplingRate*removalProb;
        transformedRates[EpidemicEvent.PSI_SAMPLE_NOREMOVE] = psiSamplingRate*(1.0-removalProb);
    }

    /**
     * Update model event list and reaction rate caches.
     */
    public void update() {
        if (!ratesDirty)
            return;

        updateModelEventList();

        // Map to store untransformed rates as we iterate over model events

        Map<ModelEvent.RateVariableType,Double> currentRates = new EnumMap<>(ModelEvent.RateVariableType.class);

        // Fill currentRates with rates at start of epidemic process
        // (Note that the epidemic process may start _after_ a model event!)

        currentRates.put(ModelEvent.RateVariableType.INFECTION_RATE,
                getRateAtTime(infectionRateInput.get(),
                              infectionRateShiftTimesInput.get(),
                              infectionRateShiftTimesBackwardInput.get(), 0.0));
        currentRates.put(ModelEvent.RateVariableType.RECOVERY_RATE,
                getRateAtTime(recoveryRateInput.get(),
                              recoveryRateShiftTimesInput.get(),
                              recoveryRateShiftTimesBackwardInput.get(), 0.0));
        currentRates.put(ModelEvent.RateVariableType.PSI_SAMPLING_VARIABLE,
                psiSamplingVariableInput.get() != null
                        ? getRateAtTime(psiSamplingVariableInput.get(),
                                        psiSamplingVariableShiftTimesInput.get(),
                                        psiSamplingVariableShiftTimesBackwardInput.get(), 0.0)
                        : 0.0);
        currentRates.put(ModelEvent.RateVariableType.REMOVAL_PROB,
                getRateAtTime(removalProbInput.get(),
                              removalProbShiftTimesInput.get(),
                              removalProbShiftTimesBackwardInput.get(), 0.0));

        // Initialize transformed rate caches

        if (rateCache.isEmpty()) {
            for (int i=0; i < modelEventList.size()+1; i++)
                rateCache.add(new Double[EpidemicEvent.nTypes]);
        }

        if (removalProbCache == null) {
            removalProbCache = new double[modelEventList.size()+1];
        }

        // Store initial transformed rates

        transformRates(rateCache.get(0), currentRates);
        removalProbCache[0] = currentRates.get(ModelEvent.RateVariableType.REMOVAL_PROB);

        // Iterate over model events, caching transformed rates along the way

        for (int i=0; i<modelEventList.size(); i++) {
            ModelEvent modelEvent = modelEventList.get(i);

            if (modelEvent.time < 0.0)
                continue; // Skip events prior to origin

            if (modelEvent.type == ModelEvent.Type.RATE_CHANGE)
                currentRates.put(modelEvent.rateVariableType, modelEvent.newRateVariableValue);

            transformRates(rateCache.get(i+1), currentRates);
            removalProbCache[i+1] = currentRates.get(ModelEvent.RateVariableType.REMOVAL_PROB);
        }

        ratesDirty = false;
    }


    /**
     * Assemble list of model events.
     */
    public void updateModelEventList() {

        modelEventList.clear();

        // Add rho sampling events
        if (rhoSamplingProbInput.get() != null) {
            for (int i = 0; i < rhoSamplingProbInput.get().getDimension(); i++) {
                if (rhoSamplingProbInput.get().getArrayValue(i) <= 0.0)
                    continue; // Do not add zero-probability rho sampling events

                ModelEvent event = new ModelEvent();
                event.type = ModelEvent.Type.RHO_SAMPLING;
                event.rho = rhoSamplingProbInput.get().getArrayValue(i);
                event.time = getForwardTime(rhoSamplingTimeInput.get(), i,
                        rhoSamplingTimesBackwardInput.get());
                modelEventList.add(event);
            }
        }

        addRateShiftEvents(infectionRateInput.get(),
                infectionRateShiftTimesInput.get(),
                infectionRateShiftTimesBackwardInput.get(),
                ModelEvent.RateVariableType.INFECTION_RATE);

        addRateShiftEvents(recoveryRateInput.get(),
                recoveryRateShiftTimesInput.get(),
                recoveryRateShiftTimesBackwardInput.get(),
                ModelEvent.RateVariableType.RECOVERY_RATE);

        addRateShiftEvents(psiSamplingVariableInput.get(),
                psiSamplingVariableShiftTimesInput.get(),
                psiSamplingVariableShiftTimesBackwardInput.get(),
                ModelEvent.RateVariableType.PSI_SAMPLING_VARIABLE);

        addRateShiftEvents(removalProbInput.get(),
                removalProbShiftTimesInput.get(),
                removalProbShiftTimesBackwardInput.get(),
                ModelEvent.RateVariableType.REMOVAL_PROB);

        Collections.sort(modelEventList);
    }

    public void addRateShiftEvents(Function rateShiftParam,
                                   Function rateShiftTimesParam,
                                   boolean rateShiftTimesBackward,
                                   ModelEvent.RateVariableType rateVariableType) {

        if (rateShiftTimesParam != null) {
            for (int i=0; i<rateShiftTimesParam.getDimension(); i++) {
                ModelEvent event = new ModelEvent();
                event.type = ModelEvent.Type.RATE_CHANGE;
                event.rateVariableType = rateVariableType;
                event.time = getForwardTime(rateShiftTimesParam, i, rateShiftTimesBackward);
                event.newRateVariableValue = getRateInInterval(rateShiftParam, i+1, rateShiftTimesBackward);
                modelEventList.add(event);
            }
        }
    }

    public double getNextModelEventTime(EpidemicState state) {
        update();

        if (state.modelIntervalIdx <modelEventList.size())
            return modelEventList.get(state.modelIntervalIdx).time;
        else
            return Double.POSITIVE_INFINITY;
    }

    public ModelEvent getNextModelEvent(EpidemicState state) {
        update();

        if (state.modelIntervalIdx < modelEventList.size())
            return modelEventList.get(state.modelIntervalIdx);
        else
            return null;
    }

    /**
     * @return list of model events
     */
    public List<ModelEvent> getModelEventList() {
        update();

        return modelEventList;
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

    /**
     * Test whether timeA is less than or equal to timeB within
     * the model's tolerance.
     *
     * @param timeA first time
     * @param timeB second time
     * @return true if timeA <= timeB
     */
    public boolean timesLEQ(double timeA, double timeB) {
        return timeA < timeB || timesEqual(timeA, timeB);
    }

    /**
     * Main method for debugging.
     *
     * @param args
     */
    public static void main(String[] args) throws Exception {
       RealParameter shiftTimesParam = new RealParameter(new Double[] {5.0});
        RealParameter rateParam = new RealParameter(new Double[] {1.0, 0.0});

        SISModel model = new SISModel();
        model.initByName(
                "origin", new RealParameter("10.0"),
                "S0", new IntegerParameter("200"),
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.1"));

        double T = 10.0;
        int nSteps = 5;
        double dt = T/(nSteps-1);
        System.out.println("t\trate");
        for (int i=0; i<nSteps; i++) {
            double t = i*dt;
            System.out.println(t +  "\t" + model.getRateAtTime(rateParam, shiftTimesParam, true, t));
        }

    }

}
