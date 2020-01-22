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

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import epiinf.models.EpidemicModel;
import epiinf.models.SISModel;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simulate an epidemic trajectory under a stochastic SIR model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate an epidemic trajectory under a stochastic SIR model.")
public class SimulatedTrajectory extends EpidemicTrajectory {
    
    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Optional name of file to write simulated trajectory to.");

    public Input<Integer> nStepsInput = new Input<>(
            "nSteps",
            "If nSteps > 0, tau leaping with this many steps will be used to" +
                    " approximate the stochastic integral.", 0);

    public Input<Integer> minSampleCountInput = new Input<>(
            "minSampleCount",
            "Minimum number of samples to accept.",
            0);

    public Input<RealParameter> conditionedSamplingTimesInput = new Input<>(
            "conditionedSamplingTimes",
            "Times at which to force psi-sampling events");

    EpidemicModel model;
    int nSteps, minSampleCount;
    double[] conditionedSamplingTimes;

    public SimulatedTrajectory() { }

    public SimulatedTrajectory(EpidemicModel model, double origin, int nSteps, int minSampleCount,
                               double[] conditionedSamplingTimes) {
        this.model = model;
        this.origin = origin;
        this.nSteps = nSteps;
        this.minSampleCount = minSampleCount;

        this.conditionedSamplingTimes = conditionedSamplingTimes;

        simulationLoop();
    }

    public SimulatedTrajectory(EpidemicModel model, double origin, int nSteps, int minSampleCount) {
        this(model, origin, nSteps, minSampleCount, null);
    }
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        
        model = modelInput.get();
        origin = modelInput.get().getOrigin();
        nSteps = nStepsInput.get();
        minSampleCount = minSampleCountInput.get();

        conditionedSamplingTimes = conditionedSamplingTimesInput.get().getDoubleValues();

        simulationLoop();

        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                    dumpTrajectory(ps);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(SimulatedTrajectory.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    private void simulationLoop() {

        boolean success;
        do {
            eventList = new ArrayList<>();
            stateList = new ArrayList<>();

            success = nSteps > 0
                    ? simulateTL()
                    : simulate();
        } while (!success || eventList.stream().filter(EpidemicEvent::isSample).count()<minSampleCount);
    }


    /**
     * Simulate an epidemic using the provided model.
     *
     */
    public boolean simulate() {

        List<Double> remainingConditionedSamplingTimes = new ArrayList<>();
        if (conditionedSamplingTimes != null) {
            for (double sampTime : conditionedSamplingTimes)
                remainingConditionedSamplingTimes.add(sampTime);
        }

        double endTime;
        endTime = model.getOrigin();

        eventList.clear();
        stateList.clear();

        EpidemicState thisState = model.getInitialState();
        stateList.add(model.getInitialState());

        thisState.time = 0;

        while (true) {
            model.calculatePropensities(thisState);

            double totalPropensity = model.propensities[EpidemicEvent.INFECTION]
                    + model.propensities[EpidemicEvent.RECOVERY]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE];

            double dt;
            if (totalPropensity>0.0)
                dt = Randomizer.nextExponential(totalPropensity);
            else
                dt = Double.POSITIVE_INFINITY;

            thisState.time += dt;

            EpidemicEvent nextEvent = new EpidemicEvent();

            double nextModelEventTime = model.getNextModelEventTime(thisState);
            double nextConditionedSamplingTime = remainingConditionedSamplingTimes.isEmpty()
                    ? Double.POSITIVE_INFINITY
                    : remainingConditionedSamplingTimes.get(0);
            double nextModelEventOrSamplingTime = Math.min(nextModelEventTime, nextConditionedSamplingTime);

            if (nextModelEventOrSamplingTime <= endTime && thisState.time > nextModelEventOrSamplingTime) {
                if (nextModelEventOrSamplingTime == nextConditionedSamplingTime) {
                    if (thisState.I <= 0.0)
                        return false;

                    if (model.currentRemovalProb == 1.0 || Randomizer.nextDouble() < model.currentRemovalProb)
                        nextEvent.type = EpidemicEvent.PSI_SAMPLE_REMOVE;
                    else
                        nextEvent.type = EpidemicEvent.PSI_SAMPLE_NOREMOVE;

                    nextEvent.time = nextConditionedSamplingTime;
                    model.incrementState(thisState, nextEvent);
                    eventList.add(nextEvent);
                    stateList.add(thisState.copy());

                    remainingConditionedSamplingTimes.remove(0);

                } else {
                    ModelEvent event = model.getNextModelEvent(thisState);

                    if (event.type == ModelEvent.Type.RHO_SAMPLING) {
                        nextEvent.type = EpidemicEvent.RHO_SAMPLE;

                        // Got to be a better way of sampling from a binomial distribution
                        nextEvent.multiplicity = 0;
                        for (int i = 0; i < thisState.I; i++) {
                            if (Randomizer.nextDouble() < event.rho)
                                nextEvent.multiplicity += 1;
                        }


                        nextEvent.time = event.time;
                        thisState.time = event.time;

                        model.incrementState(thisState, nextEvent);
                        eventList.add(nextEvent);
                        stateList.add(thisState.copy());
                    }

                    thisState.modelIntervalIdx += 1;
                }

                continue;
            }

            if (thisState.time>=endTime)
                break;


            nextEvent.time = thisState.time;

            double u = totalPropensity*Randomizer.nextDouble();

            for (int type = 0; type<EpidemicEvent.nTypes; type++) {
                u -= model.propensities[type];

                if (u<0) {
                    nextEvent.type = type;
                    break;
                }
            }

            model.incrementState(thisState, nextEvent);
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }

        return true;
    }

    /**
     * Simulate an epidemic using the provided model with tau leaping.
     *
     */
    public boolean simulateTL() {

        List<Double> remainingConditionedSamplingTimes = new ArrayList<>();
        if (conditionedSamplingTimes != null) {
            for (double sampTime : conditionedSamplingTimes)
                remainingConditionedSamplingTimes.add(sampTime);
        }

        double endTime;
        endTime = model.getOrigin();

        double dt = endTime/(nSteps-1);

        eventList.clear();
        stateList.clear();

        EpidemicState thisState = model.getInitialState();
        stateList.add(model.getInitialState());

        thisState.time = 0;

        for (int tidx = 1; tidx<nSteps; tidx++) {
            model.calculatePropensities(thisState);

            double nextModelEventTime = model.getNextModelEventTime(thisState);
            double nextConditionedSamplingTime = remainingConditionedSamplingTimes.isEmpty()
                    ? Double.POSITIVE_INFINITY
                    : remainingConditionedSamplingTimes.get(0);
            double nextModelEventOrSamplingTime = Math.min(nextModelEventTime, nextConditionedSamplingTime);

            double trueDt = Math.min(dt, nextModelEventOrSamplingTime - thisState.time);

            EpidemicEvent infectEvent = new EpidemicEvent();
            infectEvent.type = EpidemicEvent.INFECTION;
            infectEvent.multiplicity = Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.INFECTION]));
            model.incrementState(thisState, infectEvent);
            infectEvent.time = thisState.time + trueDt;
            eventList.add(infectEvent);

            EpidemicEvent recovEvent = new EpidemicEvent();
            recovEvent.type = EpidemicEvent.RECOVERY;
            recovEvent.multiplicity = Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.RECOVERY]));
            model.incrementState(thisState, recovEvent);
            recovEvent.time = thisState.time + trueDt;
            eventList.add(recovEvent);

            EpidemicEvent psiSampRemoveEvent = new EpidemicEvent();
            psiSampRemoveEvent.type = EpidemicEvent.PSI_SAMPLE_REMOVE;
            psiSampRemoveEvent.multiplicity = Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]));
            model.incrementState(thisState, psiSampRemoveEvent);
            psiSampRemoveEvent.time = thisState.time + trueDt;
            eventList.add(psiSampRemoveEvent);

            EpidemicEvent psiSampNoRemoveEvent = new EpidemicEvent();
            psiSampNoRemoveEvent.type = EpidemicEvent.PSI_SAMPLE_NOREMOVE;
            psiSampNoRemoveEvent.multiplicity = Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]));
            psiSampNoRemoveEvent.time = thisState.time + trueDt;
            eventList.add(psiSampNoRemoveEvent);


            if (trueDt < dt) {

                if (nextModelEventOrSamplingTime == nextConditionedSamplingTime) {
                    if (thisState.I <= 0.0)
                        return false;

                    EpidemicEvent samplingEvent = new EpidemicEvent();

                    if (model.currentRemovalProb == 1.0 || Randomizer.nextDouble() < model.currentRemovalProb)
                        samplingEvent.type = EpidemicEvent.PSI_SAMPLE_REMOVE;
                    else
                        samplingEvent.type = EpidemicEvent.PSI_SAMPLE_NOREMOVE;

                    samplingEvent.time = nextConditionedSamplingTime;
                    model.incrementState(thisState, samplingEvent);
                    eventList.add(samplingEvent);
                    stateList.add(thisState.copy());

                    remainingConditionedSamplingTimes.remove(0);

                } else {
                    ModelEvent event = model.getNextModelEvent(thisState);

                    if (event.type == ModelEvent.Type.RHO_SAMPLING) {

                        EpidemicEvent rhoSampEvent = new EpidemicEvent();
                        rhoSampEvent.type = EpidemicEvent.RHO_SAMPLE;
                        rhoSampEvent.time = event.time;

                        // Got to be a better way of sampling from a binomial distribution
                        rhoSampEvent.multiplicity = 0;
                        for (int i = 0; i < thisState.I; i++) {
                            if (Randomizer.nextDouble() < event.rho)
                                rhoSampEvent.multiplicity += 1;
                        }

                        model.incrementState(thisState, rhoSampEvent);
                    }

                    thisState.modelIntervalIdx += 1;
                }
            }

            // Rough error correction - doesn't conserve number.
            if (!thisState.isValid()) {
                thisState.I = Math.max(0, thisState.I);
                thisState.S = Math.max(0, thisState.S);
                thisState.R = Math.max(0, thisState.R);
            }

            thisState.time += trueDt;

            stateList.add(thisState.copy());
        }

        return true;
    }

    /**
     * Main method for debugging.
     *
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        SISModel model = new SISModel();
        model.initByName(
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.1"),
                "S0", new IntegerParameter("199"));
        SimulatedTrajectory traj = new SimulatedTrajectory(model, 10, 1000, 0);

        System.out.println(traj);


    }
}
