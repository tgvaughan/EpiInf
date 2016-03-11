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
import beast.core.Function;
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
    
    public Input<Double> durationInput = new Input<>(
            "maxDuration", "Maximum duration of epidemic to simulate. "
                    + "Defaults to infinity.", Double.POSITIVE_INFINITY);

    public Input<Function> originInput = new Input<>(
            "origin", "Origin with respect to most recent sample in tree. " +
            "If provided, trimes will be logged as ages before most recent " +
            "sample.");

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Optional name of file to write simulated trajectory to.");

    public Input<Integer> nStepsInput = new Input<>(
            "nSteps",
            "If nSteps > 0, tau leaping with this many steps will be used to" +
                    " approximate the stochastic integral.", 0);
    
    EpidemicModel model;
    double duration;
    int nSteps;
    
    public SimulatedTrajectory() { }

    public SimulatedTrajectory(EpidemicModel model, double duration, int nSteps) {
        this.model = model;
        this.duration = duration;
        this.nSteps = nSteps;

        eventList = new ArrayList<>();
        stateList = new ArrayList<>();

        if (nSteps > 0)
            simulateTL();
        else
            simulate();

        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                dumpTrajectory(ps);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(SimulatedTrajectory.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        
        model = modelInput.get();
        duration = durationInput.get();
        nSteps = nStepsInput.get();

        if (nSteps > 0)
            simulateTL();
        else
            simulate();

        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                    dumpTrajectory(ps);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(SimulatedTrajectory.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }


    /**
     * Simulate an epidemic using the provided model.
     *
     */
    public void simulate() {

        double endTime;
        if (origin != null)
            endTime = origin;
        else
            endTime = duration;

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
            if (nextModelEventTime < endTime && thisState.time > nextModelEventTime) {
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
    }

    /**
     * Simulate an epidemic using the provided model with tau leaping.
     *
     */
    public void simulateTL() {

        double endTime;
        if (origin != null)
            endTime = origin;
        else
            endTime = duration;

        double dt = endTime/(nSteps-1);

        eventList.clear();
        stateList.clear();

        EpidemicState thisState = model.getInitialState();
        stateList.add(model.getInitialState());

        thisState.time = 0;

        for (int tidx = 1; tidx<nSteps; tidx++) {
            model.calculatePropensities(thisState);
//
//            double totalPropensity = model.propensities[EpidemicEvent.INFECTION]
//                    + model.propensities[EpidemicEvent.RECOVERY]
//                    + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
//                    + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE];

            double nextModelEventTime = model.getNextModelEventTime(thisState);
            double trueDt = Math.min(dt, nextModelEventTime - thisState.time);

            EpidemicEvent infectEvent = new EpidemicEvent();
            infectEvent.type = EpidemicEvent.INFECTION;
            infectEvent.multiplicity = (int)Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.INFECTION]));
            model.incrementState(thisState, infectEvent);
            infectEvent.time = thisState.time + trueDt;
            eventList.add(infectEvent);

            EpidemicEvent recovEvent = new EpidemicEvent();
            recovEvent.type = EpidemicEvent.RECOVERY;
            recovEvent.multiplicity = (int)Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.RECOVERY]));
            model.incrementState(thisState, recovEvent);
            recovEvent.time = thisState.time + trueDt;
            eventList.add(recovEvent);

            EpidemicEvent psiSampRemoveEvent = new EpidemicEvent();
            psiSampRemoveEvent.type = EpidemicEvent.PSI_SAMPLE_REMOVE;
            psiSampRemoveEvent.multiplicity = (int)Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]));
            model.incrementState(thisState, psiSampRemoveEvent);
            psiSampRemoveEvent.time = thisState.time + trueDt;
            eventList.add(psiSampRemoveEvent);

            EpidemicEvent psiSampNoRemoveEvent = new EpidemicEvent();
            psiSampNoRemoveEvent.type = EpidemicEvent.PSI_SAMPLE_NOREMOVE;
            psiSampNoRemoveEvent.multiplicity = (int)Math.round(
                    Randomizer.nextPoisson(trueDt*model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]));
            psiSampNoRemoveEvent.time = thisState.time + trueDt;
            eventList.add(psiSampNoRemoveEvent);


            if (trueDt < dt) {
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

            // Rough error correction - doesn't conserve number.
            if (!thisState.isValid()) {
                thisState.I = Math.max(0, thisState.I);
                thisState.S = Math.max(0, thisState.S);
                thisState.R = Math.max(0, thisState.R);
            }

            thisState.time += trueDt;

            stateList.add(thisState.copy());
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {
        if (originInput.get() != null)
            origin = originInput.get().getArrayValue();

        simulate();

        super.log(nSample, out);
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
        SimulatedTrajectory traj = new SimulatedTrajectory(model, 10, 1000);

        System.out.println(traj);


    }
}
