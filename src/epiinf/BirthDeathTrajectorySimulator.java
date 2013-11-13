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
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.util.Randomizer;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simulate an epidemic trajectory under a stochastic SIR model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate an epidemic trajectory under a stochastic SIR model.")
public class BirthDeathTrajectorySimulator extends EpidemicTrajectory implements StateNodeInitialiser {
    
    public Input<Double> birthRateInput = new Input<Double>(
            "birthRate", "Rate of birth (per infected).",
            Validate.REQUIRED);

    public Input<Double> deathRateInput = new Input<Double>(
            "deathRate", "Rate of death (per infected)",
            Validate.REQUIRED);
    
    public Input<Double> durationInput = new Input<Double>(
            "duration", "Duration of epidemic.",
            Validate.REQUIRED);
    
    public Input<String> fileNameInput = new Input<String>(
            "fileName",
            "Optional name of file to write simulated trajectory to.");
    
    double birthRate, deathRate;
    double duration;
    
    public BirthDeathTrajectorySimulator() { }
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        birthRate = birthRateInput.get();
        deathRate = deathRateInput.get();
        duration = durationInput.get();
        
        simulate();
        
        if (fileNameInput.get() != null) {
            try {
                PrintStream ps = new PrintStream(fileNameInput.get());
                dumpTrajectory(ps);
                ps.close();
            } catch (FileNotFoundException ex) {
                Logger.getLogger(BirthDeathTrajectorySimulator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    private void simulate() {
        eventList.clear();
        stateList.clear();
        
        EpidemicState thisState = new EpidemicState(0, 1, 0);
        stateList.add(thisState.copy());
        
        double t = 0.0;

        while (true) {
            
            double infProp = thisState.I*birthRate;
            double recProp = thisState.I*deathRate;
            double totalProp = infProp + recProp;

            t += Randomizer.nextExponential(totalProp);
            if (t>duration) {
                t = duration;
                break;
            }
            
            EpidemicEvent nextEvent = new EpidemicEvent();
            nextEvent.time = t;

            if (Randomizer.nextDouble()*totalProp < infProp) {
                nextEvent.type = EpidemicEvent.EventType.INFECTION;
                thisState.I += 1;
            } else {
                nextEvent.type = EpidemicEvent.EventType.RECOVERY;
                thisState.I -= 1;
            }
            
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }
    }
    
    @Override
    public void initStateNodes() throws Exception {
        //simulate();
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
