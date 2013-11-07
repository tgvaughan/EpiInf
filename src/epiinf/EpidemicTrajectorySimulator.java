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

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.util.Randomizer;
import java.util.List;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpidemicTrajectorySimulator extends EpidemicTrajectory implements StateNodeInitialiser {
    
    public Input<Integer> S0Input = new Input<Integer>(
            "S0", "Initial susceptible count.",
            Validate.REQUIRED);

    public Input<Integer> I0Input = new Input<Integer>(
            "I0", "Initial susceptible count.",
            Validate.REQUIRED);
    
    public Input<Integer> R0Input = new Input<Integer>(
            "R0", "Initial susceptible count.",
            Validate.REQUIRED);
    
    public Input<Double> infectionRateInput = new Input<Double>(
            "infectionRate", "Rate of infection (per susceptible per infected)",
            Validate.REQUIRED);

    public Input<Double> recoveryRateInput = new Input<Double>(
            "recoveryRate", "Rate of recovery (per infected)",
            Validate.REQUIRED);
    
    double infectionRate, recoveryRate;
    
    public EpidemicTrajectorySimulator() { }
    
    @Override
    public void initAndValidate() {
        infectionRate = infectionRateInput.get();
        recoveryRate = recoveryRateInput.get();
    }
    
    private void simulate() {
        eventList.clear();
        
        EpidemicState thisState = initialState.copy();

        while (thisState.I>0) {
            
            double infProp = thisState.S*thisState.I*infectionRate;
            double recProp = thisState.I*recoveryRate;
            double totalProp = infProp*recProp;

            Event nextEvent = new Event();
            nextEvent.waitingTime = Randomizer.nextExponential(totalProp);

            if (Randomizer.nextDouble()<infProp) {
                nextEvent.type = EventType.INFECTION;
                thisState.S -= 1;
                thisState.I += 1;
            } else {
                nextEvent.type = EventType.RECOVERY;
                thisState.I -= 1;
                thisState.R += 1;
            }
            
            eventList.add(nextEvent);
        }
    }
    
    @Override
    public void initStateNodes() throws Exception {
        initialState = new EpidemicState(
                S0Input.get(), I0Input.get(), R0Input.get());
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
