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

package epiinf.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.TreeEventList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which proposes new trajectories conditional on "
        + "the transmission tree and the tree origin.")
public class TrajectoryOperator extends Operator {
    
    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate", "Rate of infection.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "recoveryRate", "Rate of recovery.", Validate.REQUIRED);
    
    public Input<IntegerParameter> S0Input = new Input<IntegerParameter>(
            "S0", "Initial number of susceptibles.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory.", Validate.REQUIRED);

    
    @Override
    public double proposal() {
        double logHR = 0.0;
        
        EpidemicTrajectory traj = trajInput.get();
        double infectionRate = infectionRateInput.get().getValue();
        double recoveryRate = recoveryRateInput.get().getValue();
        
        List<EpidemicEvent> eventList = Lists.newArrayList();
        List<TreeEventList.TreeEvent> treeEventList = treeEventListInput.get().getEventList();

        EpidemicState thisState = new EpidemicState(S0Input.get().getValue(), 1, 0);
        traj.setInitialState(thisState);
        
        double t = 0.0;
        for (TreeEventList.TreeEvent treeEvent : treeEventList) {
            while (true) {
                double infProp = infectionRate*thisState.S*thisState.I;
                double recProp = recoveryRate*thisState.I;
                double totalProp = infProp + recProp;
                
                t += Randomizer.nextExponential(totalProp);
                
                if (t > treeEvent.time) {
                    t = treeEvent.time;
                    break;
                }
                
                EpidemicEvent newEvent = new EpidemicEvent();
                newEvent.time = t;
                
                if (Randomizer.nextDouble()*totalProp < infProp) {
                    newEvent.type = EpidemicEvent.EventType.INFECTION;
                    thisState.S -= 1;
                    thisState.I += 1;
                } else {
                    newEvent.type = EpidemicEvent.EventType.RECOVERY;
                    thisState.I -= 1;
                    thisState.R += 1;
                }
                
                eventList.add(newEvent);
            }
            
            EpidemicEvent newEvent = new EpidemicEvent();
            newEvent.time = treeEvent.time;
            if (treeEvent.type == TreeEventList.TreeEventType.COALESCENCE) {
                newEvent.type = EpidemicEvent.EventType.INFECTION;
                thisState.S -= 1;
                thisState.I += 1;
                
            } else {
                newEvent.type = EpidemicEvent.EventType.RECOVERY;
                thisState.I -= 1;
                thisState.R += 1;
            }
        }
        
        return logHR;
    }
    
    /**
     * Obtain probability density of proposing existing trajectory given
     * current tree.
     * 
     * @return log probability density
     */
    private double getReverseMoveDensity() {
        
        List<TreeEventList.TreeEvent> treeEventList = treeEventListInput.get().getEventList();
        List<EpidemicEvent> eventList = trajInput.get().getEventList();
        List<EpidemicState> stateList = trajInput.get().getStateList();
        
        double infRate = infectionRateInput.get().getValue();
        double recRate = recoveryRateInput.get().getValue();
        
        double logP = 0.0;
        
        int idx=0;
        double lastTime = 0.0;
        for (TreeEventList.TreeEvent treeEvent : treeEventList) {
            
            while (idx<eventList.size() &&
                    !TreeEventList.eventsMatch(treeEvent, eventList.get(idx))) {

                EpidemicState thisState = stateList.get(idx);
                EpidemicEvent thisEvent = eventList.get(idx);
                
                double infProp = infRate*thisState.S*thisState.I;
                double recProp = recRate*thisState.I;
                double totalProp = infProp + recProp;
                
                logP += -totalProp*(thisEvent.time - lastTime);
                
                if (thisEvent.type == EpidemicEvent.EventType.INFECTION)
                    logP += Math.log(infProp);
                else
                    logP += Math.log(recProp);

                lastTime = thisEvent.time;
                idx += 1;
            }
            
            if (idx>=eventList.size())
                return Double.NEGATIVE_INFINITY;
            
            EpidemicState thisState = stateList.get(idx);
            double infProp = infRate*thisState.S*thisState.I;
            double recProp = recRate*thisState.I;
            double totalProp = infProp + recProp;
            
            logP += -totalProp*(treeEvent.time - lastTime);
            
            lastTime = treeEvent.time;
            idx += 1;
        }
        
        return logP;
    }
    
}
