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
public class BirthDeathTrajectoryOperator extends Operator {
    
    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<RealParameter> birthRateInput = new Input<RealParameter>(
            "birthRate", "Birth rate (per infected).", Validate.REQUIRED);
    
    public Input<RealParameter> deathRateInput = new Input<RealParameter>(
            "deathRate", "Death rate (per infected).", Validate.REQUIRED);
    
    public Input<EpidemicTrajectory> trajInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory.", Validate.REQUIRED);

    @Override
    public void initAndValidate() { }
    
    @Override
    public double proposal() {
        
        // Initialize HR with reverse move density
        double logHR = getReverseMoveDensity();
        
        EpidemicTrajectory traj = trajInput.get();
        double birthRate = birthRateInput.get().getValue();
        double deathRate = deathRateInput.get().getValue();
        
        List<EpidemicEvent> eventList = Lists.newArrayList();
        List<EpidemicState> stateList = Lists.newArrayList();
        List<TreeEventList.TreeEvent> treeEventList = treeEventListInput.get().getEventList();

        EpidemicState thisState = new EpidemicState(0, 1, 0);
        stateList.add(thisState.copy());
        
        double t = 0.0;
        for (TreeEventList.TreeEvent treeEvent : treeEventList) {
            
            while (true) {
                double birthProp = birthRate*thisState.I;
                double deathProp = deathRate*thisState.I;
                double totalProp = birthProp + deathProp;
                
                double dt = Randomizer.nextExponential(totalProp);
                logHR -= -totalProp*(Math.min(t+dt, treeEvent.time)-t);
                
                t += dt;
                
                if (t > treeEvent.time) {
                    t = treeEvent.time;
                    break;
                }
                
                EpidemicEvent newEvent = new EpidemicEvent();
                newEvent.time = t;
                
                if (Randomizer.nextDouble()*totalProp < birthProp) {
                    newEvent.type = EpidemicEvent.EventType.INFECTION;
                    thisState.I += 1;
                    logHR -= Math.log(birthProp);
                } else {
                    newEvent.type = EpidemicEvent.EventType.RECOVERY;
                    thisState.I -= 1;
                    logHR -= Math.log(deathProp);
                }
                
                if (!thisState.isValid())
                    return Double.NEGATIVE_INFINITY;
                
                eventList.add(newEvent);
                stateList.add(thisState.copy());
            }
            
            EpidemicEvent newEvent = new EpidemicEvent();
            newEvent.time = treeEvent.time;
            if (treeEvent.type == TreeEventList.TreeEventType.COALESCENCE) {
                newEvent.type = EpidemicEvent.EventType.INFECTION;
                thisState.I += 1;
                
            } else {
                newEvent.type = EpidemicEvent.EventType.RECOVERY;
                thisState.I -= 1;
            }
            
            if (!thisState.isValid())
                return Double.NEGATIVE_INFINITY;
            
            eventList.add(newEvent);
            stateList.add(thisState.copy());
        }
        
        traj.setEventAndStateList(eventList, stateList);
        
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
        
        double birthRate = birthRateInput.get().getValue();
        double deathRate = deathRateInput.get().getValue();
        
        double logP = 0.0;
        
        int idx=0;
        double lastTime = 0.0;
        for (TreeEventList.TreeEvent treeEvent : treeEventList) {
            
            while (idx<eventList.size() &&
                    !TreeEventList.eventsMatch(treeEvent, eventList.get(idx))) {

                EpidemicState thisState = stateList.get(idx);
                EpidemicEvent thisEvent = eventList.get(idx);
                
                double birthProp = birthRate*thisState.I;
                double deathProp = deathRate*thisState.I;
                double totalProp = birthProp + deathProp;
                
                logP += -totalProp*(thisEvent.time - lastTime);
                
                if (thisEvent.type == EpidemicEvent.EventType.INFECTION)
                    logP += Math.log(birthProp);
                else
                    logP += Math.log(deathProp);

                lastTime = thisEvent.time;
                idx += 1;
            }
            
            if (idx>=eventList.size())
                return Double.NEGATIVE_INFINITY;
            
            EpidemicState thisState = stateList.get(idx);
            double birthProp = birthRate*thisState.I;
            double deathProp = deathRate*thisState.I;
            double totalProp = birthProp + deathProp;
            
            logP += -totalProp*(treeEvent.time - lastTime);
            
            lastTime = treeEvent.time;
            idx += 1;
        }
        
        return logP;
    }
    
}
