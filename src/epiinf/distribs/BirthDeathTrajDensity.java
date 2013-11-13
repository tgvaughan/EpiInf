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

package epiinf.distribs;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.SIRTrajectorySimulator;
import java.io.PrintStream;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BirthDeathTrajDensity extends Distribution {

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory",
            "Epidemic trajectory state object.", Validate.REQUIRED);
    
    public Input<RealParameter> birthRateInput = new Input<RealParameter>(
            "birthRate",
            "Birth rate (per infected).", Validate.REQUIRED);

    public Input<RealParameter> deathRateInput = new Input<RealParameter>(
            "deathRate",
            "Recovery rate (per susceptible per infected).", Validate.REQUIRED);

    
    private EpidemicTrajectory trajectory;
    private RealParameter birthRate, deathRate;
    
    public BirthDeathTrajDensity() { }
    
    @Override
    public void initAndValidate() {
        trajectory = trajectoryInput.get();
        birthRate = birthRateInput.get();
        deathRate = deathRateInput.get();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
        List<EpidemicState> stateList = trajectory.getStateList();
        
        double lastEventTime = 0.0;
                
        for (int i=0; i<trajectory.getEventList().size(); i++) {
            EpidemicEvent thisEvent = trajectory.getEventList().get(i);
            EpidemicState thisState = trajectory.getStateList().get(i);
            
            double birthProp = thisState.I*birthRate.getValue();
            double deathProp = thisState.I*deathRate.getValue();
            double totalProp = birthProp + deathProp;
            
            logP += -totalProp*(thisEvent.time-lastEventTime);
            lastEventTime = thisEvent.time;
            
            if (thisEvent.type == EpidemicEvent.EventType.INFECTION)
                logP += Math.log(birthProp);
            else
                logP += Math.log(deathProp);
        }
        
        return logP;
    }

    @Override
    public boolean requiresRecalculation() {
        return super.requiresRecalculation(); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) { }
    
}
