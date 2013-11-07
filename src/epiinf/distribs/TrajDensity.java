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
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TrajDensity extends Distribution {

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory",
            "Epidemic trajectory state object.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate",
            "Infection rate (per susceptible per infected).", Validate.REQUIRED);

    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "infectionRate",
            "Infection rate (per susceptible per infected).", Validate.REQUIRED);

    
    private EpidemicTrajectory trajectory;
    private RealParameter infectionRate, recoveryRate;
    
    public TrajDensity() { }
    
    @Override
    public void initAndValidate() {
        trajectory = trajectoryInput.get();
        infectionRate = infectionRateInput.get();
        recoveryRate = recoveryRateInput.get();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
        EpidemicState thisState = trajectory.getInitialState();
        double t = 0.0;

        for (EpidemicEvent event : trajectory.getEventList()) {
            
            double infectionProp = thisState.S*thisState.I*infectionRate.getValue();
            double recoveryProp = thisState.I*recoveryRate.getValue();
            double totalProp = infectionProp + recoveryProp;
            
            logP += -totalProp*(event.time-t);
            t = event.time;
            
            if (event.type == EpidemicEvent.EventType.INFECTION)
                logP += Math.log(infectionProp);
            else
                logP += Math.log(recoveryProp);
            
        }
        
        return logP;
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
