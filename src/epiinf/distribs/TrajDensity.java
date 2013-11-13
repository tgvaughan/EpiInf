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
public class TrajDensity extends Distribution {

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory",
            "Epidemic trajectory state object.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate",
            "Infection rate (per susceptible per infected).", Validate.REQUIRED);

    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "recoveryRate",
            "Recovery rate (per susceptible per infected).", Validate.REQUIRED);

    
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
        
        List<EpidemicState> stateList = trajectory.getStateList();
        
        double lastEventTime = 0.0;
                
        for (int i=0; i<trajectory.getEventList().size(); i++) {
            EpidemicEvent thisEvent = trajectory.getEventList().get(i);
            EpidemicState thisState = trajectory.getStateList().get(i);
            
            double infectionProp = thisState.S*thisState.I*infectionRate.getValue();
            double recoveryProp = thisState.I*recoveryRate.getValue();
            double totalProp = infectionProp + recoveryProp;
            
            logP += -totalProp*(thisEvent.time-lastEventTime);
            lastEventTime = thisEvent.time;
            
            if (thisEvent.type == EpidemicEvent.EventType.INFECTION)
                logP += Math.log(infectionProp);
            else
                logP += Math.log(recoveryProp);
            
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

    /**
     * Main method for debugging only.  Generates cross-sections of the
     * parameter likelihood surface and writes these to files for further
     * analysis.
     * 
     * @param args 
     * @throws java.lang.Exception 
     */
    public static void main (String [] args) throws Exception {
               Randomizer.setSeed(42);
        
        SIRTrajectorySimulator trajSim = new SIRTrajectorySimulator();
        trajSim.initByName(
                "S0", 1000,
                "I0", 1,
                "R0", 0,
                "infectionRate", 0.001,
                "recoveryRate", 0.2);
        
        trajSim.initStateNodes();
        
        PrintStream ps = new PrintStream("infectionRateLikelihood.txt");
        ps.println("infectionRate logP");
        
        TrajDensity trajDensity = new TrajDensity();
        for (int i=0; i<100; i++) {
            double infectionRate = Math.pow(10, -4 + (i/99.0)*(-2 - -4));
            trajDensity.initByName(
                "epidemicTrajectory", trajSim,
                "infectionRate", new RealParameter(String.valueOf(infectionRate)),
                "recoveryRate", new RealParameter("0.2"));
            
            ps.println(infectionRate + " " + trajDensity.calculateLogP());
        }
        
        ps.close();
        
        ps = new PrintStream("recoveryRateLikelihood.txt");
        ps.println("recoveryRate logP");
        
        for (int i=0; i<100; i++) {
            double recoveryRate = Math.pow(10, Math.log10(0.2)-1 + (i/99.0)*2.0);
            trajDensity.initByName(
                "epidemicTrajectory", trajSim,
                "infectionRate", new RealParameter("0.001"),
                "recoveryRate", new RealParameter(String.valueOf(recoveryRate)));
            
            ps.println(recoveryRate + " " + trajDensity.calculateLogP());
        }
        
        ps.close();
    }
    
}
