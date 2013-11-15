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
import epiinf.EpidemicModel;
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
    
    public Input<EpidemicModel> modelInput = new Input<EpidemicModel>(
            "model",
            "Epidemic model.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory",
            "Epidemic trajectory state object.", Validate.REQUIRED);
    
    private EpidemicModel model;
    private EpidemicTrajectory trajectory;
    
    public TrajDensity() { }
    
    @Override
    public void initAndValidate() {
        trajectory = trajectoryInput.get();
        model = modelInput.get();
    }
    
    @Override
    public double calculateLogP() {
        
        return model.getPathProbability(0, trajectory.getDuration(),
                model.getInitialState(), trajectory.getEventList());
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
