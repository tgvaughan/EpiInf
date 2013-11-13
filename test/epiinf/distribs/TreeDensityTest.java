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

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import epiinf.SIRTrajectorySimulator;
import epiinf.TransmissionTreeSimulator;
import epiinf.TreeEventList;
import java.io.PrintStream;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TreeDensityTest {
    
    public TreeDensityTest() {
    }

    @Test
    public void test() throws Exception {
        Randomizer.setSeed(42);
        
        SIRTrajectorySimulator trajSim = new SIRTrajectorySimulator();
        trajSim.initByName(
                "S0", 1000,
                "infectionRate", 0.001,
                "recoveryRate", 0.2);
        
        trajSim.initStateNodes();
        
        Tree tree = new Tree();
        RealParameter treeOrigin = new RealParameter();
        
        TransmissionTreeSimulator treeSim = new TransmissionTreeSimulator();
        treeSim.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin,
                "epidemicTrajectory", trajSim,
                "nLeaves", 100);
        treeSim.initStateNodes();
        
        TreeEventList treeEventList = new TreeEventList();
        treeEventList.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin);
        
        TreeDensity treeDensity = new TreeDensity();
        treeDensity.initByName(
                "treeEventList", treeEventList,
                "epidemicTrajectory", trajSim);

        double logP = treeDensity.calculateLogP();
        double logPtruth = -284.84151224202884;
        
        assertTrue(2.0*Math.abs(logP-logPtruth)/(logP+logPtruth)<1e-10);
    }
}
