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

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import epiinf.models.SIRModel;
import epiinf.TrajectorySimulator;
import epiinf.TransmissionTreeSimulator;
import epiinf.TreeEventList;
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
        Randomizer.setSeed(2785);
        
        SIRModel model = new SIRModel();
        model.initByName(
                "S0", new IntegerParameter("999"),
                "infectionRate", new RealParameter("0.001"),
                "recoveryRate", new RealParameter("0.2"));
        
        TrajectorySimulator trajSim = new TrajectorySimulator();
        trajSim.initByName("model", model);
        
        Tree tree = new Tree();
        RealParameter treeOrigin = new RealParameter();
        
        TransmissionTreeSimulator treeSim = new TransmissionTreeSimulator();
        treeSim.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin,
                "epidemicTrajectory", trajSim,
                "nLeaves", 100,
                "model", model);
        treeSim.initStateNodes();
        
        TreeEventList treeEventList = new TreeEventList();
        treeEventList.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin);
        
        TreeProb treeDensity = new TreeProb();
        treeDensity.initByName(
                "treeEventList", treeEventList,
                "epidemicTrajectory", trajSim,
                "model", model);

        double logP = treeDensity.calculateLogP();
        double logPtruth = -281.686741693071;
        
        assertTrue(Math.abs(logP-logPtruth)<1e-10);
    }
}
