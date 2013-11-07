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

import beast.util.Randomizer;
import java.io.PrintStream;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpidemicTrajectorySimulatorTest {
    
    public EpidemicTrajectorySimulatorTest() {
    }
    
    @Test
    public void test() throws Exception {
        
        Randomizer.setSeed(42);
        
        EpidemicTrajectorySimulator trajSim = new EpidemicTrajectorySimulator();
        trajSim.initByName(
                "S0", 1000,
                "I0", 1,
                "R0", 0,
                "infectionRate", 0.001,
                "recoveryRate", 0.2);
        
        trajSim.initStateNodes();
        
        // Assert some arbitary features of trajectory:
        assertEquals(trajSim.getEventList().size(), 1983);
        assertEquals(trajSim.getStateList().get(0).S, 1000);
        assertEquals(trajSim.getStateList().get(0).I, 1);
        assertEquals(trajSim.getStateList().get(0).R, 0);
        
        EpidemicState finalState = trajSim.getStateList()
                .get(trajSim.getStateList().size()-1);
        
        assertEquals(finalState.I, 0);
        assertEquals(finalState.R, 1001-finalState.S);
    }
}
