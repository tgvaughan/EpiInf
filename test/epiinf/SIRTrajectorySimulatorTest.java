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
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SIRTrajectorySimulatorTest {
    
    public SIRTrajectorySimulatorTest() {
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
        
        // Assert some arbitary features of trajectory:
        assertEquals(trajSim.getEventList().size(), 1983);
        assertEquals((long)trajSim.getStateList().get(0).S[0], 1000);
        assertEquals((long)trajSim.getStateList().get(0).I[0], 1);
        assertEquals((long)trajSim.getStateList().get(0).R[0], 0);
        
        EpidemicState finalState = trajSim.getStateList()
                .get(trajSim.getStateList().size()-1);
        
        assertEquals((long)finalState.I[0], 0);
        assertEquals((long)finalState.R[0], 1001-(long)finalState.S[0]);
    }
}
