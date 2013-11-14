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
public class BirthDeathTrajectorySimulatorTest {
    
    public BirthDeathTrajectorySimulatorTest() {
    }
    
    @Test
    public void test() throws Exception {
        
        Randomizer.setSeed(42);
        
        BirthDeathTrajectorySimulator trajSim = new BirthDeathTrajectorySimulator();
        trajSim.initByName(
                "birthRate", 1.0,
                "deathRate", 0.1,
                "duration", 5.0);
        
        trajSim.initStateNodes();
        
        assertEquals(trajSim.getEventList().size(), 34);
        assertEquals((long)trajSim.getStateList().get(0).I[0], 1);
        assertEquals((long)trajSim.getStateList().get(34).I[0], 27);

    }
}
