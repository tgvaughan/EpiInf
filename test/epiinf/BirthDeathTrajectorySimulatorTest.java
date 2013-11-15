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

import beast.core.parameter.RealParameter;
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
        
        Randomizer.setSeed(1);
        
        BirthDeathModel model = new BirthDeathModel();
        model.initByName(
                "birthRate", new RealParameter("2.0"),
                "deathRate", new RealParameter("1.0"));
        
        TrajectorySimulator trajSim = new TrajectorySimulator();
        trajSim.initByName(
                "model", model,
                "maxDuration", 5.0);
        
        trajSim.initStateNodes();
        
        assertEquals(trajSim.getEventList().size(), 151);
        assertEquals((long)trajSim.getStateList().get(0).I, 1);
        assertEquals((long)trajSim.getStateList().get(151).I, 48);

    }
}
