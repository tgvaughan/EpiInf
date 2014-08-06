/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.core.Input;
import beast.core.Input.Validate;
import epiinf.models.EpidemicModel;

/**
 * Simulate an epidemic trajectory under the same model used for inference.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TrajectorySimulator extends beast.core.Runnable {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<Double> durationInput = new Input<>(
            "maxDuration", "Maximum duration of epidemic to simulate. "
                    + "Defaults to infinity.", Double.POSITIVE_INFINITY);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to write simulated trajectory to.", Validate.REQUIRED);
    
    @Override
    public void initAndValidate() { }
    
    @Override
    public void run() throws Exception {
        
        SimulatedTrajectory traj = new SimulatedTrajectory();
        traj.initByName(
                "model", modelInput.get(),
                "maxDuration", durationInput.get(),
                "fileName", fileNameInput.get());
    }
}
