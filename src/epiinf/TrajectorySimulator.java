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
import beast.core.parameter.RealParameter;
import epiinf.models.EpidemicModel;

/**
 * Simulate an epidemic trajectory under the same model used for inference.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TrajectorySimulator extends beast.core.Runnable {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to write simulated trajectory to.", Validate.REQUIRED);

    public Input<Integer> nStepsInput = new Input<>(
            "nSteps",
            "If nSteps > 0, tau leaping with this many steps will be used to" +
                    " approximate the stochastic integral.", 0);

    public Input<Integer> minSampleCountInput = new Input<>(
            "minSampleCount",
            "Minimum number of samples to accept.",
            0);

    @Override
    public void initAndValidate() { }
    
    @Override
    public void run() throws Exception {
        
        SimulatedTrajectory traj = new SimulatedTrajectory();
        traj.initByName(
                "model", modelInput.get(),
                "fileName", fileNameInput.get(),
                "nSteps", nStepsInput.get(),
                "minSampleCount", minSampleCountInput.get());
    }
}
