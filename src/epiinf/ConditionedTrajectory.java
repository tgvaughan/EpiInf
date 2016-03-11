/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
import epiinf.distribs.TrajectoryRecorder;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ConditionedTrajectory extends EpidemicTrajectory {

    public Input<TrajectoryRecorder> treeDensityInput = new Input<>("treeDensity",
            "SMC Tree density from which to log trajectories.",
            Input.Validate.REQUIRED);

    TrajectoryRecorder treeDensity;

    public ConditionedTrajectory() {
    }

    @Override
    public void initAndValidate() {
        treeDensity = treeDensityInput.get();
        super.initAndValidate();
    }

    @Override
    public void log(int nSample, PrintStream out) {
        assignFrom(treeDensity.getConditionedTrajectory());
        super.log(nSample, out);
    }
}
