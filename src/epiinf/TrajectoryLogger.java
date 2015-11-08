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

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import epiinf.distribs.SMCTreeDensity;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TrajectoryLogger extends BEASTObject implements Loggable {

    public Input<SMCTreeDensity> treeDensityInput = new Input<>("treeDensity",
            "SMC Tree density from which to log trajectories.",
            Input.Validate.REQUIRED);

    SMCTreeDensity treeDensity;

    public TrajectoryLogger() { }

    @Override
    public void initAndValidate() throws Exception {
        treeDensity = treeDensityInput.get();
    }

    @Override
    public void init(PrintStream out) throws Exception {
        out.print("trajectory" + "\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {

        treeDensity.calculateLogP(true);
        if (treeDensity.getRecordedTrajectory().isEmpty()) {
            out.print("NA\t");
            return;
        }

        boolean isFirst = true;
        for (EpidemicState state : treeDensity.getRecordedTrajectory()) {
            if (!isFirst)
                out.print(",");
            else
                isFirst = false;

            out.print((treeDensity.getRecordedOrigin() - state.time)
                    + ":" + state.S + ":" + state.I + ":" + state.R);
        }

        out.print("\t");
    }

    @Override
    public void close(PrintStream out) { }
}
