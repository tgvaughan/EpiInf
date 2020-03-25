/*
 * Copyright (C) 2017 Tim Vaughan <tgvaughan@gmail.com>
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
import epiinf.distribs.SMCTreeDensity;
import epiinf.models.EpidemicModel;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TrajectoryLogger extends BEASTObject implements Loggable{

    public Input<SMCTreeDensity> treeDensityInput = new Input<>("treeDensity",
            "SMC Tree density from which to log trajectories.",
            Input.Validate.REQUIRED);

    public Input<Boolean> logMostRecentTrajectoryInput = new Input<>("logMostRecentTrajectory",
            "Force density recalculation and logs most recent trajectory. Do not use in a standard " +
                    "PMMH analysis!  Defaults to false.", false);

    SMCTreeDensity treeDensity;
    EpidemicModel model;

    public TrajectoryLogger() { }

    @Override
    public void initAndValidate() {
        treeDensity = treeDensityInput.get();
        model = treeDensity.modelInput.get();
    }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("trajectory\t");
        else
            out.print(getID() + "\t");
    }

    @Override
    public void log(long nSample, PrintStream out) {

        EpidemicTrajectory traj;
        if (logMostRecentTrajectoryInput.get()) {
            treeDensity.calculateLogP();
            traj = treeDensity.getMostRecentTrajectory();
        } else {
            traj = treeDensity.getConditionedTrajectory();
        }

        if (traj == null || traj.getStateList().isEmpty()) {
            out.print("NA\t");
            return;
        }

        boolean isFirst = true;
        for (EpidemicState state : traj.getStateList()) {
            if (!isFirst)
                out.print(",");
            else
                isFirst = false;

            out.print(treeDensity.getModel().getOrigin() - state.time);

            out.print(":" + state.S + ":" + state.I + ":" + state.R + ":" + state.algorithm);

            model.calculatePropensities(state);
            double reproductiveNumber = model.propensities[EpidemicEvent.INFECTION]/
                    (model.propensities[EpidemicEvent.RECOVERY] + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]);

            double incidence = model.propensities[EpidemicEvent.INFECTION];

            out.print(":" + incidence + ":" + reproductiveNumber);
        }

        out.print("\t");
    }


    @Override
    public void close(PrintStream out) { }
}
