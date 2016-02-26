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
import beast.core.parameter.RealParameter;

import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class OriginFromTrajectory extends RealParameter {
    public Input<EpidemicTrajectory> trajectoryInput = new Input<>("trajectory",
            "Epidemic trajectory from which to compute origin.",
            Input.Validate.REQUIRED);

    public OriginFromTrajectory() {
        valuesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() throws Exception {
        List<EpidemicEvent> eventList = trajectoryInput.get().getEventList();
        for (int i=0; i<eventList.size(); i++) {
            if (eventList.get(i).isSample()) {
                valuesInput.setValue(eventList.get(i).time, this);
                super.initAndValidate();
                return;
            }
        }

        throw new IllegalArgumentException("Origin undefined for epidemic " +
                "trajectories without a sampling event.");
    }
}
