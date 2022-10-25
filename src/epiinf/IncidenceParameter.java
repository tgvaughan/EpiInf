/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IncidenceParameter extends RealParameter {

    public Input<IncidenceData> incidenceDataInput = new Input<>(
            "incidenceData", "Incidence data.", Input.Validate.REQUIRED);

    IncidenceData incidenceData;

    public IncidenceParameter() {
        valuesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        incidenceData = incidenceDataInput.get();

        dimensionInput.setValue(String.valueOf(incidenceData.getAges().size()), this);

        StringBuilder valueStringBuilder = new StringBuilder();
        for (double t : incidenceData.getAges()) {
            if (incidenceData.getError()>0.0)
                t += Randomizer.nextDouble()*incidenceData.getError() - 0.5*incidenceData.getError();

            t = Math.max(t, 0.0); // Samples times can't be earlier than the start of the sampling period

            valueStringBuilder.append(" ").append(t);
        }

        valuesInput.setValue(valueStringBuilder.toString(), this);

        super.initAndValidate();
    }


}
