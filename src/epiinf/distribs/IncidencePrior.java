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

package epiinf.distribs;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import epiinf.IncidenceData;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IncidencePrior extends Distribution {

    public Input<IncidenceData> incidenceDataInput = new Input<>(
            "incidenceData", "Incidence data.", Input.Validate.REQUIRED);

    public Input<RealParameter> incidenceParamInput = new Input<>(
            "incidenceParam", "Incidence parameter", Input.Validate.REQUIRED);

    List<Double> lower, upper;

    IncidenceData incidenceData;
    RealParameter incidenceParameter;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        incidenceData = incidenceDataInput.get();
        incidenceParameter = incidenceParamInput.get();

        lower = new ArrayList<>();
        upper = new ArrayList<>();
        for (int i=0; i<incidenceData.getAges().size(); i++) {
            lower.add(incidenceData.getAges().get(i) - incidenceData.getError());
            upper.add(incidenceData.getAges().get(i) + incidenceData.getError());
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int i=0; i<incidenceData.getAges().size(); i++) {
            if (incidenceParameter.getValue(i) < lower.get(i)
                    || incidenceParameter.getValue(i) > upper.get(i)) {
                logP = Double.NEGATIVE_INFINITY;
                break;
            }
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Sampling not supported.");
    }
}
