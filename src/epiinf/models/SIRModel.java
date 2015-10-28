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

package epiinf.models;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;

/**
 * General stochastic SIR model of an epidemic.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SIRModel extends EpidemicModel {
    
    public Input<IntegerParameter> S0Input = new Input<>(
            "S0", "Initial size of susceptible population.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<>(
            "infectionRate", "Infection rate.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<>(
            "recoveryRate", "Recovery rate.", Validate.REQUIRED);

    
    @Override
    public EpidemicState getInitialState() {
        return new EpidemicState(S0Input.get().getValue(), 1, 0);
    }
    
    @Override
    public void calculatePropensities(EpidemicState state) {
        propensities.put(EpidemicEvent.Type.INFECTION,
                infectionRateInput.get().getValue()*state.S*state.I);

        propensities.put(EpidemicEvent.Type.RECOVERY,
                recoveryRateInput.get().getValue()*state.I);

        totalPropensity = propensities.get(EpidemicEvent.Type.INFECTION)
                + propensities.get(EpidemicEvent.Type.RECOVERY);
    }

    @Override
    public void incrementState(EpidemicState state, EpidemicEvent event) {
        switch(event.type) {
            case INFECTION:
                state.S -= event.multiplicity;
                state.I += event.multiplicity;
                break;
            case RECOVERY:
            case RHO_SAMPLE:
            case PSI_SAMPLE:
            case OTHER_SAMPLE:
                state.I -= event.multiplicity;
                state.R += event.multiplicity;
                break;
            default:
                break;
        }
    }
}
