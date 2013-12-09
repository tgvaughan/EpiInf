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
    
    public Input<IntegerParameter> S0Input = new Input<IntegerParameter>(
            "S0", "Initial size of susceptible population.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate", "Infection rate.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
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
    public void incrementState(EpidemicState state, EpidemicEvent.Type type) {
        switch(type) {
            case INFECTION:
                state.S -= 1;
                state.I += 1;
                break;
            case RECOVERY:
                state.I -= 1;
                state.R += 1;
                break;
            default:
                break;
        }
    }

    @Override
    public EpidemicEvent.Type getCoalescenceEventType() {
        return EpidemicEvent.Type.INFECTION;
    }

    @Override
    public EpidemicEvent.Type getLeafEventType() {
        return EpidemicEvent.Type.RECOVERY;
    }

    @Override
    public double getProbCoalescence(EpidemicState state, int lineages) {
        double N = state.I;
        
        if (N<2)
            return 0.0;
        else
            return (lineages)*(lineages-1)/(N*(N-1));
    }

    @Override
    public double getProbNoCoalescence(EpidemicState state, int lineages) {
        return 1.0 - getProbCoalescence(state, lineages);
    }

    @Override
    public double getProbLeaf() {
        return 1.0;
    }

    @Override
    public double getProbNoLeaf() {
        return 1.0;
    }


}
