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
import beast.core.parameter.RealParameter;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;

/**
 * General birth-death model of an epidemic.  Suitable only for the
 * exponential growth period of the epidemic.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BirthDeathModel extends EpidemicModel {
    
    public Input<RealParameter> birthRateInput = new Input<RealParameter>(
            "birthRate", "Lineage birth rate.", Validate.REQUIRED);
    
    public Input<RealParameter> deathRateInput = new Input<RealParameter>(
            "deathRate", "Lineage death rate.", Validate.REQUIRED);
    
    public Input<RealParameter> durationInput = new Input<RealParameter>(
            "duration", "Parameter specifying duration of process.  Usually"
            + "set to the origin of the tree.", Validate.REQUIRED);

    @Override
    public EpidemicState getInitialState() {
        return new EpidemicState(0, 1, 0);
    }
    
    @Override
    public void calculatePropensities(EpidemicState state) {
        propensities.put(EpidemicEvent.Type.INFECTION,
                birthRateInput.get().getValue()*state.I);

        propensities.put(EpidemicEvent.Type.RECOVERY,
                deathRateInput.get().getValue()*state.I);

        totalPropensity = propensities.get(EpidemicEvent.Type.INFECTION)
                + propensities.get(EpidemicEvent.Type.RECOVERY);
    }

    @Override
    public void incrementState(EpidemicState state, EpidemicEvent.Type type) {
        switch(type) {
            case INFECTION:
                state.I += 1;
                break;
            case RECOVERY:
                state.I -= 1;
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
        return lineages*(lineages-1)/(N*(N-1));
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
