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

package epiinf.inference;

import beast.core.Input;
import beast.core.MCMC;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ProdableMCMC extends MCMC {

    public Input<Integer> sspInput = new Input<>("stepsPerProd",
            "Number of steps between subsequent prods. " +
                    "A value of zero disables prodding.", 0);

    int ssp;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        ssp = sspInput.get();
    }

    @Override
    protected void callUserFunction(int iSample) {

            if (ssp>0 && iSample % ssp == 0) {
                try {
                    oldLogLikelihood = state.robustlyCalcPosterior(posterior);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
    }
}
