/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import epiinf.models.EpidemicModel;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TransmissionTreeSimulator extends beast.core.Runnable {

    public Input<EpidemicTrajectory> trajInput = new Input<>(
            "epidemicTrajectory",
            "Epidemic trajectory object.",
            Validate.REQUIRED);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName", "Name of file to save Newick representation of tree to.",
            Validate.REQUIRED);

    public Input<Boolean> truncateTrajectoryInput = new Input<>(
            "truncateTrajectory",
            "Truncate trajectory at most recent sample. (Default true.)",
            false);

    public Input<RealParameter> finalTreeSampleOffsetInput = new Input<>(
            "finalTreeSampleOffsetParam",
            "Parameter in which to store final tree sample offset.",
            Validate.REQUIRED);

    public Input<Double> leafSampleFracInput = new Input<>(
            "leafSampleFrac",
            "Fraction of samples that correspond to tree leaves. " +
                    "(The rest are simply used as incidence data.)", 1.0);

    public Input<Boolean> deterministicLeafSampleSelectionInput = new Input<>(
            "deterministicLeafSampleSelection",
            "Deterministically select which samples to associate with " +
                    "tree leaves, instead of probabilistically.", false);

    public Input<Boolean> ensureFinalSampleIsLeafInput = new Input<>(
            "ensureFinalSampleIsLeaf",
            "Ensure the final sample is a tree leaf, even when " +
                    "leafSampleFrac=0.", false);

    public Input<Boolean> measureOriginFromFinalSampleInput = new Input<>(
            "measureOriginFromFinalSample",
            "Adjusts trajectory origin so that age of final sample " +
                    "is zero.  This is for comparing with other likelihoods " +
                    "that do not allow for final sample age to be random.",
            false);

    public Input<RealParameter> incidenceParamInput = new Input<>(
            "incidenceParam",
            "Parameter containing incidence event ages prior to end " +
                    "of observation period.");

    public Input<String> incidenceFileNameInput = new Input<>(
            "incidenceFileName",
            "Name of file to write incidence times to.");
    
    @Override
    public void initAndValidate() { }
    
    @Override
    public void run() {
        
        (new SimulatedTransmissionTree()).initByName(
                "epidemicTrajectory", trajInput.get(),
                "finalTreeSampleOffsetParam", new RealParameter("0.0"),
                "fileName", fileNameInput.get(),
                "truncateTrajectory", truncateTrajectoryInput.get(),
                "finalTreeSampleOffsetParam", finalTreeSampleOffsetInput.get(),
                "leafSampleFrac", leafSampleFracInput.get(),
                "deterministicLeafSampleSelection", deterministicLeafSampleSelectionInput.get(),
                "ensureFinalSampleIsLeaf", ensureFinalSampleIsLeafInput.get(),
                "measureOriginFromFinalSample", measureOriginFromFinalSampleInput.get(),
                "incidenceParam", incidenceParamInput.get(),
                "incidenceFileName", incidenceFileNameInput.get());

    }

    
}
