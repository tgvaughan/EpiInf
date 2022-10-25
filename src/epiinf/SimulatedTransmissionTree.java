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

package epiinf;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 * State node initializer which simulates a transmission tree compatible with
 * the given epidemic trajectory.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate transmission tree conditional on epidemic trajectory.")
public class SimulatedTransmissionTree extends Tree {

    public Input<EpidemicTrajectory> trajInput = new Input<>(
            "epidemicTrajectory",
            "Epidemic trajectory object.",
            Validate.REQUIRED);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to save Newick representation of tree to.");
    
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

    public SimulatedTransmissionTree() { }
    
    @Override
    public void initAndValidate() {
        EpidemicTrajectory traj = trajInput.get();
        boolean truncateTrajectory = truncateTrajectoryInput.get();
        double leafFrac = leafSampleFracInput.get();
        boolean useDetLeafSel = deterministicLeafSampleSelectionInput.get();
        boolean measureOriginFromFinalSample = measureOriginFromFinalSampleInput.get();

        if (leafFrac<1.0 && incidenceParamInput.get() == null)
            throw new IllegalArgumentException("Must provide incidenceParam " +
                    "if leafSampleFrac below 1.0.");

        TreeSet<EpidemicEvent> sequencedSamplingEvents = new TreeSet<>();
        TreeSet<EpidemicEvent> unsequencedSamplingEvents = new TreeSet<>();

        int nLeafSamples = 0, nNonleafSamples = 0;
        double cumulativeLeafFrac = 0;

        for (EpidemicEvent event : traj.getEventList()) {
            if (event.isSample()) {

                if (event.type == EpidemicEvent.RHO_SAMPLE
                        || event.type == EpidemicEvent.OTHER_SAMPLE
                        || leafFrac == 1.0
                        || (useDetLeafSel && cumulativeLeafFrac < leafFrac)
                        || Randomizer.nextDouble()<leafFrac) {
                    sequencedSamplingEvents.add(event);
                    nLeafSamples += event.multiplicity;
                } else {
                    unsequencedSamplingEvents.add(event);
                    nNonleafSamples += event.multiplicity;
                }
                cumulativeLeafFrac = nLeafSamples/(float)(nLeafSamples + nNonleafSamples);
            }
        }

        if (ensureFinalSampleIsLeafInput.get() && nNonleafSamples>0) {
            if (nLeafSamples == 0) {
                EpidemicEvent lastUnseqEvent = unsequencedSamplingEvents.last();
                unsequencedSamplingEvents.remove(lastUnseqEvent);
                sequencedSamplingEvents.add(lastUnseqEvent);
                nNonleafSamples -= lastUnseqEvent.multiplicity;
                nLeafSamples += lastUnseqEvent.multiplicity;

            } else {

                if (sequencedSamplingEvents.last().time < unsequencedSamplingEvents.last().time) {
                    EpidemicEvent lastSeqEvent = sequencedSamplingEvents.last();
                    EpidemicEvent lastUnseqEvent = unsequencedSamplingEvents.last();

                    sequencedSamplingEvents.remove(lastSeqEvent);
                    unsequencedSamplingEvents.remove(lastUnseqEvent);
                    sequencedSamplingEvents.add(lastUnseqEvent);
                    unsequencedSamplingEvents.add(lastSeqEvent);

                    int leafCountDelta = lastUnseqEvent.multiplicity - lastSeqEvent.multiplicity;
                    nLeafSamples += leafCountDelta;
                    nNonleafSamples -= leafCountDelta;
                }
            }
        }

        if (nLeafSamples + nNonleafSamples == 0)
            throw new NoSamplesException();

        // Adjust origin if requested
        if (measureOriginFromFinalSample) {
            double finalSampleTime = 0;
            if (nLeafSamples > 0)
                finalSampleTime = Math.max(finalSampleTime, sequencedSamplingEvents.last().time);
            if (nNonleafSamples > 0)
                finalSampleTime = Math.max(finalSampleTime, unsequencedSamplingEvents.last().time);

            traj.origin = finalSampleTime;
        }


        // Store ages (relative to end of process) of unsequenced samples
        // to incidence parameter
        if (incidenceParamInput.get() != null && leafFrac < 1.0) {
            Double[] incidenceTimes = new Double[unsequencedSamplingEvents.size()];
            int idx = 0;
            for (EpidemicEvent event : unsequencedSamplingEvents)
                incidenceTimes[idx++] = traj.origin - event.time;

            incidenceParamInput.get().assignFromWithoutID(new RealParameter(incidenceTimes));
        }

        // Simulate tree

        double youngestSequencedSampTime = !sequencedSamplingEvents.isEmpty()
                ? sequencedSamplingEvents.last().time
                : Double.NaN;

        int nextLeafNr = 0;
        int nextInternalID = nLeafSamples;

        List<Node> activeNodes = new ArrayList<>();

        List<EpidemicEvent> revEventList = traj.getEventList();
        Collections.reverse(revEventList);

        List<EpidemicState> revStateList = traj.getStateList();
        Collections.reverse(revStateList);
        
        int nSequencedSamplingEventsSeen = 0;
        for (int eidx=0; eidx<revEventList.size(); eidx++) {
            
            EpidemicEvent epidemicEvent = revEventList.get(eidx);
            EpidemicState epidemicState = revStateList.get(eidx);

            if (sequencedSamplingEvents.contains(epidemicEvent))
                nSequencedSamplingEventsSeen += 1;

            int k = activeNodes.size();
            double N = epidemicState.I;

            for (int i=0; i<epidemicEvent.multiplicity; i++) {

                Node leaf, parent;

                switch (epidemicEvent.type) {
                    case EpidemicEvent.INFECTION:
                        double pCoalesce = k * (k - 1) / (N * (N - 1));

                        if (Randomizer.nextDouble() < pCoalesce) {
                            int childIdx = Randomizer.nextInt(k);
                            Node child1 = activeNodes.get(childIdx);
                            activeNodes.remove(childIdx);

                            childIdx = Randomizer.nextInt(k - 1);
                            Node child2 = activeNodes.get(childIdx);
                            activeNodes.remove(childIdx);

                            parent = new Node();
                            parent.addChild(child1);
                            parent.addChild(child2);
                            parent.setHeight(youngestSequencedSampTime - epidemicEvent.time);
                            parent.setNr(nextInternalID++);
                            activeNodes.add(parent);
                        }
                        break;

                    case EpidemicEvent.RHO_SAMPLE:
                    case EpidemicEvent.PSI_SAMPLE_REMOVE:
                        if (!sequencedSamplingEvents.contains(epidemicEvent))
                            break;

                        leaf = new Node();
                        leaf.setHeight(youngestSequencedSampTime - epidemicEvent.time);
                        leaf.setNr(nextLeafNr);
                        leaf.setID("t" + nextLeafNr);
                        nextLeafNr += 1;
                        activeNodes.add(leaf);
                        break;

                    case EpidemicEvent.PSI_SAMPLE_NOREMOVE:
                        if (!sequencedSamplingEvents.contains(epidemicEvent))
                            break;

                        leaf = new Node();
                        leaf.setHeight(youngestSequencedSampTime - epidemicEvent.time);
                        leaf.setNr(nextLeafNr);
                        leaf.setID("t" + nextLeafNr);
                        nextLeafNr += 1;

                        double u = Randomizer.nextDouble()*N;
                        if (u < k) {
                            Node sibling = activeNodes.get((int)u);
                            activeNodes.remove(sibling);
                            parent = new Node();
                            parent.addChild(sibling);
                            parent.addChild(leaf);
                            parent.setHeight(youngestSequencedSampTime - epidemicEvent.time);
                            parent.setNr(nextInternalID++);
                            activeNodes.add(parent);
                        } else {
                            activeNodes.add(leaf);
                        }
                        break;

                    default:
                }
            }

            // Stop when we reach the MRCA of all sampled events.
            if (nSequencedSamplingEventsSeen==sequencedSamplingEvents.size() && activeNodes.size() == 1)
                break;
        }

        // Record final tree sample offset
        double offset = traj.origin - revEventList.get(0).time;
        if (nLeafSamples>0) {
            finalTreeSampleOffsetInput.get().assignFromWithoutID(
                    new RealParameter(String.valueOf(offset)));
        }

        // Truncate trajectory at most recent sample if requested:
        if (truncateTrajectory) {
            while (revEventList.get(0).time>youngestSequencedSampTime) {
                revEventList.remove(0);
                revStateList.remove(0);
            }
        }

        // Initialise state nodes
        if (!sequencedSamplingEvents.isEmpty())
            assignFromWithoutID(new Tree(activeNodes.get(0)));

        initArrays();
        
        // Write tree to disk if requested:

        if (!sequencedSamplingEvents.isEmpty() && fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                String newick = root.toNewick().concat(";");
                ps.println(newick.replace(":0.0;", ":" + (youngestSequencedSampTime-root.getHeight()) + ";"));
            } catch(FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + fileNameInput.get() + ".");
            }
        }

        // Write incidence times to disk if requested
        if (incidenceFileNameInput.get() != null && !unsequencedSamplingEvents.isEmpty()) {
            try (PrintStream ps = new PrintStream(incidenceFileNameInput.get())) {
                ps.println("time age");
                for (EpidemicEvent event : unsequencedSamplingEvents) {
                    ps.println(event.time + " " + (youngestSequencedSampTime - event.time));
                }
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Error writing to file "
                        + incidenceFileNameInput.get() + ".");
            }
        }
    }
    
    /**
     * Exception thrown when available trajectory contains no samples to
     * generate a tree from.
     */
    public class NoSamplesException extends RuntimeException {

        @Override
        public String getMessage() {
            return "Tried to simulate tree from trajectory "
                    + "containing no samples.";
        }
    }
}
