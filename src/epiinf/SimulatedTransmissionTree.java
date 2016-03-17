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

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 * State node initializer which simulates a transmission tree compatible with
 * the given epidemic trajectory.  The tree origin parameter is also initialized.
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
            "fileName", "Name of file to save Newick representation of tree to.");
    
    public Input<Boolean> truncateTrajectoryInput = new Input<>(
            "truncateTrajectory",
            "Truncate trajectory at most recent sample. (Default true.)", true);

    public Input<Double> sequencingFractionInput = new Input<>(
            "sequencingFraction", "Fraction of samples that are sequenced. " +
            "(The rest are simply used as incidence data.)", 1.0);

    public Input<RealParameter> incidenceParamInput = new Input<>(
            "incidenceParam", "Incidence times parameter");

    public Input<String> incidenceFileNameInput = new Input<>(
            "incidenceFileName", "Name of file to write incidence times to.");

    public SimulatedTransmissionTree() { }
    
    @Override
    public void initAndValidate() {
        EpidemicTrajectory traj = trajInput.get();
        boolean truncateTrajectory = truncateTrajectoryInput.get();
        double seqFrac = sequencingFractionInput.get();

        if (seqFrac<1.0 && incidenceParamInput.get() == null)
            throw new IllegalArgumentException("Must provide incidenceParam " +
                    "if sequencingFraction below 1.0.");

        TreeSet<EpidemicEvent> sequencedSamples = new TreeSet<>();
        TreeSet<EpidemicEvent> unsequencedSamples = new TreeSet<>();

        for (EpidemicEvent event : traj.getEventList()) {
            if (event.isSample()) {

                if (event.type == EpidemicEvent.RHO_SAMPLE
                        || event.type == EpidemicEvent.OTHER_SAMPLE
                        || seqFrac == 1.0 || Randomizer.nextDouble()<seqFrac) {
                    sequencedSamples.add(event);
                } else {
                    unsequencedSamples.add(event);
                }
            }
        }

        if (sequencedSamples.isEmpty())
            throw new NoSamplesException();

        double youngestSequencedSampTime = sequencedSamples.last().time;

        // Store times of unsequenced samples to incidence parameter
        if (incidenceParamInput.get() != null && seqFrac < 1.0) {
            Double[] incidenceTimes = new Double[unsequencedSamples.size()];
            int idx = 0;
            for (EpidemicEvent event : unsequencedSamples)
                incidenceTimes[idx++] = youngestSequencedSampTime - event.time;

            incidenceParamInput.get().assignFromWithoutID(new RealParameter(incidenceTimes));
        }

        // Simulate tree

        int nextLeafNr = 0;
        int nextInternalID = sequencedSamples.size();
        
        List<Node> activeNodes = new ArrayList<>();

        List<EpidemicEvent> revEventList = traj.getEventList();
        Collections.reverse(revEventList);

        List<EpidemicState> revStateList = traj.getStateList();
        Collections.reverse(revStateList);
        
        int sequencedSamplesSeen = 0;
        for (int eidx=0; eidx<revEventList.size(); eidx++) {
            
            EpidemicEvent epidemicEvent = revEventList.get(eidx);
            EpidemicState epidemicState = revStateList.get(eidx);

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
                        if (!sequencedSamples.contains(epidemicEvent))
                            break;

                        leaf = new Node();
                        leaf.setHeight(youngestSequencedSampTime - epidemicEvent.time);
                        leaf.setNr(nextLeafNr);
                        leaf.setID("t" + nextLeafNr);
                        nextLeafNr += 1;
                        activeNodes.add(leaf);

                        sequencedSamplesSeen += 1;
                        break;

                    case EpidemicEvent.PSI_SAMPLE_NOREMOVE:
                        if (!sequencedSamples.contains(epidemicEvent))
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

                        sequencedSamplesSeen += 1;
                        break;

                    default:
                }
            }

            // Stop when we reach the MRCA of all sampled events.
            if (sequencedSamplesSeen==sequencedSamples.size() && activeNodes.size() == 1)
                break;
        }
        
        // Truncate trajectory at most recent sample if requested:
        if (truncateTrajectory) {
            while (revEventList.get(0).time>youngestSequencedSampTime) {
                revEventList.remove(0);
                revStateList.remove(0);
            }
        }

        // Initialise state nodes
        assignFromWithoutID(new Tree(activeNodes.get(0)));
        
        // Write tree to disk if requested:
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                String newick = root.toNewick().concat(";");
                ps.println(newick.replace(":0.0;", ":" + (youngestSequencedSampTime-root.getHeight()) + ";"));
            } catch(FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + fileNameInput.get() + ".");
            }
        }

        // Write incidence times to disk if requested
        if (incidenceFileNameInput.get() != null && !unsequencedSamples.isEmpty()) {
            try (PrintStream ps = new PrintStream(incidenceFileNameInput.get())) {
                ps.println("time age");
                for (EpidemicEvent event : unsequencedSamples) {
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
