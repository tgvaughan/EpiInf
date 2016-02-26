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
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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

    public SimulatedTransmissionTree() { }
    
    @Override
    public void initAndValidate() throws Exception {
        EpidemicTrajectory traj = trajInput.get();
        boolean truncateTrajectory = truncateTrajectoryInput.get();

        int nSamples = 0;
        double youngestSamp = 0.0;
        for (EpidemicEvent event : traj.getEventList()) {
            if (event.isSample()) {
                nSamples += event.multiplicity;
                youngestSamp = Math.max(event.time, youngestSamp);
            }
        }
        
        if (nSamples==0)
            throw new NoSamplesException();

        int nextLeafNr = 0;
        int nextInternalID = nSamples;
        
        List<Node> activeNodes = new ArrayList<>();

        List<EpidemicEvent> revEventList = traj.getEventList();
        Collections.reverse(revEventList);

        List<EpidemicState> revStateList = traj.getStateList();
        Collections.reverse(revStateList);
        
        int samplesSeen = 0;
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
                            parent.setHeight(youngestSamp - epidemicEvent.time);
                            parent.setNr(nextInternalID++);
                            activeNodes.add(parent);
                        }
                        break;

                    case EpidemicEvent.RHO_SAMPLE:
                    case EpidemicEvent.PSI_SAMPLE_REMOVE:
                        leaf = new Node();
                        leaf.setHeight(youngestSamp - epidemicEvent.time);
                        leaf.setNr(nextLeafNr);
                        leaf.setID("t" + nextLeafNr);
                        nextLeafNr += 1;
                        activeNodes.add(leaf);
                        samplesSeen += 1;
                        break;

                    case EpidemicEvent.PSI_SAMPLE_NOREMOVE:
                        leaf = new Node();
                        leaf.setHeight(youngestSamp - epidemicEvent.time);
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
                            parent.setHeight(youngestSamp - epidemicEvent.time);
                            parent.setNr(nextInternalID++);
                            activeNodes.add(parent);
                        } else {
                            activeNodes.add(leaf);
                        }

                        samplesSeen += 1;

                        break;

                    default:
                }
            }

            // Stop when we reach the MRCA of all sampled events.
            if (samplesSeen==nSamples && activeNodes.size()<2)
                break;
        }
        
        // Truncate trajectory at most recent sample if requested:
        if (truncateTrajectory) {
            while (revEventList.get(0).time>youngestSamp) {
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
                ps.println(newick.replace(":0.0;", ":" + (youngestSamp-root.getHeight()) + ";"));
            }
        }
    }
    
    /**
     * Exception thrown when available trajectory contains no samples to
     * generate a tree from.
     */
    public class NoSamplesException extends Exception {

        @Override
        public String getMessage() {
            return "Tried to simulate tree from trajectory "
                    + "containing no samples.";
        }
    }
}
