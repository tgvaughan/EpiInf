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
import beast.core.Function;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Maintains a sorted list of events in the tree.")
public class TreeEventList {

    private TreeInterface tree;
    private Function treeOrigin;
    private List<TreeEvent> eventList;
    private List<Integer> lineageCounts;

    /**
     * Tolerance for deviation between tree node ages and trajectory events.
     */
    public static final double tolerance = 1e-10;
    
    private boolean dirty;
    
    /**
     * Default constructor.
     */
    public TreeEventList(TreeInterface tree, Function treeOrigin) {
        this.tree = tree;
        this.treeOrigin = treeOrigin;

        eventList = new ArrayList<>();
        lineageCounts = new ArrayList<>();

        dirty = true;
    }

    /**
     * Ensure list of tree events is up to date.
     */
    public void updateEventList() {
        if (!dirty)
            return;

        eventList.clear();

        // Assemble event list
        for (Node node : tree.getNodesAsArray()) {

            if (node.isDirectAncestor()
                || (!node.isLeaf() &&
                    (node.getLeft().isDirectAncestor()
                    || node.getRight().isDirectAncestor())))
                continue;


            TreeEvent event = new TreeEvent();
            if (node.isLeaf()) {
                if (node.isDirectAncestor())
                    event.type = TreeEvent.Type.SAMPLED_ANCESTOR;
                else
                    event.type = TreeEvent.Type.LEAF;
            } else {
                if (!node.isFake())
                    event.type = TreeEvent.Type.COALESCENCE;
            }
       
            event.time = getTimeFromHeight(node.getHeight());
            
            eventList.add(event);
        }

        // Sort events in order of absolute time
        Collections.sort(eventList, (TreeEvent e1, TreeEvent e2) -> {
            if (e1.time < e2.time)
                return -1;
            
            if (e1.time > e2.time)
                return 1;
            
            return 0;
        });
        
        // Collate concurrent events
        for (int i=1; i<eventList.size(); i++) {
            if (Math.abs(eventList.get(i).time-eventList.get(i-1).time)<tolerance
                    && eventList.get(i).type == eventList.get(i-1).type) {
                eventList.get(i-1).multiplicity += 1;
                eventList.remove(i);
                i -= 1;
            }
        }

        // Compute lineage counts
        lineageCounts.clear();
        int k = 0;
        for (TreeEvent event : eventList) {
            switch (event.type) {
                case COALESCENCE:
                    k += event.multiplicity;
                    break;
                case LEAF:
                    k -= event.multiplicity;
                    break;
                default:
            }

            lineageCounts.add(k);
        }

        dirty = false;
    }
    
    /**
     * Obtain absolute epidemic time corresponding to height on tree.
     * 
     * @param height height to convert
     * @return time
     */
    public double getTimeFromHeight (double height) {
        return treeOrigin.getArrayValue() - height;
    }

    /**
     * @return The time before the most recent sample that the epidemic began.
     */
    public double getOrigin() {
        return treeOrigin.getArrayValue();
    }
    
    /**
     * Retrieve list of events on tree.
     * 
     * @return event list
     */
    public List<TreeEvent> getEventList() {
        updateEventList();
        return eventList;
    }

    /**
     * @return list of lineage counts following each tree event
     */
    public List<Integer> getLineageCounts() {
        updateEventList();
        return lineageCounts;
    }

    /**
     * Cause the event list to be updated before it is next queried
     */
    public void makeDirty() {
        dirty = true;
    }
}
