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

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Maintains a sorted list of events in the tree.")
public class TreeEventList extends CalculationNode {
    
    public Input<Tree> treeInput = new Input<>("tree",
            "Transmission tree.", Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<>(
            "treeOrigin", "Age of the epidemic measured relative to the"
                    + "most recent sample on the tree.", Validate.REQUIRED);
    
    private Tree tree;
    private RealParameter treeOrigin;
    private List<TreeEvent> eventList, eventListStored;

    
    /**
     * Tolerance for deviation between tree node ages and trajectory events.
     */
    public static final double tolerance = 1e-10;
    
    private boolean dirty;
    
    /**
     * Default constructor.
     */
    public TreeEventList() { }
    
    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        treeOrigin = treeOriginInput.get();
        
        eventList = Lists.newArrayList();
        eventListStored = Lists.newArrayList();
        
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

        dirty = false;
    }
    
    /**
     * Obtain absolute epidemic time corresponding to height on tree.
     * 
     * @param height
     * @return time
     */
    public double getTimeFromHeight (double height) {
        return treeOrigin.getValue() - height;
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
    
    @Override
    protected boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    protected void store() {
        eventListStored.clear();
        eventListStored.addAll(eventList);
        
        super.store();
    }

    @Override
    protected void restore() {
        eventList.clear();
        eventList.addAll(eventListStored);
        dirty = false;
        
        super.restore();
    }
}
