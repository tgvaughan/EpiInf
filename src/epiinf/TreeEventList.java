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
import java.util.Comparator;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Maintains a sorted list of events in the tree.")
public class TreeEventList extends CalculationNode {
    
    public Input<Tree> treeInput = new Input<Tree>("tree",
            "Transmission tree.", Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<RealParameter>(
            "treeOrigin", "Difference between time of MRCA and start of "
                    + "epidemic.", Validate.REQUIRED);
    
    public enum TreeEventType { SAMPLE, COALESCENCE }
    public class TreeEvent {
        public Node node;
        public TreeEventType type;
        public double time;    
    }
    
    private Tree tree;
    private RealParameter treeOrigin;
    private List<TreeEvent> eventList, eventListStored;
    
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
        
        updateEventList();
    }
    
    /**
     * Ensure list of tree events is up to date.
     */
    private void updateEventList() {
        eventList.clear();

        // Assemble event list
        for (Node node : tree.getNodesAsArray()) {
            TreeEvent event = new TreeEvent();
            if (node.isLeaf())
                event.type = TreeEventType.SAMPLE;
            else
                event.type = TreeEventType.COALESCENCE;
       
            event.time = getTimeFromHeight(node.getHeight());
            event.node = node;
            
            eventList.add(event);
        }
        
        // Sort events in time-reversed order
        Collections.sort(eventList, new Comparator<TreeEvent>() {

            @Override
            public int compare(TreeEvent e1, TreeEvent e2) {
                if (e1.time < e2.time)
                    return 1;
                
                if (e1.time > e2.time)
                    return -1;
                
                return 0;
            }
            
        });

        dirty = false;
    }
    
    /**
     * Obtain absolute epidemic time corresponding to height on tree.
     * 
     * @param height
     * @return time
     */
    private double getTimeFromHeight (double height) {
        return treeOrigin.getValue() + tree.getRoot().getHeight() - height;
    }

    /**
     * Retrieve list of events on tree.
     * 
     * @return event list
     */
    public List<TreeEvent> getEventList() {
        if (dirty)
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
