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
import epiinf.util.IncidenceData;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Maintains a sorted list of observed events.")
public class ObservedEventsList {

    private TreeInterface tree;
    private IncidenceData incidenceData;
    private Function treeOrigin;
    private List<ObservedEvent> eventList;

    /**
     * Tolerance for deviation between tree node ages and trajectory events.
     */
    public static final double tolerance = 1e-10;
    
    private boolean dirty;
    
    /**
     * Default constructor.
     */
    public ObservedEventsList(TreeInterface tree, IncidenceData incidenceData,
                              Function treeOrigin) {
        this.tree = tree;
        this.treeOrigin = treeOrigin;
        this.incidenceData = incidenceData;

        eventList = new ArrayList<>();

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

            if (node.isFake())
                continue;

            ObservedEvent event = new ObservedEvent();
            if (node.isLeaf()) {
                if (node.isDirectAncestor())
                    event.type = ObservedEvent.Type.SAMPLED_ANCESTOR;
                else
                    event.type = ObservedEvent.Type.LEAF;
            } else {
                event.type = ObservedEvent.Type.COALESCENCE;
            }
       
            event.time = getTimeFromHeight(node.getHeight());
            
            eventList.add(event);
        }

        if (incidenceData != null) {
            for (int i = 0; i < incidenceData.nCounts; i++) {
                ObservedEvent event = new ObservedEvent();
                event.type = ObservedEvent.Type.UNSEQUENCED_SAMPLE;
                event.time = getTimeFromHeight(incidenceData.ages[i]);
                event.multiplicity = incidenceData.counts[i];
                eventList.add(event);
            }
        }

        // Sort events in order of absolute time
        Collections.sort(eventList);

        // Collate concurrent events
        int i=1;
        while (i<eventList.size()) {
            if (Math.abs(eventList.get(i).time-eventList.get(i-1).time)<tolerance
                    && eventList.get(i).type == eventList.get(i-1).type) {
                eventList.get(i-1).multiplicity += 1;
                eventList.remove(i);
            } else
                i += 1;
        }

        // Mark final event:
        eventList.get(eventList.size()-1).isFinal = true;

        // Compute lineage counts
        int k = 1;
        for (ObservedEvent event : eventList) {
            event.lineages = k;

            switch (event.type) {
                case COALESCENCE:
                    k += event.multiplicity;
                    break;
                case LEAF:
                    k -= event.multiplicity;
                    break;
                default:
                    break;
            }
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
    public List<ObservedEvent> getEventList() {
        updateEventList();
        return eventList;
    }

    public ObservedEvent getNextObservedEvent(EpidemicState state) {
        updateEventList();

        if (state.observedEventIdx<eventList.size())
            return eventList.get(state.observedEventIdx);
        else
            return null;
    }

    public double getNextObservedEventTime(EpidemicState state) {
        updateEventList();

        if (state.observedEventIdx<eventList.size())
            return eventList.get(state.observedEventIdx).time;
        else
            return Double.POSITIVE_INFINITY;
    }

    public int getCurrentLineageCount(EpidemicState state) {
        updateEventList();

        if (state.observedEventIdx<eventList.size())
            return eventList.get(state.observedEventIdx).lineages;
        else
            return 0;
    }

    /**
     * Cause the event list to be updated before it is next queried
     */
    public void makeDirty() {
        dirty = true;
    }
}
