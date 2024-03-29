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
import beast.base.core.Function;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import epiinf.models.EpidemicModel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Maintains a sorted list of observed events.")
public class ObservedEventsList {

    private TreeInterface tree;
    private Function incidenceAges;
    private IncidenceData incidenceData;
    private Function finalTreeSampleOffset;
    private List<ObservedEvent> eventList, rhoEventsToAdd;

    private EpidemicModel model;

    /**
     * Tolerance for deviation between tree node ages and trajectory events.
     */
    public static final double tolerance = 1e-10;
    
    private boolean dirty;
    
    public ObservedEventsList(TreeInterface tree,
                              Function incidenceAges,
                              IncidenceData incidenceData,
                              EpidemicModel model,
                              Function finalTreeSampleOffset) {
        this.tree = tree;
        this.incidenceAges = incidenceAges;
        this.incidenceData = incidenceData;
        this.model = model;
        this.finalTreeSampleOffset = finalTreeSampleOffset;

        eventList = new ArrayList<>();
        rhoEventsToAdd = new ArrayList<>();

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

        if (tree != null) {
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

                event.time = getTimeFromAge(node.getHeight() + finalTreeSampleOffset.getArrayValue());

                eventList.add(event);
            }
        }

        if (incidenceAges != null) {
            for (int i = 0; i < incidenceAges.getDimension(); i++) {
                ObservedEvent event = new ObservedEvent();
                event.type = ObservedEvent.Type.UNSEQUENCED_SAMPLE;
                event.time = getTimeFromAge(incidenceAges.getArrayValue(i));
                eventList.add(event);
            }
        }

        if (incidenceData != null) {
            List<Double> uniqueAges = incidenceData.getUniqueAges();
            List<Integer> ageCounts = incidenceData.getAgeCounts();
            for (int i=0; i<uniqueAges.size(); i++) {
                ObservedEvent event = new ObservedEvent();
                event.type = ObservedEvent.Type.UNSEQUENCED_SAMPLE;
                event.time = getTimeFromAge(uniqueAges.get(i));
                event.multiplicity = ageCounts.get(i);
                eventList.add(event);
            }
        }

        // Sort events in order of absolute time
        Collections.sort(eventList);

        // Add multiplicity 0 LEAF events for each rho sampling time.
        // Required to ensure a lack of samples resulting from a rho sampling
        // event is still considered an observation.
        eventList.addAll(model.getModelEventList().stream()
                .filter(e -> e.type == ModelEvent.Type.RHO_SAMPLING)
                .map(e -> {
                    ObservedEvent event = new ObservedEvent();
                    event.time = e.time;
                    event.type = ObservedEvent.Type.LEAF;
                    event.multiplicity = 0;

                    return event;
                })
                .collect(Collectors.toList()));

        // Include end-of-observation event
        // (This is always the last event in the list, even when a rho sampling event occurs at
        // exactly the same time.)
        ObservedEvent endOfObservationEvent = new ObservedEvent();
        endOfObservationEvent.type = ObservedEvent.Type.OBSERVATION_END;
        endOfObservationEvent.time = model.getOrigin();
        eventList.add(endOfObservationEvent);

        // Collate concurrent events
        int i=1;
        while (i<eventList.size()) {
            if (Math.abs(eventList.get(i).time-eventList.get(i-1).time)<tolerance
                    && eventList.get(i).type == eventList.get(i-1).type) {
                eventList.get(i-1).multiplicity += eventList.get(i).multiplicity;
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
     * Obtain absolute epidemic time corresponding to given age prior to end of observation period.
     * 
     * @param age age to convert
     * @return time
     */
    public double getTimeFromAge(double age) {
        return model.getOrigin() - age;
    }

    /**
     * @return The time before the end of the observation period that the epidemic began.
     */
    public double getOrigin() {
        return model.getOrigin();
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

    public void dump(PrintStream ps) {
        dump(ps, null);
    }


    public void dump(PrintStream ps, ObservedEvent finalEvent) {
        updateEventList();

        ps.println("time multiplicity type isFinal");
        for (ObservedEvent event : eventList) {
            ps.println(event.time + " " + event.multiplicity + " " + event.type
                    + " " + (event == finalEvent));
        }
    }
}
