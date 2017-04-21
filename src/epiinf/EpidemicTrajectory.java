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

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Loggable;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Object representing a complete epidemic trajectory.")
public class EpidemicTrajectory extends BEASTObject implements Loggable {

    protected List<EpidemicEvent> eventList;
    protected List<EpidemicState> stateList;
    protected Double origin;
    
    
    public EpidemicTrajectory () {
        eventList = new ArrayList<>();
        stateList = new ArrayList<>();
    }

    public EpidemicTrajectory (List<EpidemicEvent> eventList,
                               List<EpidemicState> stateList,
                               double origin) {
        this.eventList = eventList;
        this.stateList = stateList;
        this.origin = origin;
    }

    public void assignFrom(EpidemicTrajectory otherTrajectory) {
        this.eventList = otherTrajectory.getEventList();
        this.stateList = otherTrajectory.getStateList();
        this.origin = otherTrajectory.getOrigin();
    }

    
    @Override
    public void initAndValidate() {
    }

    /**
     * Retrieve list of states corresponding to complete epidemic trajectory
     * or null if unavailable.
     * 
     * @return state list
     */
    public List<EpidemicState> getStateList() {
        return stateList;
    }

    /**
     * Retrieve list of events corresponding to complete epidemic trajectory
     * or null if unavailable.
     * 
     * @return event list
     */
    public List<EpidemicEvent> getEventList() {
        return eventList;
    }

    /**
     * Retrieve origin parameter or null if unavailable.
     *
     * @return origin parameter
     */
    public Double getOrigin() {
        return origin;
    }
    
    /**
     * Write trajectory state sequence to PrintStream.
     * 
     * @param ps where to send output
     */
    public void dumpTrajectory(PrintStream ps) {
        ps.println("t " + EpidemicState.getHeader() + " eventType eventMultiplicity");
        ps.println("0.0 " + getStateList().get(0).getRecord() + " START 0");
        for (int i=0; i<eventList.size(); i++) {
            double t = eventList.get(i).time;
            ps.format("%g %s %s\n", t,
                    getStateList().get(i+1).getRecord(),
                    eventList.get(i).getRecord());
        }
    }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("trajectory\t");
        else
            out.print(getID() + "\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        if (stateList.isEmpty()) {
            out.print("NA\t");
            return;
        }

        boolean isFirst = true;
        for (EpidemicState state : stateList) {
            if (!isFirst)
                out.print(",");
            else
                isFirst = false;

            if (origin != null)
                out.print(origin - state.time);
            else
                out.print(state.time);

            out.print(":" + state.S + ":" + state.I + ":" + state.R + ":" + state.algorithm);
        }

        out.print("\t");
    }

    @Override
    public void close(PrintStream out) { }

}
