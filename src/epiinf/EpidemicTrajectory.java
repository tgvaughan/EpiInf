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
import beast.core.StateNode;
import com.google.common.collect.Lists;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import org.w3c.dom.Node;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("StateNode representing a complete epidemic trajectory.")
public class EpidemicTrajectory extends StateNode {
    
    protected List<EpidemicEvent> eventList, storedEventList;
    protected List<EpidemicState> stateList, storedStateList;
    
    
    public EpidemicTrajectory () { };
    
    @Override
    public void initAndValidate() {
        eventList = Lists.newArrayList();
        storedEventList = Lists.newArrayList();

        stateList = Lists.newArrayList();
        storedStateList = Lists.newArrayList();
    }
    
    /**
     * Retrieve list of states corresponding to complete epidemic trajectory.
     * 
     * @return state list
     */
    public List<EpidemicState> getStateList() {
        return stateList;
    }
    
    /**
     * Retrieve the initial state of the epidemic.
     * 
     * @return initial state
     */
    public EpidemicState getInitialState() {
        return stateList.get(0);
    }

    public void setEventAndStateList(List<EpidemicEvent> eventList,
            List<EpidemicState> stateList) {
        startEditing();
        this.eventList = eventList;
        this.stateList = stateList;
    }
    
    /**
     * Retrieve list of events corresponding to complete epidemic trajectory.
     * 
     * @return event list
     */
    public List<EpidemicEvent> getEventList() {
        return eventList;
    }
    
    /**
     * Retrieve length of time between start of epidemic and last recorded
     * event.
     * 
     * @return duration of recorded epidemic
     */
    public double getDuration() {
        return eventList.get(eventList.size()-1).time;
    }
    
    /**
     * Write trajectory state sequence to PrintStream.
     * 
     * @param ps where to send output
     */
    public void dumpTrajectory(PrintStream ps) {
        ps.println("t " + EpidemicState.getHeader());
        ps.println("0.0 " + getStateList().get(0).getRecord());
        for (int i=0; i<eventList.size(); i++) {
            double t = eventList.get(i).time;
            ps.format("%g %s\n", t, getStateList().get(i+1).getRecord());
        }
    }

    /**
     * If state exists, notify that state that this statenode has changed.
     */
    protected void startEditing() {
        if (getState() != null)
            startEditing(null);
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public StateNode copy() {
        EpidemicTrajectory trajCopy = new EpidemicTrajectory();
        
        trajCopy.initAndValidate();
        trajCopy.eventList.addAll(eventList);
        trajCopy.storedEventList.addAll(storedEventList);
        trajCopy.stateList.addAll(stateList);
        trajCopy.storedStateList.addAll(storedStateList);
        
        return trajCopy;
    }

    @Override
    public void assignTo(StateNode other) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void assignFrom(StateNode other) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void assignFromFragile(StateNode other) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void fromXML(Node node) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int scale(double fScale) throws Exception {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    protected void store() {
        storedEventList.clear();
        storedEventList.addAll(eventList);
        storedStateList.clear();
        storedStateList.addAll(stateList);
    }

    @Override
    public void restore() {
        eventList.clear();
        eventList.addAll(storedEventList);
        stateList.clear();
        stateList.addAll(storedStateList);

        hasStartedEditing = false;
    }

    @Override
    public void init(PrintStream out) throws Exception {
        out.print(getID() + ".duration\t");
        out.print(getID() + ".timeAtPeak\t");
        out.print(getID() + ".infected\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(eventList.get(eventList.size()-1).time + "\t");
        
        double tpeak = 0.0;
        double Ipeak = 0;
        for (int idx=0; idx<getEventList().size(); idx++) {
            EpidemicEvent epiEvent = getEventList().get(idx);
            EpidemicState epiState = getStateList().get(idx);
            if (epiState.I>Ipeak) {
                Ipeak = epiState.I;
                tpeak = epiEvent.time;
            }
        }
        out.print(tpeak + "\t");
        
        EpidemicState finalState = getStateList().get(getStateList().size()-1);
        out.print((finalState.I+finalState.R) + "\t");
    }

    @Override
    public void close(PrintStream out) { }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return 0.0;
    }

    @Override
    public double getArrayValue(int iDim) {
        return 0.0;
    }
    
}
