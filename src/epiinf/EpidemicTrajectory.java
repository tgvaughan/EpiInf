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
    protected EpidemicState initialState, storedInitialState;
    protected List<EpidemicState> stateList;
    
    protected boolean stateListDirty;
    
    public EpidemicTrajectory () { };
    
    @Override
    public void initAndValidate() {
        eventList = new ArrayList<EpidemicEvent>();
        storedEventList = new ArrayList<EpidemicEvent>();

        stateList = new ArrayList<EpidemicState>();
        stateListDirty = false;
    }
    
    private void updateStateList() {
        if (!stateListDirty)
            return;
        
        stateList.clear();
        
        stateList.add(initialState);
        
        for (EpidemicEvent event : eventList) {
            EpidemicState nextState = stateList.get(stateList.size()-1).copy();
            
            if (event.type == EpidemicEvent.EventType.INFECTION) {
                nextState.S -= 1; 
                nextState.I += 1;
            } else {
                nextState.I -= 1;
                nextState.R += 1;
            }
            
            stateList.add(nextState);
        }
    }
    
    /**
     * Retrieve list of states corresponding to complete epidemic trajectory.
     * 
     * @return state list
     */
    public List<EpidemicState> getStateList() {
        updateStateList();
        return stateList;
    }
    
    /**
     * Set the initial state of the epidemic.
     * 
     * @param state 
     */
    public void setInitialState(EpidemicState state) {
        startEditing();
        initialState = state;
        stateListDirty = true;
    }
    
    /**
     * Retrieve the initial state of the epidemic.
     * 
     * @return initial state
     */
    public EpidemicState getInitialState() {
        return initialState;
    }

    public void setEventList(List<EpidemicEvent> eventList) {
        startEditing();
        this.eventList = eventList;
        stateListDirty = true;
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
     * Write trajectory state sequence to PrintStream.
     * 
     * @param ps where to send output
     */
    public void dumpTrajectory(PrintStream ps) {
        ps.println("t S I R");
        ps.format("0.0 %d %d %d\n",
                initialState.S,
                initialState.I,
                initialState.R);
        
        for (int i=0; i<eventList.size(); i++) {
            double t = eventList.get(i).time;
            int S = getStateList().get(i+1).S;
            int I = getStateList().get(i+1).I;
            int R = getStateList().get(i+1).R;
            
            ps.format("%g %d %d %d\n", t, S, I, R);
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
        trajCopy.initialState = initialState.copy();
        trajCopy.eventList.addAll(eventList);
        trajCopy.storedEventList.addAll(storedEventList);
        trajCopy.stateListDirty = true;
        
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
        storedInitialState = initialState.copy();
        storedEventList.clear();
        storedEventList.addAll(eventList);
    }

    @Override
    public void restore() {
        initialState = storedInitialState.copy();
        eventList.clear();
        eventList.addAll(storedEventList);
        stateListDirty = true;
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
        int Ipeak = 0;
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
