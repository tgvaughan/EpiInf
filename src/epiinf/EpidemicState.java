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

/**
 * A state of an SIR epidemic trajectory.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpidemicState {
    public double[] S, I, R;
    int dim;
    
    public EpidemicState() {
        dim = 1;
        S = new double[1];
        I = new double[1];
        R = new double[1];
    }
    
    public EpidemicState(double S, double I, double R) {
        this.S = new double[1];
        this.S[0] = S;
        this.I = new double[1];
        this.I[0] = I;
        this.R = new double[1];
        this.R[0] = R;
    }
    
    /**
     * Test whether state is valid or not.
     * 
     * @return true if state is valid.
     */
    public boolean isValid() {
        for (int i=0; i<dim; i++) 
            if (this.S[i]>=0 && this.I[i]>=0 && this.R[i]>=0)
                return false;
        
        return true;
    }
    
    public EpidemicState copy() {
        EpidemicState stateCopy = new EpidemicState();
        stateCopy.dim = dim;
        stateCopy.S = new double[dim];
        stateCopy.I = new double[dim];
        stateCopy.R = new double[dim];
        
        for (int i=0; i<dim; i++) {
            stateCopy.S[i] = S[i];
            stateCopy.I[i] = I[i];
            stateCopy.R[i] = R[i];
        }
        
        return stateCopy;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        
        for (int i=0; i<dim; i++) {
            if (i>0)
                sb.append(", ");
            sb.append("S[").append(i).append("]:").append(S[i])
                    .append("I[").append(i).append("]:").append(I[i])
                    .append("R[").append(i).append("]:").append(R[i]);
        }

        return sb.toString();
    }
    
    /**
     * Retrieve header suitable for dumping rows of states to an R-compatible
     * file.
     * 
     * @return R input file header
     */
    public String getHeader() {
        StringBuilder sb = new StringBuilder();
        
        for (int i=0; i<dim; i++) {
            if (i>0)
                sb.append(" ");
            sb.append("S").append(i)
                    .append(" I").append(i)
                    .append(" R").append(i);
        }
        
        return sb.toString();
    }
    
    /**
     * Obtain this state as a record in an R-compatible text file.
     * 
     * @return R input file record.
     */
    public String getRecord() {
        StringBuilder sb = new StringBuilder();
        
        for (int i=0; i<dim; i++) {
            if (i>0)
                sb.append(" ");
            sb.append(S[i])
                    .append(" ").append(I[i])
                    .append(" ").append(R[i]);
        }
        
        return sb.toString();
    }
}
