/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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

package beast.app.beauti;

import epiinf.EpidemicState;
import epiinf.SimulatedTrajectory;
import epiinf.models.EpidemicModel;
import org.knowm.xchart.*;
import org.knowm.xchart.internal.chartpart.Chart;
import org.knowm.xchart.internal.style.Styler;
import org.knowm.xchart.internal.style.markers.SeriesMarkers;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiTrajPanel extends JPanel {

    EpidemicModel epidemicModel;

    Chart_XY chart;

    List<List<Number>> times, prevalences;
    int nTraj;

    public EpiTrajPanel(EpidemicModel epidemicModel, int nTraj) {
        this.epidemicModel = epidemicModel;
        this.nTraj = nTraj;

        // Set up chart

        ChartBuilder_XY builder = new ChartBuilder_XY();
        builder.title("Example prevalence trajectories (tau leaping approx.)");
        builder.xAxisTitle("Time before most recent sample");
        builder.yAxisTitle("Prevalence");
        builder.height(getFontMetrics(getFont()).getHeight()*20);
        builder.width(getFontMetrics(getFont()).getHeight()*40);

        chart = builder.build();

        Styler styler = chart.getStyler();
        styler.setChartBackgroundColor(getBackground());
        styler.setLegendVisible(false);

        // Initialise arrays used in simulation.
        times = new ArrayList<>();
        prevalences = new ArrayList<>();
        for (int i=0; i<nTraj; i++) {
            times.add(new ArrayList<>());
            prevalences.add(new ArrayList<>());
        }
    }

    public void updateChart() {
        new TrajWorker().execute();
    }

    class TrajWorker extends SwingWorker<Void, Void> {

        @Override
        protected Void doInBackground() throws Exception {
            double origin = epidemicModel.treeOriginInput.get().getArrayValue();

            final int nSamples = 1000;

            for (int i=0; i<nTraj; i++) {
                List<Number> theseTimes = times.get(i);
                List<Number> thesePrevs = prevalences.get(i);

                SimulatedTrajectory traj = new SimulatedTrajectory(epidemicModel, origin, nSamples);

                theseTimes.clear();
                thesePrevs.clear();

                int sIdx = 0;
                for (EpidemicState state : traj.getStateList()) {
                    theseTimes.add(origin - state.time);
                    thesePrevs.add(state.I);
                }
            }

            return null;
        }

        @Override
        protected void done() {
            for (int i=0; i<nTraj; i++) {
                Series_XY series = chart.addSeries("traj" + i,
                        times.get(i), prevalences.get(i));
                series.setMarker(SeriesMarkers.NONE);
            }

            repaint();
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        chart.paint((Graphics2D)g);
    }

    @Override
    public Dimension getPreferredSize() {
        return new Dimension(chart.getWidth(), chart.getHeight());
    }
}
