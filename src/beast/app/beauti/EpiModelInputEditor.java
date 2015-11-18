/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import epiinf.models.BirthDeathModel;
import epiinf.models.EpidemicModel;
import epiinf.models.SIRModel;
import epiinf.models.SISModel;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiModelInputEditor extends InputEditor.Base {

    EpidemicModel epidemicModel;
    String[] modelNames = {"Linear birth-death", "SIS", "SIR"};

    IntegerParameter S0;
    RealParameter infectionRate, recoveryRate, psiSamplingRate, rhoSamplingProb;
    RealParameter infectionRateChangeTimes, recoveryRateChangeTimes, psiSamplingRateChangeTimes;

    ComboBoxModel<String> emSelectorModel;
    DefaultTableModel infectionRateModel, infectionRateChangeTimesModel;
    JTable infectionRateTable, infectionRateChangeTimesTable;

    DefaultTableModel recoveryRateModel, recoveryRateChangeTimesModel;
    JTable recoveryRateTable, recoveryRateChangeTimesTable;

    DefaultTableModel psiSamplingRateModel, psiSamplingRateChangeTimesModel;
    JTable psiSamplingRateTable, psiSamplingRateChangeTimesTable;

    SpinnerNumberModel nInfectionRateShiftsModel,
            nRecoveryRateShiftsModel,
            nPsiSamplingRateShiftsModel;
    JSpinner nInfectionRateShiftsSpinner,
            nRecoveryRateShiftsSpinner,
            nPsiSamplingRateShiftsSpinner;

    JTextField S0TextField, rhoSamplingProbTextField;

    JCheckBox estimateS0, estimateInfectionRate, estimateRecoveryRate;
    JCheckBox estimateInfectionRateShiftTimes, estimateRecoveryRateShiftTimes;
    JCheckBox estimatePsiSamplingRate, estimatePsiSamplingRateShiftTimes;
    JCheckBox estimateRhoSamplingProb;

    public EpiModelInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return EpidemicModel.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr,
                     ExpandOption bExpandOption, boolean bAddButtons) {

        // Set up fields
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_plugin = plugin;
        this.itemNr = itemNr;
        epidemicModel = (EpidemicModel) input.get();

        // Adds label to the left of input editor
        addInputLabel();

        // Create component models and fill them with data from input
        emSelectorModel = new DefaultComboBoxModel<>(modelNames);

        // Create and lay out GUI components
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
        panel.setBorder(new EtchedBorder());

        Box box = Box.createHorizontalBox();
        box.add(new JLabel("Epidemic model to use:"));

        JComboBox<String> emSelector = new JComboBox<>(emSelectorModel);
        box.add(emSelector);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Infection rate:"));
        infectionRateModel = new DefaultTableModel(1,1);
        infectionRateTable = new JTable(infectionRateModel);
        infectionRateTable.setShowGrid(true);
        box.add(infectionRateTable);
        estimateInfectionRate = new JCheckBox("estimate");
        box.add(estimateInfectionRate);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nInfectionRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nInfectionRateShiftsSpinner = new JSpinner(nInfectionRateShiftsModel);
        box.add(nInfectionRateShiftsSpinner);

        box.add(new JLabel("Change times:"));
        infectionRateChangeTimesModel = new DefaultTableModel(1,1);
        infectionRateChangeTimesTable = new JTable(infectionRateChangeTimesModel);
        infectionRateChangeTimesTable.setShowGrid(true);
        box.add(infectionRateChangeTimesTable);
        estimateInfectionRateShiftTimes = new JCheckBox("estimate");
        box.add(estimateInfectionRateShiftTimes);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Recovery rate:"));
        recoveryRateModel = new DefaultTableModel(1,1);
        recoveryRateTable = new JTable(recoveryRateModel);
        recoveryRateTable.setShowGrid(true);
        box.add(recoveryRateTable);
        estimateRecoveryRate = new JCheckBox("estimate");
        box.add(estimateRecoveryRate);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nRecoveryRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nRecoveryRateShiftsSpinner = new JSpinner(nRecoveryRateShiftsModel);
        box.add(nRecoveryRateShiftsSpinner);

        box.add(new JLabel("Change times:"));
        recoveryRateChangeTimesModel = new DefaultTableModel(1,1);
        recoveryRateChangeTimesTable = new JTable(recoveryRateChangeTimesModel);
        recoveryRateChangeTimesTable.setShowGrid(true);
        box.add(recoveryRateChangeTimesTable);
        estimateRecoveryRateShiftTimes = new JCheckBox("estimate");
        box.add(estimateRecoveryRateShiftTimes);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Psi sampling rate:"));
        psiSamplingRateModel = new DefaultTableModel(1,1);
        psiSamplingRateTable = new JTable(psiSamplingRateModel);
        psiSamplingRateTable.setShowGrid(true);
        box.add(psiSamplingRateTable);
        estimatePsiSamplingRate = new JCheckBox("estimate");
        box.add(estimatePsiSamplingRate);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nPsiSamplingRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nPsiSamplingRateShiftsSpinner = new JSpinner(nPsiSamplingRateShiftsModel);
        box.add(nPsiSamplingRateShiftsSpinner);

        box.add(new JLabel("Change times:"));
        psiSamplingRateChangeTimesModel = new DefaultTableModel(1,1);
        psiSamplingRateChangeTimesTable = new JTable(psiSamplingRateChangeTimesModel);
        psiSamplingRateChangeTimesTable.setShowGrid(true);
        box.add(psiSamplingRateChangeTimesTable);
        estimatePsiSamplingRateShiftTimes = new JCheckBox("estimate");
        box.add(estimatePsiSamplingRateShiftTimes);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Rho sampling probability:"));
        rhoSamplingProbTextField = new JTextField();
        box.add(rhoSamplingProbTextField);
        estimateRhoSamplingProb = new JCheckBox("estimate");
        box.add(estimateRhoSamplingProb);
        panel.add(box);

        add(panel);

        loadFromModel();

        // Event handlers

        emSelector.addActionListener(e -> saveToModel());
    }

    public void loadFromModel() {

        infectionRate = (RealParameter)epidemicModel.infectionRateInput.get();
        infectionRateModel.setColumnCount(infectionRate.getDimension());
        for (int i=0; i<infectionRate.getDimension(); i++) {
            infectionRateModel.setValueAt(infectionRate.getValue(i), 0, i);
        }

        nInfectionRateShiftsModel.setValue(infectionRate.getDimension()-1);
        infectionRateChangeTimes = epidemicModel.infectionRateShiftTimesInput.get();
        if (infectionRate.getDimension()>1) {
            infectionRateChangeTimesTable.setEnabled(true);
            infectionRateChangeTimesModel.setColumnCount(infectionRate.getDimension()-1);
            for (int i=0; i<infectionRate.getDimension()-1; i++) {
                infectionRateChangeTimesModel.setValueAt(
                        infectionRateChangeTimes.getValue(i), 0, i);
            }
        } else {
            infectionRateChangeTimesTable.setEnabled(false);
        }

        recoveryRate = (RealParameter)epidemicModel.recoveryRateInput.get();
        recoveryRateModel.setColumnCount(recoveryRate.getDimension());
        for (int i=0; i<recoveryRate.getDimension(); i++) {
            recoveryRateModel.setValueAt(recoveryRate.getValue(i), 0, i);
        }

        nRecoveryRateShiftsModel.setValue(recoveryRate.getDimension()-1);
        recoveryRateChangeTimes = epidemicModel.recoveryRateShiftTimesInput.get();
        if (recoveryRate.getDimension()>1) {
            recoveryRateChangeTimesTable.setEnabled(true);
            recoveryRateChangeTimesModel.setColumnCount(recoveryRate.getDimension()-1);
            for (int i=0; i<recoveryRate.getDimension()-1; i++) {
                recoveryRateChangeTimesModel.setValueAt(
                        recoveryRateChangeTimes.getValue(i), 0, i);
            }
        } else {
            recoveryRateChangeTimesTable.setEnabled(false);
        }

        psiSamplingRate = (RealParameter)epidemicModel.psiSamplingRateInput.get();
        psiSamplingRateModel.setColumnCount(psiSamplingRate.getDimension());
        for (int i=0; i<psiSamplingRate.getDimension(); i++) {
            psiSamplingRateModel.setValueAt(psiSamplingRate.getValue(i), 0, i);
        }

        nPsiSamplingRateShiftsModel.setValue(psiSamplingRate.getDimension()-1);
        psiSamplingRateChangeTimes = epidemicModel.psiSamplingRateShiftTimesInput.get();
        if (psiSamplingRate.getDimension()>1) {
            psiSamplingRateChangeTimesTable.setEnabled(true);
            psiSamplingRateChangeTimesModel.setColumnCount(psiSamplingRate.getDimension()-1);
            for (int i=0; i<psiSamplingRate.getDimension()-1; i++) {
                psiSamplingRateChangeTimesModel.setValueAt(
                        psiSamplingRateChangeTimes.getValue(i), 0, i);
            }
        } else {
            psiSamplingRateChangeTimesTable.setEnabled(false);
        }

        if (epidemicModel instanceof SISModel) {
            SISModel sisModel = (SISModel)epidemicModel;

            emSelectorModel.setSelectedItem("SIS");
            S0 = sisModel.S0Input.get();

        } else if (epidemicModel instanceof SIRModel){
            SIRModel sirModel = (SIRModel)epidemicModel;

            emSelectorModel.setSelectedItem("SIR");
            S0 = sirModel.S0Input.get();

        } else {
            BirthDeathModel bdModel = (BirthDeathModel)epidemicModel;

            emSelectorModel.setSelectedItem("Linear birth-death");
            S0 = null;
        }
    }

    public void saveToModel() {
        EpidemicModel newEpiModel;
        switch ((String)emSelectorModel.getSelectedItem()) {
            case "SIS":
                break;
            case "SIR":
                break;
            case "Linear birth-death":
                break;
        }

    }

}
