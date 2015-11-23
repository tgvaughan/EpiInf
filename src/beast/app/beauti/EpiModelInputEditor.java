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
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiModelInputEditor extends InputEditor.Base {

    EpidemicModel epidemicModel;
    String[] modelNames = {"Linear birth-death", "SIS", "SIR"};

    IntegerParameter S0;
    RealParameter infectionRate, recoveryRate, psiSamplingRate, rhoSamplingProb, rhoSamplingTime;
    RealParameter infectionRateChangeTimes, recoveryRateChangeTimes, psiSamplingRateChangeTimes;

    ComboBoxModel<String> emSelectorModel;
    DefaultTableModel infectionRateModel, infectionRateChangeTimesModel;
    JTable infectionRateTable, infectionRateChangeTimesTable;
    JLabel infectionRateChangeTimesLabel;

    DefaultTableModel recoveryRateModel, recoveryRateChangeTimesModel;
    JTable recoveryRateTable, recoveryRateChangeTimesTable;
    JLabel recoveryRateChangeTimesLabel;

    DefaultTableModel psiSamplingRateModel, psiSamplingRateChangeTimesModel;
    JTable psiSamplingRateTable, psiSamplingRateChangeTimesTable;
    JLabel psiSamplingRateChangeTimesLabel;

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
        infectionRateTable.setCellSelectionEnabled(false);
        box.add(infectionRateTable);
        estimateInfectionRate = new JCheckBox("estimate");
        box.add(estimateInfectionRate);
        panel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nInfectionRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nInfectionRateShiftsSpinner = new JSpinner(nInfectionRateShiftsModel);
        box.add(nInfectionRateShiftsSpinner);

        infectionRateChangeTimesLabel = new JLabel("Change times:");
        box.add(infectionRateChangeTimesLabel);
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

        recoveryRateChangeTimesLabel = new JLabel("Change times:");
        box.add(recoveryRateChangeTimesLabel);
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

        psiSamplingRateChangeTimesLabel = new JLabel("Change times:");
        box.add(psiSamplingRateChangeTimesLabel);
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

        box = Box.createHorizontalBox();
        box.add(new JLabel("Initial Susceptible population size:"));
        S0TextField = new JTextField();
        box.add(S0TextField);
        estimateS0 = new JCheckBox("estimate");
        box.add(estimateS0);
        panel.add(box);

        add(panel);

        loadFromModel();

        // Event handlers

        emSelector.addItemListener(e -> saveToModel());

        TableModelListener tableListener = e -> {
            if (e.getType() == TableModelEvent.UPDATE)
                saveToModel();
        };

        infectionRateModel.addTableModelListener(tableListener);
        estimateInfectionRate.addItemListener(e -> saveToModel());
        infectionRateChangeTimesModel.addTableModelListener(tableListener);
        estimateInfectionRateShiftTimes.addItemListener(e -> saveToModel());
        nInfectionRateShiftsModel.addChangeListener(e -> saveToModel());

        recoveryRateModel.addTableModelListener(tableListener);
        estimateRecoveryRate.addItemListener(e -> saveToModel());
        recoveryRateChangeTimesModel.addTableModelListener(tableListener);
        estimateRecoveryRateShiftTimes.addItemListener(e -> saveToModel());
        nRecoveryRateShiftsModel.addChangeListener(e -> saveToModel());

        psiSamplingRateModel.addTableModelListener(tableListener);
        estimatePsiSamplingRate.addItemListener(e -> saveToModel());
        psiSamplingRateChangeTimesModel.addTableModelListener(tableListener);
        estimatePsiSamplingRateShiftTimes.addItemListener(e -> saveToModel());
        nPsiSamplingRateShiftsModel.addChangeListener(e -> saveToModel());

        S0TextField.addActionListener(e -> saveToModel());
        estimateS0.addItemListener(e -> saveToModel());
        rhoSamplingProbTextField.addActionListener(e -> saveToModel());
        estimateRhoSamplingProb.addItemListener(e -> saveToModel());

    }

    public void loadFromModel() {

        infectionRate = (RealParameter)epidemicModel.infectionRateInput.get();
        infectionRateModel.setColumnCount(infectionRate.getDimension());
        for (int i=0; i<infectionRate.getDimension(); i++) {
            infectionRateModel.setValueAt(infectionRate.getValue(i), 0, i);
        }
        estimateInfectionRate.setSelected(infectionRate.isEstimatedInput.get());

        nInfectionRateShiftsModel.setValue(infectionRate.getDimension()-1);
        infectionRateChangeTimes = epidemicModel.infectionRateShiftTimesInput.get();
        if (infectionRate.getDimension()>1) {
            infectionRateChangeTimesTable.setEnabled(true);
            infectionRateChangeTimesLabel.setEnabled(true);
            infectionRateChangeTimesModel.setColumnCount(infectionRate.getDimension()-1);
            for (int i=0; i<infectionRate.getDimension()-1; i++) {
                infectionRateChangeTimesModel.setValueAt(
                        infectionRateChangeTimes.getValue(i), 0, i);
            }
            estimateInfectionRateShiftTimes.setEnabled(true);
            estimateInfectionRateShiftTimes.setSelected(infectionRateChangeTimes.isEstimatedInput.get());
        } else {
            infectionRateChangeTimesTable.setEnabled(false);
            infectionRateChangeTimesLabel.setEnabled(false);
            estimateInfectionRateShiftTimes.setSelected(false);
            estimateInfectionRateShiftTimes.setEnabled(false);
        }

        recoveryRate = (RealParameter)epidemicModel.recoveryRateInput.get();
        recoveryRateModel.setColumnCount(recoveryRate.getDimension());
        for (int i=0; i<recoveryRate.getDimension(); i++) {
            recoveryRateModel.setValueAt(recoveryRate.getValue(i), 0, i);
        }
        estimateRecoveryRate.setSelected(recoveryRate.isEstimatedInput.get());

        nRecoveryRateShiftsModel.setValue(recoveryRate.getDimension()-1);
        recoveryRateChangeTimes = epidemicModel.recoveryRateShiftTimesInput.get();
        if (recoveryRate.getDimension()>1) {
            recoveryRateChangeTimesTable.setEnabled(true);
            recoveryRateChangeTimesLabel.setEnabled(true);
            recoveryRateChangeTimesModel.setColumnCount(recoveryRate.getDimension()-1);
            for (int i=0; i<recoveryRate.getDimension()-1; i++) {
                recoveryRateChangeTimesModel.setValueAt(
                        recoveryRateChangeTimes.getValue(i), 0, i);
            }
            estimateRecoveryRateShiftTimes.setEnabled(true);
            estimateRecoveryRateShiftTimes.setSelected(recoveryRateChangeTimes.isEstimatedInput.get());
        } else {
            recoveryRateChangeTimesTable.setEnabled(false);
            recoveryRateChangeTimesLabel.setEnabled(false);
            estimateRecoveryRateShiftTimes.setEnabled(false);
            estimateRecoveryRateShiftTimes.setSelected(false);
        }

        psiSamplingRate = (RealParameter)epidemicModel.psiSamplingRateInput.get();
        psiSamplingRateModel.setColumnCount(psiSamplingRate.getDimension());
        for (int i=0; i<psiSamplingRate.getDimension(); i++) {
            psiSamplingRateModel.setValueAt(psiSamplingRate.getValue(i), 0, i);
        }
        estimatePsiSamplingRate.setSelected(psiSamplingRate.isEstimatedInput.get());

        nPsiSamplingRateShiftsModel.setValue(psiSamplingRate.getDimension()-1);
        psiSamplingRateChangeTimes = epidemicModel.psiSamplingRateShiftTimesInput.get();
        if (psiSamplingRate.getDimension()>1) {
            psiSamplingRateChangeTimesTable.setEnabled(true);
            psiSamplingRateChangeTimesLabel.setEnabled(true);
            psiSamplingRateChangeTimesModel.setColumnCount(psiSamplingRate.getDimension()-1);
            for (int i=0; i<psiSamplingRate.getDimension()-1; i++) {
                psiSamplingRateChangeTimesModel.setValueAt(
                        psiSamplingRateChangeTimes.getValue(i), 0, i);
            }
            estimatePsiSamplingRateShiftTimes.setEnabled(true);
            estimatePsiSamplingRateShiftTimes.setSelected(psiSamplingRateChangeTimes.isEstimatedInput.get());
        } else {
            psiSamplingRateChangeTimesTable.setEnabled(false);
            psiSamplingRateChangeTimesLabel.setEnabled(false);
            estimatePsiSamplingRateShiftTimes.setEnabled(false);
            estimatePsiSamplingRateShiftTimes.setSelected(false);
        }

        rhoSamplingProb = (RealParameter)epidemicModel.rhoSamplingProbInput.get();
        rhoSamplingProbTextField.setText(String.valueOf(rhoSamplingProb.getValue()));
        rhoSamplingTime = epidemicModel.rhoSamplingTimeInput.get();
        if (rhoSamplingProb.getValue()>0) {
            estimateRhoSamplingProb.setEnabled(true);
            estimateRhoSamplingProb.setSelected(rhoSamplingProb.isEstimatedInput.get());
        } else {
            estimateRhoSamplingProb.setEnabled(false);
            estimateRhoSamplingProb.setSelected(false);
        }

        if (epidemicModel instanceof SISModel) {
            SISModel sisModel = (SISModel)epidemicModel;

            emSelectorModel.setSelectedItem("SIS");
            S0 = sisModel.S0Input.get();
            S0TextField.setText(String.valueOf(S0.getValue()));
            S0TextField.setEnabled(true);
            estimateS0.setEnabled(true);
            estimateS0.setSelected(S0.isEstimatedInput.get());

        } else if (epidemicModel instanceof SIRModel){
            SIRModel sirModel = (SIRModel)epidemicModel;

            emSelectorModel.setSelectedItem("SIR");
            S0 = sirModel.S0Input.get();
            S0TextField.setText(String.valueOf(S0.getValue()));
            S0TextField.setEnabled(true);
            estimateS0.setEnabled(true);
            estimateS0.setSelected(S0.isEstimatedInput.get());

        } else {
            emSelectorModel.setSelectedItem("Linear birth-death");
            S0TextField.setText("");
            S0TextField.setEnabled(false);
            estimateS0.setEnabled(false);
            estimateS0.setSelected(false);
        }
    }

    public void saveToModel() {


        infectionRate.setDimension(infectionRateModel.getColumnCount());
        StringBuilder sbInfectionRate = new StringBuilder();
        for (int i=0; i<infectionRateTable.getColumnCount(); i++)
            sbInfectionRate.append(" ").append(infectionRateModel.getValueAt(0, i));
        infectionRate.valuesInput.setValue(sbInfectionRate.toString(), infectionRate);
        infectionRate.isEstimatedInput.setValue(estimateInfectionRate.isSelected(), infectionRate);

        recoveryRate.setDimension(recoveryRateModel.getColumnCount());
        StringBuilder sbRecoveryRate = new StringBuilder();
        for (int i=0; i<recoveryRateTable.getColumnCount(); i++)
            sbRecoveryRate.append(" ").append(recoveryRateModel.getValueAt(0, i));
        recoveryRate.valuesInput.setValue(sbRecoveryRate.toString(), recoveryRate);
        recoveryRate.isEstimatedInput.setValue(estimateRecoveryRate.isSelected(), recoveryRate);

        psiSamplingRate.setDimension(psiSamplingRateModel.getColumnCount());
        StringBuilder sbPsiSamplingRate = new StringBuilder();
        for (int i=0; i<psiSamplingRateTable.getColumnCount(); i++)
            sbPsiSamplingRate.append(" ").append(psiSamplingRateModel.getValueAt(0, i));
        psiSamplingRate.valuesInput.setValue(sbPsiSamplingRate.toString(), psiSamplingRate);
        psiSamplingRate.isEstimatedInput.setValue(estimatePsiSamplingRate.isSelected(), psiSamplingRate);

        rhoSamplingProb.valuesInput.setValue(rhoSamplingProbTextField.getText(), rhoSamplingProb);
        rhoSamplingProb.isEstimatedInput.setValue(estimateRhoSamplingProb.isSelected(), rhoSamplingProb);

        try {
            switch ((String)emSelectorModel.getSelectedItem()) {
                case "SIS":
                    epidemicModel = new SISModel();
                    break;
                case "SIR":
                    epidemicModel = new SIRModel();
                    break;
                case "Linear birth-death":
                    epidemicModel = new BirthDeathModel();
                    break;
                default:
                    throw new RuntimeException("Incompatible epidemic model.");
            }

            if (epidemicModel instanceof SISModel || epidemicModel instanceof SIRModel) {
                if (S0 == null) {
                    try {
                        PartitionContext pc = doc.getContextFor(m_plugin);
                        String S0id = "S0.t:" + pc.tree;
                        if (doc.pluginmap.containsKey(S0id))
                            S0 = (IntegerParameter)doc.pluginmap.get(S0id);
                        else {
                            S0 = new IntegerParameter("199");
                            S0.setID("S0.t:" + pc.tree);
                            doc.registerPlugin(S0);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
                if (!S0TextField.getText().trim().isEmpty())
                    S0.valuesInput.setValue(S0TextField.getText(), S0);

                S0.isEstimatedInput.setValue(estimateS0.isSelected(), S0);

                epidemicModel.setInputValue("S0", S0);
                S0.initAndValidate();
            } else {
//                if (S0 != null)
//                    doc.unregisterPlugin(S0);
            }

            epidemicModel.setInputValue("infectionRate", infectionRate);
            epidemicModel.setInputValue("infectionRateShiftTimes", infectionRateChangeTimes);
            epidemicModel.setInputValue("recoveryRate", recoveryRate);
            epidemicModel.setInputValue("recoveryRateShiftTimes", recoveryRateChangeTimes);
            epidemicModel.setInputValue("psiSamplingRate", psiSamplingRate);
            epidemicModel.setInputValue("psiSamplingRateShiftTimes", psiSamplingRateChangeTimes);
            epidemicModel.setInputValue("rhoSamplingProb", rhoSamplingProb);
            epidemicModel.setInputValue("rhoSamplingTime", rhoSamplingTime);

            m_input.setValue(epidemicModel, m_plugin);

            infectionRate.initAndValidate();
            recoveryRate.initAndValidate();
            psiSamplingRate.initAndValidate();
            rhoSamplingProb.initAndValidate();
            epidemicModel.initAndValidate();
        } catch (Exception e) {
            System.err.println("Error updating epidemic model.");
            e.printStackTrace();
        }
        refreshPanel();
    }
}
