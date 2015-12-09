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

import beast.app.draw.BEASTObjectPanel;
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
import javax.swing.border.TitledBorder;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

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

    boolean modelSaveInProgress = false;

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

        JPanel infectionPanel = new JPanel();
        infectionPanel.setLayout(new BoxLayout(infectionPanel, BoxLayout.PAGE_AXIS));
        infectionPanel.setBorder(new TitledBorder("Infection"));

        box = Box.createHorizontalBox();
        box.add(new JLabel("Infection rate:"));
        infectionRateModel = new DefaultTableModel(1,1);
        infectionRateTable = new JTable(infectionRateModel);
        infectionRateTable.setShowGrid(true);
        infectionRateTable.setCellSelectionEnabled(false);
        box.add(infectionRateTable);
        estimateInfectionRate = new JCheckBox("estimate");
        box.add(estimateInfectionRate);
        infectionPanel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nInfectionRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nInfectionRateShiftsSpinner = new JSpinner(nInfectionRateShiftsModel);
        box.add(nInfectionRateShiftsSpinner);

        infectionRateChangeTimesLabel = new JLabel("Change times:");
        box.add(infectionRateChangeTimesLabel);
        infectionRateChangeTimesModel = new DefaultTableModel(1,0);
        infectionRateChangeTimesTable = new JTable(infectionRateChangeTimesModel);
        infectionRateChangeTimesTable.setShowGrid(true);
        infectionRateChangeTimesTable.setCellSelectionEnabled(false);
        box.add(infectionRateChangeTimesTable);
        estimateInfectionRateShiftTimes = new JCheckBox("estimate");
        box.add(estimateInfectionRateShiftTimes);
        infectionPanel.add(box);

        panel.add(infectionPanel);

        JPanel recoveryPanel = new JPanel();
        recoveryPanel.setLayout(new BoxLayout(recoveryPanel, BoxLayout.PAGE_AXIS));
        recoveryPanel.setBorder(new TitledBorder("Recovery"));

        box = Box.createHorizontalBox();
        box.add(new JLabel("Recovery rate:"));
        recoveryRateModel = new DefaultTableModel(1,1);
        recoveryRateTable = new JTable(recoveryRateModel);
        recoveryRateTable.setShowGrid(true);
        recoveryRateTable.setCellSelectionEnabled(false);
        box.add(recoveryRateTable);
        estimateRecoveryRate = new JCheckBox("estimate");
        box.add(estimateRecoveryRate);
        recoveryPanel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nRecoveryRateShiftsModel = new SpinnerNumberModel(0, 0, Integer.MAX_VALUE, 1);
        nRecoveryRateShiftsSpinner = new JSpinner(nRecoveryRateShiftsModel);
        box.add(nRecoveryRateShiftsSpinner);

        recoveryRateChangeTimesLabel = new JLabel("Change times:");
        box.add(recoveryRateChangeTimesLabel);
        recoveryRateChangeTimesModel = new DefaultTableModel(1,0);
        recoveryRateChangeTimesTable = new JTable(recoveryRateChangeTimesModel);
        recoveryRateChangeTimesTable.setShowGrid(true);
        recoveryRateChangeTimesTable.setCellSelectionEnabled(false);
        box.add(recoveryRateChangeTimesTable);
        estimateRecoveryRateShiftTimes = new JCheckBox("estimate");
        box.add(estimateRecoveryRateShiftTimes);
        recoveryPanel.add(box);

        panel.add(recoveryPanel);

        JPanel psiSamplingPanel = new JPanel();
        psiSamplingPanel.setLayout(new BoxLayout(psiSamplingPanel, BoxLayout.PAGE_AXIS));
        psiSamplingPanel.setBorder(new TitledBorder("Psi Sampling"));

        box = Box.createHorizontalBox();
        box.add(new JLabel("Psi sampling rate:"));
        psiSamplingRateModel = new DefaultTableModel(1,1);
        psiSamplingRateTable = new JTable(psiSamplingRateModel);
        psiSamplingRateTable.setShowGrid(true);
        psiSamplingRateTable.setCellSelectionEnabled(false);
        box.add(psiSamplingRateTable);
        estimatePsiSamplingRate = new JCheckBox("estimate");
        box.add(estimatePsiSamplingRate);
        psiSamplingPanel.add(box);

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
        psiSamplingRateChangeTimesTable.setCellSelectionEnabled(false);
        box.add(psiSamplingRateChangeTimesTable);
        estimatePsiSamplingRateShiftTimes = new JCheckBox("estimate");
        box.add(estimatePsiSamplingRateShiftTimes);
        psiSamplingPanel.add(box);

        panel.add(psiSamplingPanel);

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
        infectionRateChangeTimes = (RealParameter) epidemicModel.infectionRateShiftTimesInput.get();
        loadModelRateParameters(infectionRate, infectionRateModel,
                estimateInfectionRate, nInfectionRateShiftsModel,
                infectionRateChangeTimes, infectionRateChangeTimesLabel,
                infectionRateChangeTimesModel, estimateInfectionRateShiftTimes);

        recoveryRate = (RealParameter)epidemicModel.recoveryRateInput.get();
        recoveryRateChangeTimes = (RealParameter) epidemicModel.recoveryRateShiftTimesInput.get();
        loadModelRateParameters(recoveryRate, recoveryRateModel,
                estimateRecoveryRate, nRecoveryRateShiftsModel,
                recoveryRateChangeTimes, recoveryRateChangeTimesLabel,
                recoveryRateChangeTimesModel, estimateRecoveryRateShiftTimes);

        psiSamplingRate = (RealParameter)epidemicModel.psiSamplingRateInput.get();
        psiSamplingRateChangeTimes = (RealParameter)epidemicModel.psiSamplingRateShiftTimesInput.get();
        loadModelRateParameters(psiSamplingRate, psiSamplingRateModel,
                estimatePsiSamplingRate, nPsiSamplingRateShiftsModel,
                psiSamplingRateChangeTimes, psiSamplingRateChangeTimesLabel,
                psiSamplingRateChangeTimesModel, estimatePsiSamplingRateShiftTimes);

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

    protected void loadModelRateParameters(RealParameter rate,
                                           DefaultTableModel rateModel,
                                           JCheckBox estimateRate,
                                           SpinnerModel nRateShiftsModel,
                                           RealParameter rateChangeTimes,
                                           JLabel rateChangeTimesLabel,
                                           DefaultTableModel rateChangeTimesModel,
                                           JCheckBox estimateRateShiftTimes) {

        rateModel.setColumnCount(rate.getDimension());
        for (int i=0; i<rate.getDimension(); i++) {
            rateModel.setValueAt(rate.getValue(i), 0, i);
        }
        estimateRate.setSelected(rate.isEstimatedInput.get());

        nRateShiftsModel.setValue(rate.getDimension()-1);
        if (rateChangeTimes != null) {
            rateChangeTimesLabel.setVisible(true);
            rateChangeTimesModel.setColumnCount(rate.getDimension()-1);
            for (int i=0; i<rate.getDimension()-1; i++) {
                rateChangeTimesModel.setValueAt(
                        rateChangeTimes.getValue(i), 0, i);
            }
            estimateRateShiftTimes.setEnabled(true);
            estimateRateShiftTimes.setSelected(rateChangeTimes.isEstimatedInput.get());
        } else {
            rateChangeTimesLabel.setVisible(false);
            rateChangeTimesModel.setColumnCount(0);
            estimateRateShiftTimes.setSelected(false);
            estimateRateShiftTimes.setEnabled(false);
        }
    }

    public void saveToModel() {

        if (modelSaveInProgress)
            return;

        modelSaveInProgress = true;

        PartitionContext partitionContext = doc.getContextFor(m_plugin);
        String partitionID = ".t:" + partitionContext.tree;

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
                        String S0id = "S0" + partitionID;
                        if (doc.pluginmap.containsKey(S0id))
                            S0 = (IntegerParameter)doc.pluginmap.get(S0id);
                        else {
                            S0 = new IntegerParameter("199");
                            S0.setID(S0id);
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
                if (S0 != null) {
                    S0.isEstimatedInput.setValue(false, S0);
                }
            }

            saveModelRateParameters("infectionRate",
                    infectionRate, infectionRateModel, estimateInfectionRate,
                    nInfectionRateShiftsModel,
                    infectionRateChangeTimes, infectionRateChangeTimesModel, estimateInfectionRateShiftTimes,
                    epidemicModel, partitionID);

            saveModelRateParameters("recoveryRate",
                    recoveryRate, recoveryRateModel, estimateRecoveryRate,
                    nRecoveryRateShiftsModel,
                    recoveryRateChangeTimes, recoveryRateChangeTimesModel, estimateRecoveryRateShiftTimes,
                    epidemicModel, partitionID);

            saveModelRateParameters("psiSamplingRate",
                    psiSamplingRate, psiSamplingRateModel, estimatePsiSamplingRate,
                    nPsiSamplingRateShiftsModel,
                    psiSamplingRateChangeTimes, psiSamplingRateChangeTimesModel, estimatePsiSamplingRateShiftTimes,
                    epidemicModel, partitionID);

            rhoSamplingProb.valuesInput.setValue(rhoSamplingProbTextField.getText(), rhoSamplingProb);
            rhoSamplingProb.isEstimatedInput.setValue(estimateRhoSamplingProb.isSelected(), rhoSamplingProb);

            epidemicModel.setInputValue("rhoSamplingProb", rhoSamplingProb);
            rhoSamplingProb.initAndValidate();

            epidemicModel.setInputValue("rhoSamplingTime", rhoSamplingTime);

            epidemicModel.initAndValidate();

            m_input.setValue(epidemicModel, m_plugin);
            BEASTObjectPanel.addPluginToMap(epidemicModel, doc);

        } catch (Exception e) {
            System.err.println("Error updating epidemic model.");
            e.printStackTrace();
        }

        modelSaveInProgress = false;

        refreshPanel();
    }

    /**
     * Frightening method responsible for storing contents of GUI elements to
     * the BEAST model held by BEAUti.
     *
     * @param paramName
     * @param rate
     * @param rateModel
     * @param estimateRate
     * @param nRateShiftsModel
     * @param rateChangeTimes
     * @param rateChangeTimesModel
     * @param estimateRateShiftTimes
     * @param epidemicModel
     * @param partionID
     * @throws Exception
     */
    protected void saveModelRateParameters(String paramName, RealParameter rate,
                                           DefaultTableModel rateModel,
                                           JCheckBox estimateRate,
                                           SpinnerModel nRateShiftsModel,
                                           RealParameter rateChangeTimes,
                                           DefaultTableModel rateChangeTimesModel,
                                           JCheckBox estimateRateShiftTimes,
                                           EpidemicModel epidemicModel,
                                           String partionID) throws Exception {

        rateModel.setColumnCount((int)nRateShiftsModel.getValue()+1);
        for (int i=0; i<rateModel.getColumnCount(); i++) {
            if (rateModel.getValueAt(0, i) == null) {
                // Only ever true for i>0
                rateModel.setValueAt(rateModel.getValueAt(0, i-1), 0, i);
            }
        }
        rate.setDimension(rateModel.getColumnCount());
        StringBuilder sbInfectionRate = new StringBuilder();
        for (int i=0; i<rateModel.getColumnCount(); i++) {
            sbInfectionRate.append(" ").append(rateModel.getValueAt(0, i));
        }
        rate.valuesInput.setValue(sbInfectionRate.toString(), rate);
        rate.isEstimatedInput.setValue(estimateRate.isSelected(), rate);

        rateChangeTimesModel.setColumnCount((int)nRateShiftsModel.getValue());
        for (int i=0; i<rateChangeTimesModel.getColumnCount(); i++) {
            if (rateChangeTimesModel.getValueAt(0, i) == null)
                rateChangeTimesModel.setValueAt(0.0, 0, i);
        }
        if (rateChangeTimesModel.getColumnCount()>0) {
            if (rateChangeTimes == null) {
                String irctID = paramName + "ChangeTimes" + partionID;
                if (doc.pluginmap.containsKey(irctID))
                    rateChangeTimes = (RealParameter)doc.pluginmap.get(irctID);
                else {
                    rateChangeTimes = new RealParameter("0.0");
                    rateChangeTimes.setID(irctID);
                }
            }
            rateChangeTimes.setDimension(rateChangeTimesModel.getColumnCount());

            StringBuilder sbInfectionRateChangeTimes = new StringBuilder();
            for (int i=0; i<rateChangeTimesModel.getColumnCount(); i++) {
                sbInfectionRateChangeTimes.append(" ").append(rateChangeTimesModel.getValueAt(0, i));
            }
            rateChangeTimes.valuesInput.setValue(sbInfectionRateChangeTimes.toString(), rateChangeTimes);
            rateChangeTimes.isEstimatedInput.setValue(estimateRateShiftTimes.isSelected(), rateChangeTimes);
        } else {
            if (rateChangeTimes != null)
                rateChangeTimes.isEstimatedInput.setValue(false, rateChangeTimes);
        }

        epidemicModel.setInputValue(paramName, rate);
        rate.initAndValidate();
        if (rateChangeTimesModel.getColumnCount()>0) {
            epidemicModel.setInputValue(paramName + "ShiftTimes", rateChangeTimes);
            rateChangeTimes.initAndValidate();
        }
    }
}
