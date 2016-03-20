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
import javax.swing.border.BevelBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;

/**
 * Input editor for EpidemicModels.  Special message to anyone (including
 * my future self) who has to read this code: I am truly and deeply sorry.
 *
 * Notes:
 * 1. Input editor only produces models with times specified as ages relative
 * to most recent sample.  Default behaviour for EpidemicModels is for these
 * to be specified as proper times which increase into the future.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiModelInputEditor extends InputEditor.Base {

    EpidemicModel epidemicModel;
    String[] modelNames = {"Linear birth-death", "SIS", "SIR"};

    IntegerParameter S0;
    RealParameter infectionRate, recoveryRate, psiSamplingVariable,
            removalProb, rhoSamplingProb, rhoSamplingTime;
    RealParameter infectionRateChangeTimes, recoveryRateChangeTimes,
            psiSamplingVariableChangeTimes, removalProbChangeTimes;
    RealParameter epiOrigin;

    ComboBoxModel<String> emSelectorModel;
    DefaultTableModel infectionRateModel, infectionRateChangeTimesModel;
    JTable infectionRateTable, infectionRateChangeTimesTable;
    JLabel infectionRateChangeTimesLabel;

    DefaultTableModel recoveryRateModel, recoveryRateChangeTimesModel;
    JTable recoveryRateTable, recoveryRateChangeTimesTable;
    JLabel recoveryRateChangeTimesLabel;

    DefaultTableModel psiSamplingVariableModel, psiSamplingVariableChangeTimesModel;
    JTable psiSamplingVariableTable, psiSamplingVariableChangeTimesTable;
    JLabel psiSamplingVariableChangeTimesLabel;
    JCheckBox useSamplingProportionCheckBox;

    DefaultTableModel removalProbModel, removalProbChangeTimesModel;
    JTable removalProbTable, removalProbChangeTimesTable;
    JLabel removalProbChangeTimesLabel;

    SpinnerNumberModel nInfectionRateShiftsModel,
            nRecoveryRateShiftsModel,
            nPsiSamplingVariableShiftsModel, nRemovalProbShiftsModel;
    JSpinner nInfectionRateShiftsSpinner,
            nRecoveryRateShiftsSpinner,
            nPsiSamplingVariableShiftsSpinner, nRemovalProbShiftsSpinner;

    JTextField S0TextField, rhoSamplingProbTextField, originTextField;

    JCheckBox estimateS0, estimateInfectionRate, estimateRecoveryRate;
    JCheckBox estimateInfectionRateShiftTimes, estimateRecoveryRateShiftTimes;
    JCheckBox estimatePsiSamplingVariable, estimatePsiSamplingVariableShiftTimes;
    JCheckBox estimateRemovalProb, estimateRemovalProbShiftTimes;
    JCheckBox estimateRhoSamplingProb;
    JCheckBox estimateOrigin;

    EpiTrajPanel epiTrajPanel;

    boolean modelSaveInProgress = false;

    public EpiModelInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return EpidemicModel.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption bExpandOption, boolean bAddButtons) {

        // Set up fields
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;
        epidemicModel = (EpidemicModel) input.get();

        // Adds label to the left of input editor
        addInputLabel();

        // Create and lay out GUI components
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
        panel.setBorder(new EtchedBorder());

        Box box = Box.createHorizontalBox();
        box.add(new JLabel("Epidemic model to use:"));
        emSelectorModel = new DefaultComboBoxModel<>(modelNames);
        JComboBox<String> emSelector = new JComboBox<>(emSelectorModel);
        box.add(emSelector);
        box.add(Box.createHorizontalGlue());
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
        nInfectionRateShiftsModel = new SpinnerNumberModel(0, 0, 100, 1);
        nInfectionRateShiftsSpinner = new JSpinner(nInfectionRateShiftsModel);
        nInfectionRateShiftsSpinner.setMaximumSize(nInfectionRateShiftsSpinner.getPreferredSize());
        box.add(nInfectionRateShiftsSpinner);

        box.add(Box.createHorizontalGlue());

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
        nRecoveryRateShiftsModel = new SpinnerNumberModel(0, 0, 100, 1);
        nRecoveryRateShiftsSpinner = new JSpinner(nRecoveryRateShiftsModel);
        nRecoveryRateShiftsSpinner.setMaximumSize(nRecoveryRateShiftsSpinner.getPreferredSize());
        box.add(nRecoveryRateShiftsSpinner);

        box.add(Box.createHorizontalGlue());

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
        box.add(new JLabel("Psi sampling:"));
        psiSamplingVariableModel = new DefaultTableModel(1,1);
        psiSamplingVariableTable = new JTable(psiSamplingVariableModel);
        psiSamplingVariableTable.setShowGrid(true);
        psiSamplingVariableTable.setCellSelectionEnabled(false);
        box.add(psiSamplingVariableTable);
        estimatePsiSamplingVariable = new JCheckBox("estimate");
        box.add(estimatePsiSamplingVariable);
        psiSamplingPanel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nPsiSamplingVariableShiftsModel = new SpinnerNumberModel(0, 0, 100, 1);
        nPsiSamplingVariableShiftsSpinner = new JSpinner(nPsiSamplingVariableShiftsModel);
        nPsiSamplingVariableShiftsSpinner.setMaximumSize(nPsiSamplingVariableShiftsSpinner.getPreferredSize());
        box.add(nPsiSamplingVariableShiftsSpinner);

        box.add(Box.createHorizontalGlue());

        psiSamplingVariableChangeTimesLabel = new JLabel("Change times:");
        box.add(psiSamplingVariableChangeTimesLabel);
        psiSamplingVariableChangeTimesModel = new DefaultTableModel(1,1);
        psiSamplingVariableChangeTimesTable = new JTable(psiSamplingVariableChangeTimesModel);
        psiSamplingVariableChangeTimesTable.setShowGrid(true);
        psiSamplingVariableChangeTimesTable.setCellSelectionEnabled(false);
        box.add(psiSamplingVariableChangeTimesTable);
        estimatePsiSamplingVariableShiftTimes = new JCheckBox("estimate");
        box.add(estimatePsiSamplingVariableShiftTimes);
        psiSamplingPanel.add(box);

        box = Box.createHorizontalBox();
        useSamplingProportionCheckBox = new JCheckBox(
                "Use sampling proportion instead of sampling rate");
        box.add(useSamplingProportionCheckBox);
        box.add(Box.createHorizontalGlue());
        psiSamplingPanel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Removal proportion:"));
        removalProbModel = new DefaultTableModel(1,1);
        removalProbTable = new JTable(removalProbModel);
        removalProbTable.setShowGrid(true);
        removalProbTable.setCellSelectionEnabled(false);
        box.add(removalProbTable);
        estimateRemovalProb = new JCheckBox("estimate");
        box.add(estimateRemovalProb);
        psiSamplingPanel.add(box);

        box = Box.createHorizontalBox();
        box.add(new JLabel("Num. changes:"));
        nRemovalProbShiftsModel = new SpinnerNumberModel(0, 0, 100, 1);
        nRemovalProbShiftsSpinner = new JSpinner(nRemovalProbShiftsModel);
        nRemovalProbShiftsSpinner.setMaximumSize(nRemovalProbShiftsSpinner.getPreferredSize());
        box.add(nRemovalProbShiftsSpinner);

        box.add(Box.createHorizontalGlue());

        removalProbChangeTimesLabel = new JLabel("Change times:");
        box.add(removalProbChangeTimesLabel);
        removalProbChangeTimesModel = new DefaultTableModel(1,1);
        removalProbChangeTimesTable = new JTable(removalProbChangeTimesModel);
        removalProbChangeTimesTable.setShowGrid(true);
        removalProbChangeTimesTable.setCellSelectionEnabled(false);
        box.add(removalProbChangeTimesTable);
        estimateRemovalProbShiftTimes = new JCheckBox("estimate");
        box.add(estimateRemovalProbShiftTimes);

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

        box = Box.createHorizontalBox();
        box.add(new JLabel("Time before present of epidemic origin:"));
        originTextField = new JTextField();
        box.add(originTextField);
        estimateOrigin = new JCheckBox("estimate");
        box.add(estimateOrigin);
        panel.add(box);

        // Prevalence trajectory plot
//        box = Box.createHorizontalBox();
        epiTrajPanel = new EpiTrajPanel(epidemicModel, 5);
        epiTrajPanel.setBorder(new BevelBorder(BevelBorder.LOWERED));
        panel.add(epiTrajPanel);

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

        psiSamplingVariableModel.addTableModelListener(tableListener);
        estimatePsiSamplingVariable.addItemListener(e -> saveToModel());
        psiSamplingVariableChangeTimesModel.addTableModelListener(tableListener);
        estimatePsiSamplingVariableShiftTimes.addItemListener(e -> saveToModel());
        nPsiSamplingVariableShiftsModel.addChangeListener(e -> saveToModel());
        useSamplingProportionCheckBox.addItemListener(e -> saveToModel());

        removalProbModel.addTableModelListener(tableListener);
        estimateRemovalProb.addItemListener(e -> saveToModel());
        removalProbChangeTimesModel.addTableModelListener(tableListener);
        estimateRemovalProbShiftTimes.addItemListener(e -> saveToModel());
        nRemovalProbShiftsModel.addChangeListener(e -> saveToModel());

        S0TextField.addActionListener(e -> saveToModel());
        estimateS0.addItemListener(e -> saveToModel());
        rhoSamplingProbTextField.addActionListener(e -> saveToModel());
        estimateRhoSamplingProb.addItemListener(e -> saveToModel());

        originTextField.addActionListener(e -> saveToModel());
        estimateOrigin.addItemListener(e -> saveToModel());
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

        psiSamplingVariable = (RealParameter)epidemicModel.psiSamplingVariableInput.get();
        psiSamplingVariableChangeTimes = (RealParameter)epidemicModel.psiSamplingVariableShiftTimesInput.get();
        loadModelRateParameters(psiSamplingVariable, psiSamplingVariableModel,
                estimatePsiSamplingVariable, nPsiSamplingVariableShiftsModel,
                psiSamplingVariableChangeTimes, psiSamplingVariableChangeTimesLabel,
                psiSamplingVariableChangeTimesModel, estimatePsiSamplingVariableShiftTimes);
        useSamplingProportionCheckBox.setSelected(epidemicModel.usePsiSamplingProportionInput.get());

        removalProb = (RealParameter)epidemicModel.removalProbInput.get();
        removalProbChangeTimes = (RealParameter)epidemicModel.removalProbShiftTimesInput.get();
        loadModelRateParameters(removalProb, removalProbModel,
                estimateRemovalProb, nRemovalProbShiftsModel,
                removalProbChangeTimes, removalProbChangeTimesLabel,
                removalProbChangeTimesModel, estimateRemovalProbShiftTimes);


        rhoSamplingProb = (RealParameter)epidemicModel.rhoSamplingProbInput.get();
        rhoSamplingProbTextField.setText(String.valueOf(rhoSamplingProb.getValue()));
        rhoSamplingTime = (RealParameter) epidemicModel.rhoSamplingTimeInput.get();
        if (rhoSamplingProb.getValue()>0) {
            estimateRhoSamplingProb.setEnabled(true);
            estimateRhoSamplingProb.setSelected(rhoSamplingProb.isEstimatedInput.get());
        } else {
            estimateRhoSamplingProb.setEnabled(false);
            estimateRhoSamplingProb.setSelected(false);
        }

        epiOrigin = (RealParameter)epidemicModel.originInput.get();
        originTextField.setText(String.valueOf(epiOrigin.getArrayValue()));
        estimateOrigin.setSelected(epiOrigin.isEstimatedInput.get());

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

        epiTrajPanel.updateChart();
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

        PartitionContext partitionContext = doc.getContextFor(m_beastObject);
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

            saveModelRateParameters("psiSamplingVariable",
                    psiSamplingVariable, psiSamplingVariableModel, estimatePsiSamplingVariable,
                    nPsiSamplingVariableShiftsModel,
                    psiSamplingVariableChangeTimes, psiSamplingVariableChangeTimesModel, estimatePsiSamplingVariableShiftTimes,
                    epidemicModel, partitionID);
            epidemicModel.usePsiSamplingProportionInput.setValue(useSamplingProportionCheckBox.isSelected(), epidemicModel);

            saveModelRateParameters("removalProb",
                    removalProb, removalProbModel, estimateRemovalProb,
                    nRemovalProbShiftsModel,
                    removalProbChangeTimes, removalProbChangeTimesModel, estimateRemovalProbShiftTimes,
                    epidemicModel, partitionID);

            rhoSamplingProb.valuesInput.setValue(rhoSamplingProbTextField.getText(), rhoSamplingProb);
            rhoSamplingProb.isEstimatedInput.setValue(estimateRhoSamplingProb.isSelected(), rhoSamplingProb);

            epidemicModel.setInputValue("rhoSamplingProb", rhoSamplingProb);
            rhoSamplingProb.initAndValidate();

            epidemicModel.setInputValue("rhoSamplingTime", rhoSamplingTime);
            epidemicModel.setInputValue("rhoSamplingTimesBackward", true);

            epiOrigin.valuesInput.setValue(originTextField.getText(), epiOrigin);
            epiOrigin.isEstimatedInput.setValue(estimateOrigin.isSelected(), epiOrigin);
            epidemicModel.setInputValue("origin", epiOrigin);
            epiOrigin.initAndValidate();

            epidemicModel.initAndValidate();

            m_input.setValue(epidemicModel, m_beastObject);
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
            epidemicModel.setInputValue(paramName + "ShiftTimesBackward", true);
        }
    }
}
