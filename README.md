EpiInf
======

Small package for performing inference of compartmental
epidemiological model parameters.  Uses the *exact* likelihood for the
parameters given a tree and a stochastic epidemic trajectory.

This package is under development.

Development Roadmap
-------------------

- [x] State classes
    - [x] EpidemicEvent
    - [x] EpidemicState
    - [x] TreeEvent
    - [x] TreeEventList
    - [x] TreeEventList

- [x] Model classes
    - [x] EpidemicModel
    - [x] SIRModel
    - [x] SISModel
    - [x] BirthDeathModel

- [x] Simulation classes
    - [x] TransmissionTreeSimulator
    - [x] SimulatedTransmissionTree
    - [x] EpidemicTrajectory
    - [x] TrajectorySimulator
    - [x] SimulatedTrajectory

- [x] PF/SMC likelihood class
    - [x] Draft implementation
    - [x] Passes comparison with ExpoTree
