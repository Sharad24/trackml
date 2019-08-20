# Fast track simulation extension for the Acts project

This project provides the code to run a fast track simulation on top of 
`ActsCore`. The fast track simulation (Fatras) uses the actor plug-in 
mechanism of the `Acts::Propagator` and its predictive navigation through
the `Acts::TrackingGeometry` to simulate particle trajectories through the
tracking detector.


As a fast track simulation, it uses the reconstruction geometry description,
i.e. the `Acts::TrackingGeometry` as a simulation geometry with simplified and
particle parameterised material effects:
  * Multiple coulomb scattering is done by gaussian (mixture) approximations
  * Ionisation loss is calculated using the Bethe-Bloch formalism
  * Radiation loss follows Bethe-Heitler formalism
  * Limited nuclear interaction processes are parameterised from `Geant4`


Dependencies for the Core components are:
  * `ActsCore` and consequently `Eigen` and `Boost`

Optional dependency exists for:
  * `Geant4` for optional functionality taken from the full simulation toolkit
    
