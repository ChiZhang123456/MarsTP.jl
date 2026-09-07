# Historical compatibility notice. Never recompute Boris velocities for saved-state PSD.
error("This step is retired. Run analyze_saved_probes.jl RUN NEW_OUTPUT; it uses the same piecewise-linear saved velocities as MarsTP.forward_psd, with no MHD loading or reintegration.")
