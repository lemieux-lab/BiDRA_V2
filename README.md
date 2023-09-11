# BiDRA V2
BiDRA web interface version 2.0

The current version is running at https://bidrav2.bioinfo.iric.ca/ and correspond to the code in `BiDRA_dashboard_22/`.
To launch a server locally, just run 'bin/server'.

BiDRA V2 is implemented in Julia (V1.8.3). The Bayesian model and the inference process are defined and runned with Turing.jl (V0.26.0) and MCMCChains.jl (V6.0.3). Figures are generated with CairoMakie.jl (V0.10.6).
The web interface was made with Genie.jl (V5.18.1) and HTTP.jl (V1.9.6) on the server side, and Stipple.jl and StippleUI.jl for the interactivity.
