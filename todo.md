- DiscreteFunction problem
    - Need a way of loading data at each timestep into Sundance.
- Find original data paths for other files
    - Put them in the README
    - Maybe include actually-used data in the repo and assume that structure (rather than the config spaghetti)
- Look at getting per-county weather data (at least temperature) to see if it correlates with the apparent wave dynamics

## The Strange Counties
fix in `data_to_compartments`
- Figure out what is going on with Nebraska, Missouri, Nevada
    - Missouri shock flare around 440 days
    - Nebraska underreporting ~550-650 followed by sudden correction (plus another similar event around 900-950 days)
    - See I/N plot to make more apparent
- Nevada:
    - Eureka County has an overcorrection (goes from ~300 to ~150 at a certain point, which reflects as 0 cases for most of the rest of the time)
    - Humboldt County also has an overcorrection (goes from ~4500 to ~2500)
    - Potentially some other counties to lesser extent (all east border counties later on seem to have some dip)
- Also VT spike ~790 days 
- Also: "Unassigned" for each state - look at these, perhaps distribute (weighted by total population) 
