# ROOT-like analysis using pandas

This repository shows basics example of ROOT tree analysis using standard python tools such as pandas. The TTree considered is 
not totally flat and contains vectors of floats, which can have different size from entry to entry.

### Tool comparison
Three tools were compared in term of speed to load the dataset: [root_panas](https://github.com/scikit-hep/root_pandas), [root_numpy](https://github.com/scikit-hep/root_numpy) and [uproot](https://github.com/scikit-hep/uproot). 

### Event-by-event calculations

Using `dataFrame.apply()` function, observable at the event level were calculated using user-defined objects like Jet, Muon and Electrons. The benchmark used to compare tools is the minimum delta R between two arbitrary collections of objects in the event( e.g. min dR(mu,jets)). Below, there is an example of min dR(j,j) as a function of the number of jets in the event:

![mindR_vs_Njets](https://github.com/rmadar/root-with-pandas/blob/master/basic-computation/mindR_vs_Njets.png)
