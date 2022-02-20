.. include:: ../isonum.txt
.. include:: ../isogrk1.txt
.. |Ca2+| replace:: Ca\ :sup:`2+`

===========================
Visualization 1: graph plot
===========================

Either single run or ConnectRun stores the results of simulation in the LM-format files. The LM-format files are indeed HDF container files, and users can directly analyze them; however, it is not so easy to capture `the data structure <https://github.com/Luthey-Schulten-Lab/Lattice_Microbes/blob/master/docs/HDF5FileFormat.text>`_. Also, similar demands on the analyses are shared by users.

LD thus provides functions for analyses, in particular, to connect data across multiple unit runs. The total numbers of molecules are saved in the LM files; therefore, it is easy to obtain the total concentration of target molecules. In the script '51_plot_conc.py', the 'ConnectTotalConcs' class handle this (Lines 25-27). The variable 'domain_name' contains a domain name for the calculation of the volume, and the variable 'lm_files' contains a list of target LM files. In this case, Line 10 specifies a target directory ('simulation_dir'), in which all LM files are obtained in an ascending order of numbers (Line 22). The instance variable 't.timepoints' that shows observed timepoints, and the method 't.get_concs(species)' shows the time development of concentration of the specified species (Lines 49-51).

Some users may want to obtain the concentration of molecules in the labeled regions. This is enabled by two classes: 'GetLabeledConcs' and 'ConnectLabeledConcs'. The GetLabeledConcs class calculates the numbers of all molecules in all labeled volumes (Lines 29-35). The calculated numbers are saved in the h5 files ('conc_files'; Lines 23, 32, 35). The saved h5 files are loaded by the GetLabeledConcs class, and converted into a single time series. Users can plot the time development of concentration of target molecules within target regions (timepoints and get_concs; Lines 45-47).


.. literalinclude:: ../../tutorial/1/51_plot_conc.py
   :language: python
   :linenos:
   :caption: 51_plot_conc.py

|

In the case of FRAP (Lines 7-11), we can see the fluorescence recovery after photobleaching (FRAP) of YFP in the target spine.

.. image:: imgs/profile_photobleach.png
   :scale: 50%
   :align: center

|

In the case of |Ca2+| influx via NMDA receptors (Lines 14-18), we can see the transient increase in |Ca2+| concentration in the target spine.

.. image:: imgs/profile_Ca.png
   :scale: 50%
   :align: center


