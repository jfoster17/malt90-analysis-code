# Malt90 Pipeline 

This project hosts the data reduction and analysis pipeline for the MALT90 project. The reduction pipeline is written in Python. It calls Livedata and Gridzilla to do the basic reduction and gridding of on-the-fly maps, and then uses custom Python routines to make moment maps. The pipeline reads and writes to a reduction log to keep track of the reduction status of data files. The analysis pipeline is a modification of the HOPS (The H2O southern Galactic Plane Survey) pipeline modified to fit multiple overlapping gaussian profiles.

