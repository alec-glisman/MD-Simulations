# Plots the Total Energy vs Time of the Simulation
# ARG1 - filename

set title "Total Energy vs Iteration Number" font ",20"
set xlabel "Iteration Number"                                                       # Set x-axis label
set ylabel "Energy"                                                                 # Set y-axis label
set grid                                                                            # Add grid to the plot
set term png size 1000, 800                                                         # Save file as png (1000px x 800px)
unset key                                                                           # no key
set output ARG1                                                                # Specify output format
set datafile separator ','                                                          # Used for CSV file input
plot '../test_folder/test_file.csv' using 1:2 with linespoints notitle              # Plot command

# Replot and show output in X11
set term X11
set output
replot