# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 19:18:12 2024

@author: 20213628
"""
import seaborn as sns
import matplotlib.pyplot as plt

def plot_density_for_7_lists(list_1, list_2, list_3, list_4, list_5, list_6, list_7):
    # Combine all lists into one for calculating min and max values
    all_lists = [list_1, list_2, list_3, list_4, list_5, list_6, list_7]
    
    # Calculate the minimum and maximum values of molecular weights across all lists
    min_value = min(min(lst) for lst in all_lists)
    max_value = max(max(lst) for lst in all_lists)
    
    labels = [1000, 2500, 5000, 7500, 10000, 20000, 30000]
    
    # Create a density plot for each list
    plt.figure(figsize=(10, 6))
    for i, lst in enumerate([list_1, list_2, list_3, list_4, list_5, list_6, list_7]):
        sns.kdeplot(lst, label=f'{labels[i]} ({len(lst)} Valid)', fill=False, common_norm=False)

    # Add labels to the axes
    plt.xlabel('Mol Weight')
    plt.ylabel('Density')

    # Add a legend
    plt.legend(title='Trainingset Size', loc='upper right')
    
    # Set x-limits based on calculated values
    plt.xlim(min_value, 200)

    # Show the plot
    plt.show()

# Example usage
list_1 = [50, 60, 70, 80, 90]
list_2 = [55, 65, 75, 85, 95]
list_3 = [40, 50, 60, 70, 80]
list_4 = [30, 40, 50, 60, 70]
list_5 = [70, 80, 90, 100, 110]
list_6 = [45, 55, 65, 75, 85]
list_7 = [60, 70, 80, 90, 100]

plot_density_for_7_lists(list_1, list_2, list_3, list_4, list_5, list_6, list_7)

