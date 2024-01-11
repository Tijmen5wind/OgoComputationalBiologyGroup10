# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 11:02:39 2024

@author: 20213628
"""
import matplotlib.pyplot as plt

def scatter_plot(list_x, list_y):
    # Check if the lengths of the lists are equal
    if len(list_x) != len(list_y):
        raise ValueError("The lengths of the lists must be equal.")

    # Create a scatter plot
    plt.scatter(list_x, [y * 100 for y in list_y], marker='o', color='blue')

    # Add labels to the axes
    plt.xlabel('Amount of Molecules')
    plt.ylabel('Uniqueness (%)')

    # Remove the grid
    plt.grid(False)

    # Adjust y-ticks to display exact percentages
    plt.yticks([y * 100 for y in list_y])

    # Adjust x-ticks
    plt.xticks(list_x, [f'{x // 1000}k' for x in list_x])

    # Show the plot
    plt.ylim(0, 102)
    plt.show()

# Example usage
example_list_x = [1000, 2500, 5000, 7500, 10000, 20000, 30000]
example_list_y = [8/9, 3/3, 8/17, 15/21, 26/47, 8/12, 10/16]

scatter_plot(example_list_x, example_list_y)


            
        