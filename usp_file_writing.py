# usp file writing
# What kind of data do we want to write???
import os


def list_to_tab_strings(vector):
    # function for taking an array of elements and turning them into strings separated by tabs

    new_string = ""
    for i in range(len(vector)):
        new_string = new_string + str(vector[i]) + '\t'
    new_string = new_string + '\n'
    return new_string


def to_excel_named(matrix, name):
    # toExcel where name can be inserted to auto write file without user input
    directory = os.getcwd()
    output_data_name = directory + "/" + name + ".csv"

    import csv

    try:
        with open(output_data_name,'w',newline='') as f:
            writer = csv.writer(f)
            writer.writerows(matrix)
    except:
        print('Error: Failed to write excel file.'
              ' Check if file is already open.')