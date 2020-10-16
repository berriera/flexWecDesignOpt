def read_output():
    """

    :return:
    """
    return


def convert_array_to_dict(array_name, array_values, starting_index=0):
    import numpy as np
    array_dictionary = {}
    for i in range(np.size(array_values)):
        element_value = array_values[i]
        key_name = array_name + '({})'.format(i + starting_index)
        array_dictionary[key_name] = element_value
    return array_dictionary


def write_dict_to_text_file(dictionary, file_name='design_variables.txt'):
    """"""  # TODO: add function description here
    f = open(file_name, 'w')
    for key in dictionary.keys():
        f.write('{}={}\n'.format(key, dictionary[key]))
    f.close()
