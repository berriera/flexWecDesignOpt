def read_output():
    """

    :return:
    """
    return


def write_dict_to_text_file(dictionary, file_name='design_variables.txt'):
    f = open(file_name, 'w')
    for key in dictionary.keys():
        f.write('{}={}\n'.format(key, dictionary[key]))
    f.close()
