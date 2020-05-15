def parse_input(input_file_location):
    """Opens and reads in specified parameters in the user created .yaml file

    Args:
        input_file_location (str):

    Returns
        str
    """
    import yaml
    with open(input_file_location, 'r') as f:
        file_names = yaml.safe_load(f)
    return file_names
