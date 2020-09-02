import numpy as np
from flexWecDesignOpt.substitution import substitute_variables_in_line


def test_substitute__0_var_array():
    line_text = ''
    variables = np.asarray([0])
    obs = substitute_variables_in_line(line_text, variables)
    exp = ''
    assert obs == exp


def test_substitute__1_var_array():
    line_text = '?1?'
    variables = np.asarray([0])
    obs = substitute_variables_in_line(line_text, variables)
    exp = '0'
    assert obs == exp


def test_substitute__5_var_array():
    line_text = '?1? test ?2? test ?3? test ?4? test ?5? ?4? ?3? ?2? ?1?'
    variables = np.asarray([1, 2, 3, 4, 5])
    obs = substitute_variables_in_line(line_text, variables)
    exp = '1 test 2 test 3 test 4 test 5 4 3 2 1'
    assert obs == exp


def test_substitute__0_var_dict():
    line_text = ''
    variables = {'var1': 1}
    obs = substitute_variables_in_line(line_text, variables)
    exp = ''
    assert obs == exp


def test_substitute__1_var_dict():
    line_text = '?var1?'
    variables = {'var1': 1}
    obs = substitute_variables_in_line(line_text, variables)
    exp = '1'
    assert obs == exp


def test_substitute__5_var_dict():
    line_text = '?var1? test ?var2? test ?var3? test ?var4? test ?var5? ?var4? ?var3? ?var2? ?var1?'
    variables = {'var1': 1, 'var2': 2, 'var3': 3, 'var4': 4, 'var5': 5}
    obs = substitute_variables_in_line(line_text, variables)
    exp = '1 test 2 test 3 test 4 test 5 4 3 2 1'
    assert obs == exp


#def test_substitute__string_value_error(): # TODO : detect error catching here
#    line_text = '?1?'
#    variables = 'string_insert'
#    obs = substitute_variables_in_line(line_text, variables)
#    exp = ValueError
#    assert obs == exp

def test_substitute__1_row_array():
    line_text = '?var_array?'
    variables = {'var_array': np.asarray([0, 1, 2])}
    obs = substitute_variables_in_line(line_text, variables)
    exp = '0 1 2'
    assert obs == exp

def test_substitute__2_row_matrix():
    line_text = '?var_matrix?'
    variables = {'var_matrix': np.asarray([[0, 1], [2, 3]])}
    obs = substitute_variables_in_line(line_text, variables)
    exp = '0 1\n2 3'
    assert obs == exp

def test_substitute__3_row_matrix():
    line_text = '?var_matrix?'
    variables = {'var_matrix': np.asarray([[0, 1, 2], [3, 4, 5], [6, 7, 8]])}
    obs = substitute_variables_in_line(line_text, variables)
    exp = '0 1 2\n3 4 5\n6 7 8'
    assert obs == exp
