import matlab.engine
import os
# TODO: delete file
directory_path = os.path.join('C:\\', 'Users', '13365', 'Desktop', 'optimization_output', 'case_1')

eng = matlab.engine.start_matlab()
objective = eng.power_calculation(directory_path, 0.005, 0.274, 0.2108, 0.323, 1000, 1.0, 10.0)
print(objective)
#eng.power_calculation(nargout=0)
eng.quit()
