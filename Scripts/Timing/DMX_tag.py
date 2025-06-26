import re
import os


def extract_values(par_file, mjd_file):
    dmxr1_values = {}
    dmxr2_values = {}
    dmxr1_pattern = re.compile(r'DMXR1_(\d+)')
    dmxr2_pattern = re.compile(r'DMXR2_(\d+)')

    with open(par_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        parts = line.split()
        if dmxr1_pattern.match(parts[0]):
            number = int(dmxr1_pattern.match(parts[0]).group(1))
            dmxr1_values[float(parts[1])] = (parts[0], number)
        elif dmxr2_pattern.match(parts[0]):
            number = int(dmxr2_pattern.match(parts[0]).group(1))
            dmxr2_values[float(parts[1])] = (parts[0], number)

    with open(mjd_file, 'r') as file:
        mjd_list = [float(line.strip()) for line in file.readlines()]

    mjd_list.sort()
    extracted_values = {}

    for mjd in mjd_list:
        closest_mjd1 = min(dmxr1_values.keys(), key=lambda x: abs(x - mjd))
        closest_mjd2 = min(dmxr2_values.keys(), key=lambda x: abs(x - mjd))
        extracted_values[mjd] = (dmxr1_values[closest_mjd1][0], dmxr1_values[closest_mjd1][1],
                                 closest_mjd1, dmxr2_values[closest_mjd2][0], dmxr2_values[closest_mjd2][1],
                                 closest_mjd2)

    return extracted_values


# User input for file paths
par_file_path = input("Enter the path to the par file: ")
mjd_file_path = input("Enter the path to the MJD file: ")

with open(par_file_path, 'r') as file:
	lines = file.readlines()
for line in lines:
	if line.strip().split()[0] == 'PSRJ' :
		psrname = line.strip().split()[1]

updated_par_file_path = "./"+psrname+".DMX.flag.par"

result = extract_values(par_file_path, mjd_file_path)

stored_values = []

print("DMXR1 and DMXR2 values with corresponding MJD")
print("--------------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------------")
for mjd, values in result.items():
    print(
        f"MJD: {mjd}, DMXR1: {values[0]}, MJD: {values[2]}, DMXR2: {values[3]}, MJD: {values[5]}")
    print("--------------------------------------------------------------------------------------")
    stored_values.append(values[1])  # Appending the value from index 1
    stored_values.append(values[4])  # Appending the value from index 4
    stored_values = list(set(stored_values))  # only give unique value
print("--------------------------------------------------------------------------------------")
print("Identified DMX")
print("--------------------------------------------------------------------------------------")
#print("DMX_:", stored_values)
for number in stored_values:
    # Pad the number with leading zeros to make it four digits
    formatted_number = f"{number:04d}"
    result = f"DMX_{formatted_number}"
    print(result)

# Extract DMX values from result strings
dmx_values = [f"DMX_{number:04d}" for number in stored_values]

# Open the par file for reading
with open(par_file_path, "r") as file:
    lines = file.readlines()

# Create the output directory if it doesn't exist
output_dir = os.path.dirname(updated_par_file_path)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Open the updated par file for writing
with open(updated_par_file_path, "w") as updated_file:
    for line in lines:
        # Check if the line starts with any of the DMX values
        if any(line.startswith(dmx) for dmx in dmx_values):
            # Split the line into parts
            parts = line.split()
            # Increase the size of the column to 4 and add "1" in the third place
            updated_line = f"{parts[0]} {parts[1]} 1 {' '.join(parts[2:])}\n"
            updated_file.write(updated_line)
        else:
            updated_file.write(line)

print("Your parfile is Updated Successfully, now you are ready to go:",
      updated_par_file_path)
