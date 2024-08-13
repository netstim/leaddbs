import json
import sys

def BndBoxDimensionsEdits(json_file_path):
    # load the JSON file
    with open(json_file_path, 'r') as file:
            data = json.load(file)
    
    # increase the dimensions
    data['BrainRegion']['Dimension']['x[mm]'] = 70
    data['BrainRegion']['Dimension']['y[mm]'] = 70
    data['BrainRegion']['Dimension']['z[mm]'] = 70

    # save the file
    with open(json_file_path, 'w') as file:
            json.dump(data, file, indent=4)

if __name__ == "__main__":
    # get the JSON file path from the command line
    if len(sys.argv) != 2:
        print("Usage: python BndBoxDimensionsEdits.py <path_to_json_file>")
        sys.exit(1)
    
    json_file_path = sys.argv[1]
    BndBoxDimensionsEdits(json_file_path)