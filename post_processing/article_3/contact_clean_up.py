import os
import pandas as pd

def delete_rows_sixth_column_zero(directory):
    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Check if the file is a .dou file
        if filename.endswith('.dou') and os.path.isfile(file_path):
            try:
                # Read the .dou file into a pandas DataFrame
                df = pd.read_csv(file_path)

                # Remove rows where the sixth column value is 0
                # Assuming the first column has no header and is unnamed or has a default header like 'Unnamed: 0'
                df = df[df.iloc[:, 6] != 0]

                # Write the modified DataFrame back to the CSV file
                df.to_csv(file_path, index=False)

                print(f"Processed file: {filename}")

            except Exception as e:
                print(f"An error occurred while processing {filename}: {e}")


if __name__=='__main__':


    file_dir = 'c:/Users/Axel/documents/DEM/Bertil_results/article_3/final_runs/1/swelling_electrode_calendering/contacts/'
    delete_rows_sixth_column_zero(file_dir)