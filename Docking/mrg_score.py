###여러 디렉토리에 있는 CSV 파일을 읽고 데이터를 하나의 데이터프레임으로 결합하여 새 CSV 파일로 저장
import os
import pandas as pd

# Initialize an empty dataframe to store the combined data
combined_df = pd.DataFrame()

# Function to read and append the CSV file from each directory
def append_scores_from_directory(directory):
    try:
        # Build the file path
        file_path = os.path.join(directory, 'sorted_combined_autodock_scores.csv')
        
        # Check if the file exists
        if os.path.exists(file_path):
            # Read the CSV file
            df = pd.read_csv(file_path)
            df['Directory'] = directory  # Add a column to track which directory the data came from
            return df
        else:
            print(f"File not found in directory: {directory}")
            return None
    except Exception as e:
        print(f"Error processing {directory}: {e}")
        return None

# Load the list of directories from the dir_list file
with open('dir_list', 'r') as f:
    directories = f.read().splitlines()

# Process each directory and append the data
for directory in directories:
    if os.path.exists(directory):
        df = append_scores_from_directory(directory)
        if df is not None:
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    else:
        print(f"Directory {directory} does not exist")

# Save the combined dataframe to a new CSV file
output_path = 'combined_all_scores.csv'
combined_df.to_csv(output_path, index=False)

print(f"All scores have been combined and saved to {output_path}")
