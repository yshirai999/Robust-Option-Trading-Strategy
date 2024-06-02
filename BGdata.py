import pandas as pd
 
file_path = 'bgset512.dat'
 
try:
    # Assuming data is in a structured format like CSV or similar
    df = pd.read_csv(file_path, delimiter='\t')
    # Process the DataFrame
    print(df.info)
 
except FileNotFoundError:
    print(f"File '{file_path}' not found.") 
except Exception as e:
    print(f"An error occurred: {e}")