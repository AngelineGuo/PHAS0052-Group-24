import os
import time
import requests
import re

# Define directories and get fasta files
base_dir = os.path.abspath(os.path.dirname(__file__))  # Directory containing this script
fasta_files_dir = os.path.join(base_dir, "..", "Files", "fasta")  # Directory containing the fasta files
fasta_files = []
for f in os.listdir(fasta_files_dir):
    fasta_files.append(os.path.join(fasta_files_dir, f))

# URL for form submission (action URL from the form)
url = "http://biomine.cs.vcu.edu/servers/biomine.php?name=NsitePred"

# Define your email address
your_email = "robin.carey17+binding@gmail.com"  # Replace this with your actual email address

# Function to submit the form
def submit_fasta(file_path, email):
    with open(file_path, 'r') as f:
        fasta_content = f.read()

    # Form data to simulate form submission
    form_data = {
        "seq": fasta_content,  # Input field for sequence
        "email1": email,       # Input field for email
    }

    # Send the POST request
    response = requests.post(url, data=form_data, allow_redirects=True)
    
    # Save the response to a file for debugging
    # with open(file_path.replace(".fasta", "_response.html"), "w") as f:
    #     f.write(response.text)

    # Check the response
    if response.status_code == 200:
        # Extract the processing URL from the response
        match = re.search(r'http://biomine.cs.vcu.edu/servers/biomine.php\?id=\d+', response.text)
        if match:
            processing_url = match.group(0)  # Get the matched URL
            # Follow the redirect URL to go to the processing page
            processing_response = requests.get(processing_url)
            if processing_response.status_code == 200:
                while "Your results are ready!" not in processing_response.text:
                    time.sleep(30) # Check results every 30 seconds
                    print("Waiting for results...")
                    processing_response = requests.get(processing_url)
                # Save the processing response to a file for debugging
                print("Results are ready!")
                # with open(file_path.replace(".fasta", "_processing_response.html"), "w") as f:
                #     f.write(processing_response.text)
                    
                # Extract the URL of the CSV file containing the results
                match = re.search(r'http://biomine.cs.vcu.edu/webresults/NsitePred/\d+/results.csv', processing_response.text)
                if match:
                    result_csv_url = match.group(0)
                    # Download the CSV file
                    csv_response = requests.get(result_csv_url)
                    if csv_response.status_code == 200:
                        csv_file_path = file_path.replace(".fasta", "_results.csv").replace("fasta", "Nsite_predictions")
                        with open(csv_file_path, "wb") as f:
                            f.write(csv_response.content)
                        print(f"CSV file saved for {file_path}")
                    else:
                        print(f"Failed to download CSV file for {file_path}. Status code: {csv_response.status_code}")
        elif "Your submission is now being processed" in response.text:
            print("Submission is being processed.")
        else:
            print(f"Submission successful for {file_path}")
    elif response.status_code in [301, 302]:
        redirect_url = response.headers.get("Location")
        print(f"Submission successful for {file_path}. Redirected to: {redirect_url}")
    else:
        print(f"Failed to process {file_path}. Status code: {response.status_code}")
        print(f"Response headers: {response.headers}")
        print(f"Response content: {response.text[:500]}")  # Print the first 500 characters for debugging

# Process each FASTA file
for fasta_file in fasta_files:
    submit_fasta(fasta_file, your_email)